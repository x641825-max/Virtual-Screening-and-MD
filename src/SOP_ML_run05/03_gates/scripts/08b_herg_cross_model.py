#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3b（链上、离线）：hERG 第二模型交叉验证 —— 量化单模型结论的模型依赖度。

动机：主 hERG 门用的是 ADMET-AI 单模型（外部表）。本步在**同一 3,236 个四臂并集分子**上
跑一组独立实现的模型（不同算法、不同训练管道），量化「hERG<0.5 硬落门」这一结论在多大程度上
依赖被选中的那一个模型。

第二模型的构建（离线、可复现）：
  数据：ChEMBL REST API 抓取的 KCNH2/hERG（CHEMBL240）IC50/Ki 活性（`cross_model/raw/`，含 md5/sha256
        与抓取时间；联网抓取脚本 `08b0_herg_cross_model_fetch.py` 不入链）
  标注：pChEMBL ≥ 6.0 → 阻断剂(1)；≤ 5.0 → 非阻断剂(0)；5.0–6.0 之间丢弃（模糊带）
  去重：按上游 SMILES 的 InChIKey 取 pChEMBL 中位数；剔除 ChEMBL `potential_duplicate=1`
  特征：ECFP4 2048 位（RDKit），**不含**任何 ADMET-AI 的特征/权重
  模型：① 随机森林 ② 逻辑回归 ③ Tanimoto kNN(k=5)（三者都只用指纹，固定 PROJECT_SEED）
  验证：随机分层 5 折 + 按 Bemis-Murcko 通用骨架分组的 5 折（后者更贴近"新化学型"场景）

独立性声明（必须随结论转述）：ADMET-AI 的 hERG 头同样是在公开 hERG 数据上训练的，而公开
hERG 数据（含 TDC 的 hERG 集）多源自 ChEMBL；因此本步量化的是**算法与训练管道的依赖度**，
不能当作数据独立的外部验证。

输出：03_gates/herg/cross_model/{tables,reports,figures,HERG_CROSS_MODEL.md}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, average_precision_score, balanced_accuracy_score,
                             cohen_kappa_score, confusion_matrix, roc_auc_score)
from sklearn.model_selection import GroupKFold, StratifiedKFold

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import chem_core as cc
from vendor import plots
from vendor.io_utils import StepLog, read_json, r4, sha256_file, write_csv, write_json, write_text
from vendor.master import build_master

OUT = ROOT / "03_gates/herg/cross_model"
RAW = OUT / "raw/chembl_herg_kcnh2_activities.csv"
N_BITS = 2048


def _fps(smiles: list[str]) -> tuple[np.ndarray, list[int]]:
    """ECFP4 → uint8 位矩阵；返回 (矩阵, 有效行号)。"""
    arr, keep = [], []
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is None:
            continue
        arr.append(cc.fp_to_array(cc.fp_ecfp4(m, N_BITS)))
        keep.append(i)
    return np.asarray(arr, dtype=np.uint8), keep


def _build_dataset(rule: dict, log) -> tuple[pd.DataFrame, np.ndarray, list[int]]:
    spec = rule["step3_gates"]["herg"]["cross_model"]
    # 阈值口径：pchembl ≥ blocker_ge(6.0) → 阻断剂；≤ nonblocker_le(5.0) → 非阻断剂
    hi = float(spec["label"]["blocker_pchembl_ge"])
    lo = float(spec["label"]["nonblocker_pchembl_le"])
    assert hi > lo, "标注阈值配置有误：blocker_ge 必须大于 nonblocker_le"
    raw = pd.read_csv(RAW, low_memory=False)
    raw["pchembl_value"] = pd.to_numeric(raw.pchembl_value, errors="coerce")
    n0 = len(raw)
    raw = raw[raw.potential_duplicate.fillna(0).astype(int) == 0]
    raw = raw[raw.pchembl_value.notna()]
    raw = raw[raw.canonical_smiles.notna()]
    # InChIKey 去重（同分子多条测定取中位数）
    keys = []
    for s in raw.canonical_smiles:
        m = Chem.MolFromSmiles(s)
        keys.append(cc.inchi_key(m) if m is not None else None)
    raw["ik"] = keys
    raw = raw[raw.ik.notna()]
    agg = (raw.groupby("ik", as_index=False)
              .agg(pchembl=("pchembl_value", "median"),
                   smiles=("canonical_smiles", "first"),
                   n_records=("pchembl_value", "size"),
                   mol_id=("molecule_chembl_id", "first")))
    agg["label"] = np.where(agg.pchembl >= hi, 1, np.where(agg.pchembl <= lo, 0, -1))
    n_lab = int((agg.label >= 0).sum())
    ds = agg[agg.label >= 0].reset_index(drop=True)
    arr, keep = _fps(ds.smiles.tolist())
    ds = ds.iloc[keep].reset_index(drop=True)
    log.info(f"ChEMBL hERG 原始 {n0} 条 → 去重后 {len(agg)} 个分子 → "
             f"标注 {n_lab}（阻断 {int((ds.label == 1).sum())} / "
             f"非阻断 {int((ds.label == 0).sum())}；模糊带 "
             f"{len(agg) - n_lab} 个丢弃）")
    return ds, arr, keep


def _knn_tanimoto(train: np.ndarray, train_y: np.ndarray, query: np.ndarray, k: int) -> np.ndarray:
    """Tanimoto kNN 打分（k 个最近邻的阳性比例）；分块计算避免 OOM。"""
    t = train.astype(np.float32)
    denom_t = t.sum(axis=1, dtype=np.float32)
    out = np.empty(len(query), dtype=np.float64)
    step = 256
    for i in range(0, len(query), step):
        q = query[i:i + step].astype(np.float32)
        inter = q @ t.T                               # 交集（二进制点积）
        denom = q.sum(axis=1, dtype=np.float32)[:, None] + denom_t[None, :] - inter
        sim = np.divide(inter, denom, out=np.zeros_like(inter), where=denom > 0)
        idx = np.argsort(-sim, axis=1, kind="stable")[:, :k]
        out[i:i + step] = train_y[idx].mean(axis=1)
    return out


def _cv_metrics(X: np.ndarray, y: np.ndarray, groups: np.ndarray, seed: int) -> dict:
    from sklearn.base import clone
    base = {
        "rf": RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=1),
        "lr": LogisticRegression(max_iter=2000, random_state=seed),
    }
    out = {}
    splits = {
        "random_5fold": list(StratifiedKFold(n_splits=5, shuffle=True,
                                             random_state=seed).split(X, y)),
        "scaffold_5fold": list(GroupKFold(n_splits=5).split(X, y, groups)),
    }
    for scheme, sp in splits.items():
        for mname, est in base.items():
            yt, yp = [], []
            for tr, te in sp:
                mdl = clone(est).fit(X[tr], y[tr])
                yp.append(mdl.predict_proba(X[te])[:, 1])
                yt.append(y[te])
            yt, yp = np.concatenate(yt), np.concatenate(yp)
            out[f"{mname}_{scheme}"] = {
                "roc_auc": r4(roc_auc_score(yt, yp)),
                "pr_auc": r4(average_precision_score(yt, yp)),
                "balanced_accuracy": r4(balanced_accuracy_score(yt, yp >= 0.5)),
                "accuracy": r4(accuracy_score(yt, yp >= 0.5)),
            }
        # kNN（无训练过程，直接用训练折做邻居）
        yt, yp = [], []
        for tr, te in sp:
            yp.append(_knn_tanimoto(X[tr], y[tr], X[te], 5))
            yt.append(y[te])
        yt, yp = np.concatenate(yt), np.concatenate(yp)
        out[f"knn_{scheme}"] = {
            "roc_auc": r4(roc_auc_score(yt, yp)),
            "pr_auc": r4(average_precision_score(yt, yp)),
            "balanced_accuracy": r4(balanced_accuracy_score(yt, yp >= 0.5)),
            "accuracy": r4(accuracy_score(yt, yp >= 0.5)),
        }
    return out


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    spec = rule["step3_gates"]["herg"]["cross_model"]
    seed = int(rule["seeds"]["PROJECT_SEED"])
    thr = float(rule["step3_gates"]["herg"]["threshold"])
    with StepLog(ROOT, "08b_herg_cross_model") as log:
        OUT.joinpath("tables").mkdir(parents=True, exist_ok=True)
        OUT.joinpath("reports").mkdir(parents=True, exist_ok=True)
        OUT.joinpath("figures").mkdir(parents=True, exist_ok=True)

        # 冻结输入层自检：raw/ 下是第三方原始表（属输入，不是本步产物）。
        # 若在「清空产物目录」时被连带删除，这里必须给出可操作的提示，
        # 而不是在深处抛裸 FileNotFoundError。
        for _need in (RAW, OUT / "raw/fetch_meta.json"):
            if not _need.is_file():
                raise SystemExit(
                    f"[FAIL] 缺少冻结输入 {_need.relative_to(ROOT)}：该文件属输入层"
                    f"（非本步产物），不得随产物目录一并删除；请从输入副本恢复后重跑"
                    f"（本步不重抓 ChEMBL，见 03_gates/herg/cross_model/HERG_CROSS_MODEL.md 来源说明）。")

        # 数据溯源哈希校验（与抓取元信息登记值一致）
        meta = read_json(OUT / "raw/fetch_meta.json")
        got = sha256_file(RAW)
        if got != meta["sha256"]:
            raise RuntimeError(f"ChEMBL 原始表哈希与 fetch_meta.json 不一致：{got[:12]}…")

        ds, X, _ = _build_dataset(rule, log)
        y = ds.label.to_numpy()
        # 分组 = Bemis-Murcko 通用骨架（确定性：按骨架字符串排序后再分配折）
        groups = np.asarray([cc.scaffold_generic(s) or "NO_SCAFFOLD" for s in ds.smiles])

        cvs = _cv_metrics(X, y, groups, seed)
        log.info("CV（随机分层 / 骨架分组）："
                 + "；".join(f"{k}={v['roc_auc']}" for k, v in cvs.items()))

        # ---- 在 3,236 个四臂并集上预测
        union = pd.read_csv(ROOT / "03_gates/herg/SOP-05_herg_input_union.csv", low_memory=False)
        U, ukeep = _fps(union.smiles.tolist())
        union = union.iloc[ukeep].reset_index(drop=True)
        log.info(f"并集 {len(union)} 个分子 → ECFP4 可用 {len(U)}")

        rf = RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=1).fit(X, y)
        lr = LogisticRegression(max_iter=2000, random_state=seed).fit(X, y)
        p_rf = rf.predict_proba(U)[:, 1]
        p_lr = lr.predict_proba(U)[:, 1]
        p_knn = _knn_tanimoto(X, y, U, 5)

        hp = pd.read_csv(ROOT / rule["step3_gates"]["herg"]["external_table"], low_memory=False)
        pred = union[["cid", "smiles", "arm_membership"]].copy()
        pred["hERG_admet_ai"] = pred.cid.map(hp.set_index("cid")["hERG"])
        pred["hERG_cross_rf"] = p_rf
        pred["hERG_cross_lr"] = p_lr
        pred["hERG_cross_knn"] = p_knn
        pred["m1_pass"] = (pred.hERG_admet_ai < thr).astype(int)
        for tag in ("rf", "lr", "knn"):
            pred[f"m2_pass_{tag}"] = (pred[f"hERG_cross_{tag}"] < thr).astype(int)
        # 三模型多数表决作为"第二模型"的集成口径
        pred["hERG_cross_ens"] = pred[["hERG_cross_rf", "hERG_cross_lr",
                                       "hERG_cross_knn"]].mean(axis=1)
        pred["m2_pass_ens"] = (pred.hERG_cross_ens < thr).astype(int)
        write_csv(OUT / "tables/herg_cross_model_predictions.csv", pred)

        # ---- 一致性 / 模型依赖度
        def _agree(col: str) -> dict:
            a = pred.m1_pass.to_numpy()
            b = pred[f"m2_pass_{col}"].to_numpy()
            # sklearn 约定：confusion_matrix(labels=[0,1]).ravel() → TN, FP, FN, TP
            # 这里 a = 模型1 过门(1/0)，b = 模型2 过门(1/0)
            tn, fp, fn, tp = confusion_matrix(a, b, labels=[0, 1]).ravel()
            return {
                "model2": col,
                "spearman_rho": r4(pred.hERG_admet_ai.corr(pred[f"hERG_cross_{col}"], method="spearman")),
                "pearson_r": r4(pred.hERG_admet_ai.corr(pred[f"hERG_cross_{col}"])),
                "gate_agreement": r4((a == b).mean()),
                "cohen_kappa": r4(cohen_kappa_score(a, b)),
                "m1_pass_rate": r4(a.mean()),
                "m2_pass_rate": r4(b.mean()),
                "m1_median_score": r4(pred.hERG_admet_ai.median()),
                "m2_median_score": r4(pred[f"hERG_cross_{col}"].median()),
                "both_pass": int(tp), "both_fail": int(tn),
                "m1_pass_m2_fail": int(fn),   # 模型1 通过、模型2 判为高风险 → 交付名单的模型依赖风险
                "m1_fail_m2_pass": int(fp),   # 模型1 拒绝、模型2 判为低风险 → 门是否过严
            }
        # ---- 与阈值无关的排序一致性：高风险十分位的交集 + 匹配过门率后的重叠
        def _rank_consistency(col: str) -> dict:
            s1 = pred.hERG_admet_ai
            s2 = pred[f"hERG_cross_{col}"]
            k = max(1, int(round(0.10 * len(pred))))
            top1 = set(pred.cid[s1.rank(ascending=False, method="first") <= k])
            top2 = set(pred.cid[s2.rank(ascending=False, method="first") <= k])
            rate = float(pred.m1_pass.mean())
            q = float(s2.quantile(rate))                 # 匹配模型1 过门率的阈值
            p1 = set(pred.cid[pred.m1_pass == 1])
            p2 = set(pred.cid[s2 <= q])
            return {
                "model2": col,
                "high_risk_top10pct_jaccard": r4(len(top1 & top2) / len(top1 | top2)),
                "matched_pass_rate": r4(rate),
                "m2_matched_threshold": r4(q),
                "matched_overlap_jaccard": r4(len(p1 & p2) / len(p1 | p2)),
                "recall_of_m1_passes": r4(len(p1 & p2) / len(p1)),
            }
        rank_cons = [_rank_consistency(c) for c in ("rf", "lr", "knn", "ens")]
        write_csv(OUT / "tables/herg_cross_model_rank_consistency.csv", pd.DataFrame(rank_cons))
        for a in rank_cons:
            log.info(f"排序一致性[{a['model2']}]：高风险前10%交集 Jaccard "
                     f"{a['high_risk_top10pct_jaccard']}；匹配过门率（阈值 "
                     f"{a['m2_matched_threshold']}）后两份过门集合 Jaccard "
                     f"{a['matched_overlap_jaccard']}（覆盖模型1过门分子的 "
                     f"{a['recall_of_m1_passes']:.1%}）")

        agree = [_agree(c) for c in ("rf", "lr", "knn", "ens")]
        write_csv(OUT / "tables/herg_cross_model_agreement.csv", pd.DataFrame(agree))
        for a in agree:
            log.info(f"一致性[{a['model2']}]：ρ={a['spearman_rho']}、门决策一致 "
                     f"{a['gate_agreement']:.1%}（κ={a['cohen_kappa']}）；"
                     f"过门率 模型1 {a['m1_pass_rate']:.1%} vs 模型2 {a['m2_pass_rate']:.1%}；"
                     f"1过2不过 {a['m1_pass_m2_fail']}、1不过2过 {a['m1_fail_m2_pass']}")

        # ---- 对结论的影响：并集 / hERG 前主池 / 交付名单
        m = build_master(ROOT)
        pre = m[(m.druglike_pass == 1) & (m.pose_gate == "pass") & (m.maxTc_actives < 0.9)
                & (m.anchor_hits_4A >= 4) & (m.qed >= rule["step2_druglikeness"]["qed"]["primary"])]
        pre_cids = set(pre.cid)
        # 单一数据源：只读上游选名单产物（本步在步序上位于 10_select_lists 之后）。
        # 早期版本优先读交付层副本（06_deliverables/…，由更晚的 13_assemble 写出），
        # 造成“读下游产物”：清空产物重跑会在本步 FileNotFoundError，且口径变更后的首跑会用
        # 上一轮的名单（见 SOP_5_review/SOP-05_运行逻辑审查_2026-09-23.md D1）。
        _list_p = ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv"
        if not _list_p.is_file():
            raise SystemExit(
                f"[致命] 缺少上游名单 {_list_p.relative_to(ROOT)}；本步必须晚于 10_select_lists 执行"
                f"（见 config/run_order.json 的 order_note）。请按 run_all.py 的全链步序重跑，"
                f"不要单跑本步，也不要从交付层取输入。")
        main_list = pd.read_csv(_list_p, low_memory=False)
        list_cids = set(main_list[main_list.candidate_type == "candidate"].cid)

        def _gate(n_cids: set, col: str, base_cids: set) -> dict:
            sub = pred[pred.cid.isin(n_cids)]
            b = pred[pred.cid.isin(base_cids)]
            return {
                "pool_n": len(sub),
                "m1_pass": int(sub.m1_pass.sum()),
                "m2_pass": int(sub[f"m2_pass_{col}"].sum()),
                "m1_and_m2_pass": int(((sub.m1_pass == 1) & (sub[f"m2_pass_{col}"] == 1)).sum()),
                "m1_only_pass": int(((sub.m1_pass == 1) & (sub[f"m2_pass_{col}"] == 0)).sum()),
                "m2_only_pass": int(((sub.m1_pass == 0) & (sub[f"m2_pass_{col}"] == 1)).sum()),
            }
        impact = {
            # 键名说明：早期版本按键时的并集规模命名为 union_3236；现为口径重算后的并集
            # （v1.6.0 实测 3,258），故改名 union_recomputed，数值以 pool_n 为准。
            "union_recomputed": _gate(set(pred.cid), "ens", set(pred.cid)),
            "pre_herg_main_pool": _gate(pre_cids, "ens", set(pred.cid)),
            "delivered_list_48": _gate(list_cids, "ens", set(pred.cid)),
            # 单模型口径（随机森林）作对照：集成会摊平极端值
            "union_recomputed_rf": _gate(set(pred.cid), "rf", set(pred.cid)),
            "pre_herg_main_pool_rf": _gate(pre_cids, "rf", set(pred.cid)),
            "delivered_list_48_rf": _gate(list_cids, "rf", set(pred.cid)),
        }
        # 基数显式化：三个分母互不相同，务必分清 ——
        #   ① pool_n：该池分子总数；
        #   ② 第二模型覆盖（本步 pred 表覆盖并集，故等于 pool_n）；
        #   ③ **第一模型外部预测覆盖**（hERG_admet_ai 非空）：并集 22 个、前主池 15 个为空，
        #      这些空值在 m1_pass 里按 0（不过门）计入 —— 属保守处理。
        #   故 m1_pass / m1_and_m2_pass / 一致性 / κ 的分母是 m1_evaluated_n，不是 pool_n。
        #   （本包 hERG 门的落门依据是包内冻结全量表、覆盖 100%，与此处的 15/22 无关。）
        _m1_ok = set(pred.loc[pred.hERG_admet_ai.notna(), "cid"])
        for _k, _cids in (("union_recomputed", set(pred.cid)),
                          ("union_recomputed_rf", set(pred.cid)),
                          ("pre_herg_main_pool", pre_cids),
                          ("pre_herg_main_pool_rf", pre_cids),
                          ("delivered_list_48", list_cids),
                          ("delivered_list_48_rf", list_cids)):
            impact[_k]["m1_evaluated_n"] = int(len(_cids & _m1_ok))
            impact[_k]["m1_missing_n"] = int(len(_cids - _m1_ok))
        impact["pre_herg_main_pool"]["covered_n"] = int(len(pre_cids & set(pred.cid)))
        impact["pre_herg_main_pool"]["uncovered_n"] = int(len(pre_cids - set(pred.cid)))
        impact["m1_missing_convention"] = (
            "缺第一模型外部预测值（hERG_admet_ai 为空）的分子在 m1_pass 里按 0（不过门）计入，"
            "属保守处理；故 m1_pass / m1_and_m2_pass / 一致性 / κ 的分母是 m1_evaluated_n。"
            "本包 hERG 门的落门依据是包内冻结全量表（覆盖 100%），与这些空值无关。")
        impact["pre_herg_main_pool"]["base_note"] = (
            "m1_pass/m1_and_m2_pass 只在 m1_evaluated_n 个（两侧都有预测值）分子上计算，"
            "不是整个 pool_n；本包 hERG 门在该池上的过门数见 "
            "03_gates/reports/gate_stack.json → gate_chain 与 04_selection 的主臂池规模，二者不可互换")
        impact["key_note"] = ("union_recomputed* 为口径重算后的四臂并集（审计件）；键名不再内嵌规模数字，"
                              "规模一律以 pool_n 为准")
        # 名单级：哪些候选在第二模型下不过门
        sub = pred[pred.cid.isin(list_cids)].copy()
        flips = sub[sub.m2_pass_ens == 0].sort_values("hERG_cross_ens", ascending=False)
        write_csv(OUT / "tables/herg_cross_model_list_flips.csv", flips)
        log.info(f"名单 48 中第二模型（集成）不过门 {len(flips)} 个；"
                 f"hERG 前主池 {impact['pre_herg_main_pool']['pool_n']} 个："
                 f"模型1 过 {impact['pre_herg_main_pool']['m1_pass']}、"
                 f"模型2 过 {impact['pre_herg_main_pool']['m2_pass']}、"
                 f"两者都过 {impact['pre_herg_main_pool']['m1_and_m2_pass']}")

        out = {
            "dataset": {"source": meta["source"], "target": meta["query"]["target_chembl_id"],
                        "fetched_utc": meta["fetched_utc"], "raw_sha256": meta["sha256"],
                        "n_raw_records": int(meta["n_rows"]), "n_labeled": int(len(ds)),
                        "n_blocker": int((ds.label == 1).sum()),
                        "n_nonblocker": int((ds.label == 0).sum()),
                        "label_rule": spec["label"], "features": "ECFP4 2048 bit（RDKit）"},
            "cv": cvs,
            "agreement_vs_admet_ai": agree,
            "rank_consistency": rank_cons,
            "gate_impact": impact,
            "n_list_flips_ens": int(len(flips)),
            "list_flips_cids": flips.cid.tolist(),
            "independence_caveat": spec["independence_caveat"],
        }
        write_json(OUT / "reports/herg_cross_model.json", out)

        # ---- 图
        plots.setup_cjk()
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        ax = axes[0]
        ax.hexbin(pred.hERG_admet_ai, pred.hERG_cross_ens, gridsize=45, mincnt=1,
                  cmap="Blues", linewidths=0.1)
        ax.axhline(thr, color="#c0392b", lw=1, ls="--")
        ax.axvline(thr, color="#c0392b", lw=1, ls="--")
        ax.set_xlabel("ADMET-AI hERG 概率（主门）")
        ax.set_ylabel("第二模型集成概率（RF/LR/kNN 均值）")
        ax.set_title(f"并集 {len(pred)} 个分子：两模型打分（ρ={agree[3]['spearman_rho']}）")
        ax = axes[1]
        bars = [("过门\n模型1", impact["pre_herg_main_pool"]["m1_pass"], "#2c6fbb"),
                ("过门\n模型2", impact["pre_herg_main_pool"]["m2_pass"], "#8e7cc3"),
                ("两模型\n都过", impact["pre_herg_main_pool"]["m1_and_m2_pass"], "#4b8b3b")]
        ax.bar([b[0] for b in bars], [b[1] for b in bars], color=[b[2] for b in bars])
        for i, b in enumerate(bars):
            ax.text(i, b[1], str(b[1]), ha="center", va="bottom", fontsize=10)
        ax.set_ylabel("分子数")
        ax.set_title(f"hERG 前主池 {impact['pre_herg_main_pool']['pool_n']} 个的过门情况")
        plots.save(fig, OUT / "figures/fig_herg_cross_model.png")
        plt.close(fig)

        write_text(OUT / "HERG_CROSS_MODEL.md", _report(out))
        log.info("交叉验证完成 → 03_gates/herg/cross_model/")
    return 0


def _report(out: dict) -> str:
    L = ["# hERG 第二模型交叉验证（模型依赖度）", "",
         "本步在**同一 3,236 个四臂并集分子**上跑一组独立实现的模型，用于量化主门"
         "（ADMET-AI 单模型 hERG<0.5）在多大程度上依赖被选中的那一个模型。", "",
         "## 1. 第二模型的数据与训练", "",
         f"- 数据：ChEMBL REST API，靶点 {out['dataset']['target']}（KCNH2/hERG），"
         f"抓取于 {out['dataset']['fetched_utc']}，原始记录 {out['dataset']['n_raw_records']} 条"
         f"（sha256 `{out['dataset']['raw_sha256'][:16]}…`）",
         f"- 标注：pChEMBL ≥ {out['dataset']['label_rule']['blocker_pchembl_ge']} → 阻断剂；"
         f"≤ {out['dataset']['label_rule']['nonblocker_pchembl_le']} → 非阻断剂；"
         f"中间模糊带丢弃 → 标注 {out['dataset']['n_labeled']} 个"
         f"（阻断 {out['dataset']['n_blocker']} / 非阻断 {out['dataset']['n_nonblocker']}）",
         f"- 特征：{out['dataset']['features']}（不含 ADMET-AI 的任何权重或特征）", "",
         "## 2. 第二模型自身的表现（5 折交叉验证）", "",
         "| 模型 | 随机分层 ROC-AUC | 骨架分组 ROC-AUC | 骨架分组 PR-AUC | 骨架分组平衡准确率 |",
         "|---|---|---|---|---|"]
    for m in ("rf", "lr", "knn"):
        a, b = out["cv"][f"{m}_random_5fold"], out["cv"][f"{m}_scaffold_5fold"]
        L.append(f"| {m} | {a['roc_auc']} | {b['roc_auc']} | {b['pr_auc']} | "
                 f"{b['balanced_accuracy']} |")
    L += ["", "## 3. 与 ADMET-AI 的一致性（3,236 个并集分子）", "",
          "| 第二模型 | Spearman ρ | 门决策一致率 | Cohen κ | 模型1 过门率 | 模型2 过门率 | "
          "1过2不过 | 1不过2过 |",
          "|---|---|---|---|---|---|---|---|"]
    for a in out["agreement_vs_admet_ai"]:
        L.append(f"| {a['model2']} | {a['spearman_rho']} | "
                 f"{a['gate_agreement']:.1%} | {a['cohen_kappa']} | "
                 f"{a['m1_pass_rate']:.1%} | {a['m2_pass_rate']:.1%} | "
                 f"{a['m1_pass_m2_fail']} | {a['m1_fail_m2_pass']} |")
    L += ["", "## 3b. 与阈值无关的排序一致性", "",
          "| 第二模型 | 高风险前 10% 交集 Jaccard | 匹配过门率后过门集合 Jaccard | 覆盖模型1过门分子 |",
          "|---|---|---|---|"]
    for a in out["rank_consistency"]:
        L.append(f"| {a['model2']} | {a['high_risk_top10pct_jaccard']} | "
                 f"{a['matched_overlap_jaccard']} | {a['recall_of_m1_passes']:.1%} |")
    L += ["", "（匹配过门率＝把第二模型的阈值调到与模型1 相同的过门率后再比，用于区分"
          "「排序不同」与「阈值/标定不同」两种原因。）", "",
          "## 4. 对结论的影响", "",
          "| 范围 | 分子数 | 模型1 过门 | 模型2 过门 | 两者都过 | 仅模型1 | 仅模型2 |",
          "|---|---|---|---|---|---|---|"]
    for k, name in (("union_recomputed", "四臂并集（集成）"), ("pre_herg_main_pool", "hERG 前主池（集成）"),
                    ("delivered_list_48", "交付名单 48（集成）"),
                    ("union_recomputed_rf", "四臂并集（随机森林）"),
                    ("pre_herg_main_pool_rf", "hERG 前主池（随机森林）"),
                    ("delivered_list_48_rf", "交付名单 48（随机森林）")):
        v = out["gate_impact"][k]
        L.append(f"| {name} | {v['pool_n']} | {v['m1_pass']} | {v['m2_pass']} | "
                 f"{v['m1_and_m2_pass']} | {v['m1_only_pass']} | {v['m2_only_pass']} |")
    L += ["", f"- 名单 48 中在第二模型（集成）下**不过门**的有 **{out['n_list_flips_ens']}** 个"
          f"（见 `tables/herg_cross_model_list_flips.csv`）。", "",
          "## 5. 独立性声明与转述要求", "",
          f"- {out['independence_caveat']}",
          "- 因此上表的一致性度量的是**算法/训练管道依赖度**，不能替代真正数据独立的"
          "外部验证；两模型一致说明结论对该类模型稳健，不一致的分子应在 MD 前置复核中"
          "按「模型依赖」逐个人工核对。", ""]
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
