#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 4：6 组名单矩阵（3 QED 档 × 2 锚点档）+ 约束触发率 + 决策卡 + 2 阴性对照。

选序（frozen_rule.step4_selection）：CNN_VS 降序 + 成对 ECFP4 Tanimoto < 0.65 +
骨架封顶（主 ≤1；敏感 ≤3/≤5）+ 化学型 ≤ round(8N/20) + AM1-BCC 风险档 ≤ round(3N/20)；
名单尾部固定 2 个阴性对照（NEG1 = 诱饵中 CNN_VS 最高且锚点 4/4；
NEG2 = 诱饵中 Vina 最强且锚点 4/4），两条规则均在本包内从 run01 归档重算。
输出：04_selection/tables/list_*.csv、decision_cards_*.csv、constraint_trigger_rates.csv
      04_selection/reports/selection_summary.json、figures/fig_list_overlap.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import RDLogger
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import picker as pk
from vendor import plots
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json
from vendor.master import build_master

RDLogger.DisableLog("rdApp.*")
ARMS = [("main_QED0.50_anchor4", 0.50, 4), ("sens_QED0.60_anchor4", 0.60, 4),
        ("sens_QED0.40_anchor4", 0.40, 4), ("sens_QED0.50_anchor3", 0.50, 3),
        ("sens_QED0.60_anchor3", 0.60, 3), ("sens_QED0.40_anchor3", 0.40, 3)]


def load_fps(cids: list) -> dict:
    z = np.load(ROOT / "02_druglikeness/tables/fingerprints_ecfp4.npz", allow_pickle=True)
    out = {}
    for cid, bits in zip(z["cid"], z["bits"]):
        bv = ExplicitBitVect(2048)
        bv.SetBitsFromList(np.flatnonzero(bits).tolist())
        out[str(cid)] = bv
    missing = [c for c in cids if c not in out]
    if missing:
        raise RuntimeError(f"指纹缺失 {len(missing)} 条：{missing[:5]}")
    return out


def known_iks() -> set:
    p = ROOT / "03_gates/tables/known_ligand_similarity.csv"
    if not p.is_file():
        return set()
    kk = pd.read_csv(p)
    return set(kk.loc[kk.known_exact_hit == 1, "ik"])


def pick_negative_controls(log) -> pd.DataFrame:
    """从 run01 归档重算 NEG1/NEG2（不引用包外既有结果）。"""
    rank = pd.read_csv(ROOT / "00_inputs/extracted/run01/02_enrichment_5098/scores/gnina_ranking.csv")
    ifp = ROOT / "01_pose_qc/tables/labeled5098_contact_fp.csv"
    if not ifp.is_file():
        raise RuntimeError("缺少 5,098 集 IFP（需 extract_large=true 与步骤 04）")
    lab = pd.read_csv(ifp)
    d = rank.merge(lab[["mol_id", "anchor_hits_4A"]], on="mol_id", how="left")
    dec = d[(d.activity_label == 0) & (d.anchor_hits_4A == 4)]
    neg1 = dec.sort_values(["cnn_vs", "mol_id"], ascending=[False, True]).iloc[0]
    neg2 = dec.sort_values(["vina_affinity", "mol_id"], ascending=[True, True]).iloc[0]
    out = pd.DataFrame([
        {"cid": neg1.mol_id, "control_role": "NEG1",
         "rule": "诱饵中 CNN_VS 最高且锚点 4/4（高可信假阳性）",
         "cnn_vs": r4(neg1.cnn_vs), "vina_affinity": r4(neg1.vina_affinity),
         "smiles": neg1.smiles, "activity_label": 0,
         "anchor_hits_4A": int(neg1.anchor_hits_4A)},
        {"cid": neg2.mol_id, "control_role": "NEG2",
         "rule": "诱饵中 Vina 最强且锚点 4/4（两引擎共识诱饵）",
         "cnn_vs": r4(neg2.cnn_vs), "vina_affinity": r4(neg2.vina_affinity),
         "smiles": neg2.smiles, "activity_label": 0,
         "anchor_hits_4A": int(neg2.anchor_hits_4A)},
    ])
    log.info(f"阴性对照重算：NEG1={neg1.mol_id}（CNN_VS={neg1.cnn_vs:.4f}，诱饵数 {len(dec)}）、"
             f"NEG2={neg2.mol_id}（Vina={neg2.vina_affinity:.4f}）")
    return out


def card_rows(m: pd.DataFrame, picks_df: pd.DataFrame, neg: pd.DataFrame,
              herg_status: str) -> pd.DataFrame:
    want = ["cid", "smiles", "CNN_VS", "qed", "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB",
            "cns_mpo", "chemotype", "scaffold", "scaffold_generic", "anchor_hits_4A",
            "n_contacts_4A", "pose_gate", "maxTc_actives", "nearest_active_id", "herg_prob",
            "herg_gate", "herg_status", "risk_am1bcc", "n_aromN", "n_N", "brenk_names",
            "pains", "ro5_viol", "veber_viol", "qc_class", "stereo_source",
            "protonation_note", "basic_pka", "clogd74", "clogd74_lower_80",
            "clogd74_upper_80", "pka_status", "cns_mpo_lo80", "cns_mpo_hi80",
            "cns_mpo_pi80_span"]
    cols = [c for c in want if c in m.columns]
    cards = picks_df[["cid", "pick_order"]].merge(m[cols], on="cid", how="left")
    cards["candidate_type"] = "candidate"
    ctrl = neg.copy()
    ctrl["pick_order"] = range(len(picks_df) + 1, len(picks_df) + 1 + len(ctrl))
    ctrl["candidate_type"] = "negative_control"
    ctrl["herg_status"] = herg_status
    return pd.concat([cards, ctrl], ignore_index=True, sort=False).sort_values(
        "pick_order").reset_index(drop=True)


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    sel = rule["step4_selection"]
    n_total, n_ctrl = sel["list_size"], sel["n_negative_controls"]
    with StepLog(ROOT, "10_select_lists") as log:
        m = build_master(ROOT, log=log)
        herg_status = str(m.herg_status.iloc[0])
        fps = load_fps(m.cid.tolist())
        neg = pick_negative_controls(log)
        ki = known_iks()

        base = m[(m.druglike_pass == 1) & (m.pose_gate == "pass") & (m.maxTc_actives < 0.9)]
        if herg_status == "applied":
            # 前置拦截（与规则 step3_gates.herg.missing_policy 对齐）：门列是可空 Int64，
            # pandas 对 NA 的 `== 1` 返回 <NA>，布尔索引会**丢弃**该行 —— 那就等于把
            # not_evaluated 静默当成「不合格」，与规则「缺值不写作 fail / 不过门 / 高风险」相反。
            # 全池口径（scope=pool）下缺值数应为 0；若不为 0，说明门作用域没覆盖全池，
            # 此处必须响亮失败，而不是静默降级。
            _na = int(m.herg_gate.isna().sum())
            _scope = str(rule["step3_gates"]["herg"].get("scope", ""))
            if _scope == "pool" and _na:
                raise SystemExit(
                    f"[致命] hERG 门列含 {_na} 个缺值（not_evaluated），而规则要求 scope=pool 时"
                    f"覆盖全池（missing_policy={rule['step3_gates']['herg'].get('missing_policy')}）。"
                    f"本层用 herg_gate == 1 过滤会把缺值行静默丢弃、等价于判不合格；"
                    f"请先修门作用域/预测表覆盖，或显式定义 not_evaluated 在选名单层的口径后重跑。")
            base = base[base.herg_gate == 1]
        log.info(f"候选基座（类药性 ∧ 姿态 pass ∧ maxTc<0.9 ∧ hERG[{herg_status}]）"
                 f"= {len(base)}")
        # ---- novelty 门（规则集 step3_gates.known_compound_novelty.hard_filter）
        kl = ROOT / "03_gates/tables/known_ligand_similarity.csv"
        if not kl.is_file():
            raise RuntimeError("缺少 known_ligand_similarity.csv（先运行 09_gate_stack）")
        kn = pd.read_csv(kl)
        base = base.merge(kn[["cid", "known_exact_hit", "max_tc_known"]], on="cid", how="left")
        n_known_in_base = int(base.known_exact_hit.fillna(0).sum())
        base_gated = base[base.known_exact_hit.fillna(0) == 0].copy()
        base_ungated = base.copy()          # 未施加 novelty 门（用于逐臂计数剔除量）
        log.info(f"novelty 门：基座内与 ChEMBL 已知 A2A 配体精确重合 {n_known_in_base} 个 "
                 f"→ 剔除后 {len(base_gated)}；"
                 f"{sorted(base.loc[base.known_exact_hit == 1, 'cid'])}")

        summaries, trig, arm_cids = {}, [], {}
        for name, q, a in ARMS:
            pool = base_gated[(base_gated.qed >= q) & (base_gated.anchor_hits_4A >= a)].copy()
            n_known_excl_arm = int(base_ungated[(base_ungated.qed >= q)
                                                & (base_ungated.anchor_hits_4A >= a)]
                                   .known_exact_hit.fillna(0).sum())
            res = pk.greedy_select(pool, fps, n_total - n_ctrl, score_col="CNN_VS",
                                   tc_max=sel["pairwise_tc_max"],
                                   scaffold_cap=sel["scaffold_cap"]["primary"])
            cards = card_rows(m, res["picks_df"], neg, herg_status)
            write_csv(ROOT / f"04_selection/tables/list_{name}.csv", cards)
            write_csv(ROOT / f"04_selection/tables/decision_cards_{name}.csv", cards)
            arm_cids[name] = sorted(res["picks"])

            s = cards[cards.candidate_type == "candidate"]
            tc_max = pk.max_pairwise_tc(s, fps)
            summaries[name] = {
                "arm": name, "qed_min": q, "anchor_min": a, "pool_n": int(len(pool)),
                "n_known_excluded_in_arm": n_known_excl_arm,
                "n_selected_candidates": int(len(s)), "n_controls": int(len(cards) - len(s)),
                "list_size": int(len(cards)),
                "anchor_all_4": bool((s.anchor_hits_4A == 4).all()),
                "anchor_hits_min": int(s.anchor_hits_4A.min()),
                "max_pairwise_tc": tc_max, "pairwise_ok": bool(tc_max < sel["pairwise_tc_max"]),
                "n_distinct_scaffold": int(s.scaffold.nunique()),
                "n_distinct_scaffold_generic": int(s.scaffold_generic.nunique()),
                "chemotype_max": int(s.chemotype.value_counts().max()) if len(s) else 0,
                "chem_cap": res["chem_cap"], "risk_n": int(s.risk_am1bcc.sum()),
                "risk_cap": res["risk_cap"], "scaffold_cap": res["scaffold_cap"],
                "CNN_VS_range": [r4(s.CNN_VS.min()), r4(s.CNN_VS.max())],
                "QED_range": [r4(s.qed.min()), r4(s.qed.max())],
                "maxTc_range": [r4(s.maxTc_actives.min()), r4(s.maxTc_actives.max())],
                "herg_range": ([r4(s.herg_prob.min()), r4(s.herg_prob.max())]
                               if herg_status == "applied" else "pending（未补算）"),
                "herg_status": herg_status,
                "known_exact_hits": int(s.cid.map(m.set_index("cid")["ik"]).isin(ki).sum()),
                "head20_cids": s.nsmallest(20, "pick_order").cid.tolist(),
                "head20_CNN_VS_min": r4(s.nsmallest(20, "pick_order").CNN_VS.min()),
                "head20_mean_qed": r4(float(s.nsmallest(20, "pick_order").qed.mean())),
            }
            trig.append({"arm": name, "pool_n": int(len(pool)), "selected": int(len(res["picks"])),
                         **{f"skip_{k}": int(v) for k, v in res["skips"].items()},
                         **{f"rate_{k}": r4(v / max(len(pool), 1))
                            for k, v in res["skips"].items()}})
            log.info(f"{name}：池 {len(pool)} → 选满 {len(res['picks'])}（Tc 跳过 "
                     f"{res['skips'].get('tc_pairwise', 0)}、骨架跳过 "
                     f"{res['skips'].get('scaffold_cap', 0)}、化学型跳过 "
                     f"{res['skips'].get('chem_cap', 0)}、风险档跳过 "
                     f"{res['skips'].get('risk_cap', 0)}）；成对 Tc 上限 {tc_max}；"
                     f"唯一骨架 {summaries[name]['n_distinct_scaffold']}")

        for cap in sel["scaffold_cap"]["sensitivity"]:
            pool = base_gated[(base_gated.qed >= 0.50) & (base_gated.anchor_hits_4A >= 4)].copy()
            res = pk.greedy_select(pool, fps, n_total - n_ctrl, score_col="CNN_VS",
                                   tc_max=sel["pairwise_tc_max"], scaffold_cap=cap)
            cards = card_rows(m, res["picks_df"], neg, herg_status)
            name = f"sens_scaffoldcap{cap}_QED0.50_anchor4"
            write_csv(ROOT / f"04_selection/tables/list_{name}.csv", cards)
            arm_cids[name] = sorted(res["picks"])
            s = cards[cards.candidate_type == "candidate"]
            summaries[name] = {
                "arm": name, "pool_n": int(len(pool)), "scaffold_cap": cap,
                "n_selected_candidates": int(len(s)),
                "max_pairwise_tc": pk.max_pairwise_tc(s, fps),
                "n_distinct_scaffold": int(s.scaffold.nunique()),
                "herg_status": herg_status,
                "overlap_with_main": r4(len(set(res["picks"])
                                            & set(arm_cids["main_QED0.50_anchor4"]))
                                        / max(len(res["picks"]), 1)),
                "note": "骨架封顶敏感性档（主口径池）"}
            log.info(f"{name}：选 {len(res['picks'])}，与主名单重叠 "
                     f"{summaries[name]['overlap_with_main']:.2%}")

        pool = base_gated[(base_gated.qed >= 0.50) & (base_gated.anchor_hits_4A >= 4)].copy()
        free = pool.sort_values(["CNN_VS", "cid"], ascending=[False, True]).head(n_total - n_ctrl)
        summaries["unconstrained_top48"] = {
            "pool_n": int(len(pool)), "chemotype_max": int(free.chemotype.value_counts().max()),
            "n_distinct_scaffold": int(free.scaffold.nunique()),
            "note": "无任何名单层约束时的 CNN_VS 前 48 名（量化约束咬合度）"}

        # ---- 未设 novelty 门时的替代名单（诊断：量化该门造成的变化）
        pool_ungated = base[(base.qed >= 0.50) & (base.anchor_hits_4A >= 4)]
        res_ung = pk.greedy_select(pool_ungated, fps, n_total - n_ctrl, score_col="CNN_VS",
                                   tc_max=sel["pairwise_tc_max"],
                                   scaffold_cap=sel["scaffold_cap"]["primary"])
        diff_cids = sorted(set(res_ung["picks"]) - set(arm_cids["main_QED0.50_anchor4"]))
        _rank = {c: i + 1 for i, c in enumerate(res_ung["picks"])}
        ungate = {"n_different_from_main": len(diff_cids), "cids_only_in_ungated": diff_cids,
                  "ungated_rank_of_known": {c: _rank[c] for c in diff_cids if c in _rank},
                  "n_known_in_ungated_list": int(
                      base.loc[base.cid.isin(res_ung["picks"]), "known_exact_hit"].sum()),
                  "known_molecule_in_ungated": bool(
                      base.loc[base.cid.isin(res_ung["picks"]), "known_exact_hit"].sum() > 0),
                  "note": "未设 novelty 门时的主口径名单（仅用于诊断，非交付名单）"}
        log.info(f"novelty 诊断：未设门时主名单有 {len(diff_cids)} 个分子不同"
                 f"（{diff_cids[:3]}）")

        write_csv(ROOT / "04_selection/tables/constraint_trigger_rates.csv", pd.DataFrame(trig))
        write_json(ROOT / "04_selection/reports/selection_summary.json",
                   {"herg_status": herg_status, "lists": summaries,
                    "novelty_gate": {"n_exact_known_in_base": n_known_in_base,
                                     "excluded_cids": sorted(base.loc[base.known_exact_hit == 1,
                                                                     "cid"].tolist()),
                                     "ungated_diagnostic": ungate},
                    "negative_controls": neg.to_dict("records"),
                    "chem_cap_rule": sel["chem_cap_rule"], "risk_cap_rule": sel["risk_cap_rule"],
                    "scaffold_cap_primary": sel["scaffold_cap"]["primary"]})

        labels = list(arm_cids)
        mat = [[pk.overlap_fraction(arm_cids[b], arm_cids[a]) for b in labels] for a in labels]
        plots.setup_cjk()
        plots.heatmap(mat, labels, labels, "各名单两两重叠率（候选部分）",
                      ROOT / "04_selection/figures/fig_list_overlap", cbar_label="重叠率")
        plots.grouped_barh(
            ["结构门（姿态 pass）", "锚点门 anchor=4", "活性相似度门 maxTc<0.9", "QED≥0.50", "hERG<0.5"],
            {"池层剔除数": [int((m.pose_gate != "pass").sum()),
                           int((m.anchor_hits_4A < 4).sum()),
                           int((m.maxTc_actives >= 0.9).sum()),
                           int((m.qed < 0.50).sum()),
                           (int((m.herg_gate == 0).sum()) if herg_status == "applied" else 0)]},
            "各门在池层的剔除规模", ROOT / "04_selection/figures/fig_gate_impact")
        log.info("名单矩阵、决策卡、约束触发率写出完成")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
