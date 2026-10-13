#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 5：稳健性 A–F 臂（G 臂确定性重跑见 12_determinism.py）。

A 打分扰动：CNN_VS + N(0, σ)，σ ∈ {0.05, 0.10, 0.15}，各 100 次重抽样选序 → 与主名单重叠
B 阈值矩阵：QED {0.40,0.50,0.60,0.70} × maxTc {0.70,0.80,0.90,无} 的 4×4 池与名单重叠
C 单轴敏感性：hERG / anchor / QED / maxTc 一次只动一轴
D 指纹替换：ECFP6、MACCS166 重跑「成对相似度约束」→ 只报重叠率
E 并列带宽：把第 48 名分数 −band 内的并列区随机置换，band ∈ {0.01, 0.05, 0.25}，各 100 次
F Butina 聚类：cutoff 0.35，簇代表 = 簇内 CNN_VS 最高；看主名单是否每簇 ≤1
输出：05_robustness/reports/robustness.json + tables/*.csv + figures/*.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import json
import numpy as np
import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import chem_core as cc
from vendor import picker as pk
from vendor import plots
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json
from vendor.master import build_master


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    rob = rule["step5_robustness"]
    seed = int(rule["seeds"]["PROJECT_SEED"])
    reps = int(rule["seeds"]["robustness_reps"])
    sel = rule["step4_selection"]
    n_total, n_ctrl = sel["list_size"], sel["n_negative_controls"]
    cap = sel["scaffold_cap"]["primary"]
    tcmax = sel["pairwise_tc_max"]

    with StepLog(ROOT, "11_robustness") as log:
        m = build_master(ROOT, log=log)
        herg_status = str(m.herg_status.iloc[0])
        base = m[(m.druglike_pass == 1) & (m.pose_gate == "pass") & (m.maxTc_actives < 0.9)
                 & (m.qed >= 0.50) & (m.anchor_hits_4A >= 4)].copy()
        if herg_status == "applied":
            base = base[base.herg_gate == 1]
        kn = pd.read_csv(ROOT / "03_gates/tables/known_ligand_similarity.csv")
        n_known = int(base.cid.map(kn.set_index("cid").known_exact_hit).fillna(0).sum())
        base = base[base.cid.map(kn.set_index("cid").known_exact_hit).fillna(0) == 0].copy()
        log.info(f"主口径池 = {len(base)}（novelty 门剔除已知物 {n_known} 个）")

        from rdkit import Chem, DataStructs
        fps4, fps6, fpm = {}, {}, {}
        for r in m.itertuples():      # 全池建指纹：C 轴（anchor≥3 / hERG off）会出现主池外的分子
            mol = Chem.MolFromSmiles(r.smiles)
            fps4[r.cid] = cc.fp_ecfp4(mol)
            fps6[r.cid] = cc.fp_ecfp6(mol)
            fpm[r.cid] = cc.fp_maccs(mol)

        main_res = pk.greedy_select(base, fps4, n_total - n_ctrl, tc_max=tcmax, scaffold_cap=cap)
        main_list = list(main_res["picks"])
        main_head20 = list(main_res["picks"][:20])
        log.info(f"主名单候选 {len(main_list)}（本轮独立重算）")
        out = {"pool_n": int(len(base)), "main_list_n": len(main_list),
               "herg_status": herg_status, "novelty_excluded": n_known}

        def run(sub, score: np.ndarray, fps: dict, order_hint=None, tc=tcmax, capv=cap):
            """在给定分数向量下重跑选序（分数与候选表对齐）。"""
            s = sub.copy()
            s["pert_score"] = score
            res = pk.greedy_select(s, fps, n_total - n_ctrl, score_col="pert_score",
                                   tc_max=tc, scaffold_cap=capv)
            return res["picks"]

        order = base.sort_values(["CNN_VS", "cid"], ascending=[False, True]).reset_index(drop=True)
        scores0 = order.CNN_VS.to_numpy()

        # ---------------- A 打分扰动
        rowsA, curve = [], {}
        for sigma in rob["A_score_sigma"]["sigmas"]:
            rng = np.random.default_rng(seed)
            ov, ov20 = [], []
            for _ in range(reps):
                s = pk.perturb_gaussian(scores0, sigma, rng)
                got = run(order, s, fps4)
                ov.append(len(set(got) & set(main_list)) / len(main_list))
                ov20.append(len(set(got[:20]) & set(main_head20)) / 20)
            rowsA.append({"sigma": sigma, "reps": reps,
                          "mean_overlap": r4(float(np.mean(ov))),
                          "min_overlap": r4(float(np.min(ov))),
                          "p10_overlap": r4(float(np.percentile(ov, 10))),
                          "mean_head20_overlap": r4(float(np.mean(ov20))),
                          "min_head20_overlap": r4(float(np.min(ov20)))})
            curve[str(sigma)] = r4(float(np.mean(ov)))
            log.info(f"A σ={sigma}：平均重叠 {rowsA[-1]['mean_overlap']}，"
                     f"最小 {rowsA[-1]['min_overlap']}，头部20 {rowsA[-1]['mean_head20_overlap']}")
        write_csv(ROOT / "05_robustness/tables/robA_sigma_sweep.csv", pd.DataFrame(rowsA))
        out["armA_score_perturbation"] = {
            "rows": rowsA, "primary_sigma": 0.10,
            "judge_mean_overlap_min": rob["A_score_sigma"]["judge_mean_overlap_min"],
            "pass_primary": bool([r for r in rowsA if r["sigma"] == 0.10][0]["mean_overlap"]
                                 >= rob["A_score_sigma"]["judge_mean_overlap_min"])}

        # ---------------- B 阈值矩阵
        rowsB = []
        for q in rob["B_threshold_matrix"]["qed"]:
            for mt in rob["B_threshold_matrix"]["maxTc_actives"]:
                sub = base[(base.qed >= q)]
                if mt is not None:
                    sub = sub[sub.maxTc_actives < mt]
                sub = sub.sort_values(["CNN_VS", "cid"], ascending=[False, True]).reset_index(drop=True)
                got = run(sub, sub.CNN_VS.to_numpy(), fps4)
                rowsB.append({"qed_min": q, "maxTc_upper": mt, "pool_n": int(len(sub)),
                              "list_n": len(got),
                              "overlap_vs_main": pk.overlap_fraction(got, main_list)})
                log.info(f"B QED≥{q} / maxTc<{mt}：池 {len(sub)}，选 {len(got)}，"
                         f"重叠 {rowsB[-1]['overlap_vs_main']}")
        dfB = pd.DataFrame(rowsB)
        write_csv(ROOT / "05_robustness/tables/robB_matrix.csv", dfB)
        out["armB_threshold_matrix"] = {"rows": rowsB,
                                        "overlap_min": r4(float(dfB.overlap_vs_main.min())),
                                        "overlap_max": r4(float(dfB.overlap_vs_main.max()))}

        # ---------------- C 单轴
        rowsC = []
        axes = [("hERG 门", "off→on<0.5", None, None, None),
                ("anchor", "4→≥3", None, 3, None),
                ("QED", "0.50→0.60", 0.60, None, None),
                ("maxTc", "0.90→0.70", None, None, 0.70)]
        for name, change, q, a, mt in axes:
            sub = m[(m.druglike_pass == 1) & (m.pose_gate == "pass")].copy()
            sub = sub[sub.cid.map(kn.set_index("cid").known_exact_hit).fillna(0) == 0]
            sub = sub[sub.maxTc_actives < (0.90 if mt is None else mt)]
            sub = sub[sub.qed >= (0.50 if q is None else q)]
            sub = sub[sub.anchor_hits_4A >= (4 if a is None else a)]
            if name == "hERG 门":
                if herg_status != "applied":
                    rowsC.append({"axis": name, "change": change, "status": "pending（无 hERG 表）",
                                  "pool_n": None, "list_n": None, "overlap_vs_main": None,
                                  "note": "外部表缺失 → 该轴无法实测，报告标注为未评估"})
                    log.info(f"C {name}（{change}）：pending，未评估")
                    continue
                # 关闭 hERG 门选一遍（其余口径不变），与开启后的主名单比较
                sub_off = sub.sort_values(["CNN_VS", "cid"],
                                          ascending=[False, True]).reset_index(drop=True)
                got_off = run(sub_off, sub_off.CNN_VS.to_numpy(), fps4)
                sub = sub[sub.herg_gate == 1]
                sub = sub.sort_values(["CNN_VS", "cid"], ascending=[False, True]).reset_index(drop=True)
                got = run(sub, sub.CNN_VS.to_numpy(), fps4)
                rowsC.append({"axis": name, "change": change, "status": "measured",
                              "pool_n": int(len(sub)), "pool_n_off": int(len(sub_off)),
                              "list_n": len(got),
                              "overlap_vs_main": pk.overlap_fraction(got_off, main_list),
                              "head20_overlap_vs_main": r4(len(set(got_off[:20]) & set(main_head20)) / 20),
                              "n_removed_by_gate": int(len(sub_off) - len(sub)),
                              "note": "关 hERG 门后重选（池 %d）与主名单（池 %d）比较"
                                      % (len(sub_off), len(sub))})
                log.info(f"C {name}（{change}）：关池 {len(sub_off)} → 开池 {len(sub)}，"
                         f"名单重叠 {rowsC[-1]['overlap_vs_main']}")
                continue
            sub = sub.sort_values(["CNN_VS", "cid"], ascending=[False, True]).reset_index(drop=True)
            got = run(sub, sub.CNN_VS.to_numpy(), fps4)
            rowsC.append({"axis": name, "change": change, "status": "measured",
                          "pool_n": int(len(sub)), "list_n": len(got),
                          "overlap_vs_main": pk.overlap_fraction(got, main_list),
                          "head20_overlap_vs_main": r4(len(set(got[:20]) & set(main_head20)) / 20),
                          "note": ""})
            log.info(f"C {name}（{change}）：池 {len(sub)}，重叠 {rowsC[-1]['overlap_vs_main']}")
        write_csv(ROOT / "05_robustness/tables/robC_axis.csv", pd.DataFrame(rowsC))
        out["armC_one_axis"] = rowsC

        # ---------------- D 指纹替换（只报重叠）
        rowsD = []
        for label, fdict in (("ECFP6", fps6), ("MACCS166", fpm)):
            got = run(order, scores0, fdict)
            rowsD.append({"fingerprint": label, "list_n": len(got),
                          "overlap_vs_main": pk.overlap_fraction(got, main_list),
                          "head20_overlap_vs_main": r4(len(set(got[:20]) & set(main_head20)) / 20),
                          "note": "仅替换成对相似度约束所用指纹，其余口径不变"})
            log.info(f"D {label}：重叠 {rowsD[-1]['overlap_vs_main']}")
        write_csv(ROOT / "05_robustness/tables/robD_fingerprint.csv", pd.DataFrame(rowsD))
        out["armD_fingerprint"] = rowsD

        # ---------------- E 并列带宽
        rowsE = []
        for band in rob["E_tie_band"]["bands"]:
            rng = np.random.default_rng(seed + 1)
            ov, ov20 = [], []
            for _ in range(reps):
                s = pk.perturb_tieband(scores0, band, len(main_list), rng)
                got = run(order, s, fps4)
                ov.append(len(set(got) & set(main_list)) / len(main_list))
                ov20.append(len(set(got[:20]) & set(main_head20)) / 20)
            rowsE.append({"band": band, "reps": reps, "mean_overlap": r4(float(np.mean(ov))),
                          "min_overlap": r4(float(np.min(ov))),
                          "mean_head20_overlap": r4(float(np.mean(ov20))),
                          "n_in_band": int((scores0 >= np.sort(scores0)[::-1][len(main_list) - 1]
                                            - band).sum())})
            log.info(f"E band={band}：平均重叠 {rowsE[-1]['mean_overlap']}，"
                     f"带宽内 {rowsE[-1]['n_in_band']} 个")
        write_csv(ROOT / "05_robustness/tables/robE_band.csv", pd.DataFrame(rowsE))
        out["armE_tie_band"] = {"rows": rowsE, "primary": rob["E_tie_band"]["primary"]}

        # ---------------- F Butina 聚类
        pool_sorted = order
        fps_pool = [fps4[c] for c in pool_sorted.cid]
        cutoff = rob["F_butina"]["cutoff"]
        clusters = pk.cluster_butina_fast(fps_pool, cutoff)   # 等价于官方 Butina，但避免 5k 级 OOM
        cs = pk.cluster_summary(clusters, pool_sorted.cid.tolist(),
                                pool_sorted.CNN_VS.to_numpy())
        member2cl = {cid: i for i, cl in enumerate(clusters) for cid in
                     [pool_sorted.cid[j] for j in cl]}
        main_cl = [member2cl[c] for c in main_list if c in member2cl]
        rep_table = pd.DataFrame(cs["clusters"])
        rep_table.to_csv(ROOT / "05_robustness/tables/robF_clusters.csv", index=False)
        top50_rep = rep_table.merge(pool_sorted[["cid", "CNN_VS"]], left_on="representative",
                                    right_on="cid", how="left").nlargest(50, "CNN_VS")
        top50_rep.to_csv(ROOT / "05_robustness/tables/robF_rep_top50.csv", index=False)
        out["armF_butina"] = {
            "cutoff": cutoff, "n_clusters": cs["n_clusters"],
            "largest_cluster": cs["largest_cluster"], "n_singletons": cs["n_singletons"],
            "n_members": cs["n_members"],
            "expected": {"n_clusters": rob["F_butina"]["expected_clusters"],
                          "largest": rob["F_butina"]["expected_largest"],
                          "singletons": rob["F_butina"]["expected_singletons"]},
            "main_list_cluster_multiplicity_max": int(pd.Series(main_cl).value_counts().max()),
            "overlap_rep_top50_vs_main": pk.overlap_fraction(
                top50_rep.representative.tolist(), main_list),
            "note": "代表 = 簇内 CNN_VS 最高的成员"}
        log.info(f"F Butina {cutoff}：{cs['n_clusters']} 簇（最大 {cs['largest_cluster']}、"
                 f"单例 {cs['n_singletons']}），主名单簇重复度上限 "
                 f"{out['armF_butina']['main_list_cluster_multiplicity_max']}")

        # 全池对照（去重后池整体聚类，用于与既有 642 簇口径对账）
        des = pd.read_csv(ROOT / "02_druglikeness/tables/descriptors.csv", low_memory=False)
        fps_all, cids_all = [], []
        for r in des.sort_values("cid").itertuples():
            mol = Chem.MolFromSmiles(r.smiles)
            fps_all.append(cc.fp_ecfp4(mol))
            cids_all.append(r.cid)
        cl_all = pk.cluster_butina_fast(fps_all, cutoff)
        cs_all = pk.cluster_summary(cl_all, cids_all, des.sort_values("cid").CNN_VS.to_numpy())
        _k = min(600, len(fps_all))
        same = ({frozenset(c) for c in pk.cluster_butina(fps_all[:_k], cutoff)}
                == {frozenset(c) for c in pk.cluster_butina_fast(fps_all[:_k], cutoff)})
        out["armF_butina_full_pool"] = {
            "impl_equivalence_check": {"n": _k, "identical_partition_vs_rdkit_official": bool(same)},
            "pool_n": len(cids_all), "n_clusters": cs_all["n_clusters"],
            "largest_cluster": cs_all["largest_cluster"],
            "n_singletons": cs_all["n_singletons"],
            "reference_expected": {"n_clusters": 642, "largest": 20, "singletons": 446},
            "note": "对去重后全池（≈4,990）聚类，与既有 642 簇口径可比"}
        log.info(f"F（全池 {len(cids_all)}）：{cs_all['n_clusters']} 簇，最大 "
                 f"{cs_all['largest_cluster']}，单例 {cs_all['n_singletons']}"
                 f"（既有口径 642 / 20 / 446，本包未能复现，见 scan）")
        # 口径扫描：把「距离阈值」与「Tanimoto 阈值」两种常见写法都给出，便于与外部数字对账
        scan = []
        for cut in (0.20, 0.30, 0.50, 0.65):
            c = pk.cluster_butina_fast(fps_all, cut)
            cs_s = pk.cluster_summary(c, cids_all, des.sort_values("cid").CNN_VS.to_numpy())
            scan.append({"dist_cutoff": cut, "tc_neighbor_min": round(1 - cut, 2),
                         "n_clusters": cs_s["n_clusters"],
                         "largest_cluster": cs_s["largest_cluster"],
                         "n_singletons": cs_s["n_singletons"]})
            log.info(f"F 扫描：距离阈值 {cut}（Tc≥{1-cut:.2f}）→ {cs_s['n_clusters']} 簇，"
                     f"最大 {cs_s['largest_cluster']}，单例 {cs_s['n_singletons']}")
        out["armF_butina_full_pool"]["cutoff_scan"] = scan
        (ROOT / "05_robustness/reports/robF_cutoff_scan.json").write_text(
            json.dumps(scan, ensure_ascii=False, indent=1), encoding="utf-8")
        pd.DataFrame(cs_all["clusters"]).to_csv(
            ROOT / "05_robustness/tables/robF_clusters_full_pool.csv", index=False)

        write_json(ROOT / "05_robustness/reports/robustness.json", out)

        plots.setup_cjk()
        plots.line_ci([r["sigma"] for r in rowsA],
                      {"平均重叠": [r["mean_overlap"] for r in rowsA],
                       "最小重叠": [r["min_overlap"] for r in rowsA],
                       "头部20 平均重叠": [r["mean_head20_overlap"] for r in rowsA]},
                      "A 打分扰动 σ 扫描（100 次重抽样）",
                      ROOT / "05_robustness/figures/fig_robA_sigma", "σ（CNN_VS 噪声）", "重叠率",
                      refs={"判据 0.8": 0.8})
        plots.heatmap(dfB.pivot(index="qed_min", columns="maxTc_upper",
                                values="overlap_vs_main").to_numpy(),
                      [f"QED≥{v}" for v in sorted(dfB.qed_min.unique())],
                      [f"maxTc<{v}" if pd.notna(v) else "maxTc 无上限"
                       for v in sorted(dfB.maxTc_upper.unique(), key=lambda x: (pd.isna(x), x))],
                      "B 阈值矩阵：与主名单重叠率",
                      ROOT / "05_robustness/figures/fig_robB_matrix", cbar_label="重叠率")
        plots.line_ci([r["band"] for r in rowsE],
                      {"平均重叠": [r["mean_overlap"] for r in rowsE],
                       "头部20 平均重叠": [r["mean_head20_overlap"] for r in rowsE]},
                      "E 并列带宽扫描（100 次重抽样）",
                      ROOT / "05_robustness/figures/fig_robE_band", "带宽（分数）", "重叠率",
                      refs={"判据 0.8": 0.8})
        plots.hist([len(c) for c in clusters], f"F Butina 0.35 簇大小分布（{len(clusters)} 簇）",
                   ROOT / "05_robustness/figures/fig_robF_clusters", "簇大小", bins=24,
                   color=plots.AMBER)
        log.info("稳健性 A–F 完成 → 05_robustness/reports/robustness.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
