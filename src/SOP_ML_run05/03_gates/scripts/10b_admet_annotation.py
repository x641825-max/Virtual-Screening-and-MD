#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3d（分析步 10b）：ADMET-AI 注释层 + 描述符独立复算交叉核验。

性质：**分析步** —— 只产出标注列与证据，不改门、不改池、不改名单（与 08b 同规格）。

输入（全部在包内；冻结表哈希登记于 frozen_rule.step3_gates.admet_layer）：
  - 03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv   f3 库外部冻结表（49,963 × 124）
  - 02_druglikeness/tables/descriptors.csv              本包 RDKit 重算口径
  - 03_gates/tables/candidate_universe.csv              池与门叠加结果（4,990）
  - 04_selection/tables/list_main_QED0.50_anchor4.csv   主名单（用于 MD 前置复核清单）

产出：
  - 03_gates/admet/tables/admet_annotation_pool.csv       池 4,990 注释全表（105 原始列 + 派生列）
  - 03_gates/admet/tables/admet_annotation_cards.csv      名单层精选注释（交付层追加列用）
  - 03_gates/admet/tables/admet_descriptor_crosscheck.csv 逐分子描述符对账
  - 03_gates/admet/tables/admet_md_precheck.csv           MD 前置复核清单（交付分子中带 ≥1 项警告者）
  - 03_gates/admet/reports/admet_annotation.json
  - 03_gates/admet/reports/admet_descriptor_crosscheck.json
  - 03_gates/admet/ADMET_LAYER.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.io_utils import (StepLog, counts_dict, r4, r6, read_json, sha256_file,
                             write_csv, write_json, write_text)

OUT = ROOT / "03_gates/admet"

# 冻结表的 105 个 ADMET 列 → 交付注释列的短名（其余列保留原名存于池全表）
RENAME = {
    "ADMET_hERG": ("admet_hERG", None),
    "ADMET_hERG_drugbank_approved_percentile": ("admet_herg_drugbank_pct", None),
    "ADMET_hydrogen_bond_donors": ("admet_HBD", None),
    "ADMET_hydrogen_bond_acceptors": ("admet_HBA", None),
    "ADMET_DILI": ("admet_DILI", None),
    "ADMET_AMES": ("admet_AMES", None),
    "ADMET_ClinTox": ("admet_ClinTox", None),
    "ADMET_Carcinogens_Lagunin": ("admet_Carcinogens", None),
    "ADMET_BBB_Martins": ("admet_BBB", None),
    "ADMET_Pgp_Broccatelli": ("admet_Pgp", None),
    "ADMET_Solubility_AqSolDB": ("admet_Solubility", None),
    "ADMET_Bioavailability_Ma": ("admet_Bioavailability", None),
    "ADMET_HIA_Hou": ("admet_HIA", None),
    "ADMET_Caco2_Wang": ("admet_Caco2", None),
    "ADMET_PPBR_AZ": ("admet_PPBR", None),
    "ADMET_Clearance_Hepatocyte_AZ": ("admet_Clearance_Hepatocyte", None),
    "ADMET_Half_Life_Obach": ("admet_Half_Life", None),
    "ADMET_VDss_Lombardo": ("admet_VDss", None),
    "ADMET_CYP2D6_Veith": ("admet_CYP2D6_Veith", None),
    "ADMET_CYP3A4_Veith": ("admet_CYP3A4_Veith", None),
    "ADMET_CYP2D6_Substrate_CarbonMangels": ("admet_CYP2D6_Substrate", None),
    "ADMET_CYP3A4_Substrate_CarbonMangels": ("admet_CYP3A4_Substrate", None),
    "ADMET_Lipinski": ("admet_Lipinski", None),
    "ADMET_QED": ("admet_QED", None),
    "ADMET_molecular_weight": ("admet_MW", None),
    "ADMET_logP": ("admet_logP", None),
    "ADMET_tpsa": ("admet_tpsa", None),
    "ADMET_PAINS_alert": ("admet_PAINS_alert", None),
    "ADMET_BRENK_alert": ("admet_BRENK_alert", None),
    "ADMET_NIH_alert": ("admet_NIH_alert", None),
}
DERIVED = ["admet_coverage", "admet_solubility_gate", "admet_flag_herg", "admet_flag_DILI",
           "admet_flag_AMES", "admet_flag_ClinTox", "admet_flag_Carc", "admet_struct_alert",
           "admet_liability_count", "admet_dev_ok", "admet_tier", "admet_herg_band",
           "admet_herg_margin", "admet_herg_margin_tight", "admet_md_precheck"]


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    spec = rule["step3_gates"]["admet_layer"]
    thr = rule["step3_gates"]["herg_threshold_analysis"]
    tol = float(spec["crosscheck"]["tolerance"])
    with StepLog(ROOT, "10b_admet_annotation") as log:
        # ---------- 输入与哈希
        raw_p = ROOT / spec["raw_table"]
        got = sha256_file(raw_p)
        if got != spec["raw_sha256"]:
            raise RuntimeError(f"ADMET 冻结表哈希不符：{got} != {spec['raw_sha256']}")
        adt = pd.read_csv(raw_p, low_memory=False)
        n_admet = sum(1 for c in adt.columns if c.startswith("ADMET_"))
        log.info(f"冻结表载入：{raw_p.name} {adt.shape[0]} 行 × {adt.shape[1]} 列"
                 f"（ADMET 列 {n_admet}）sha256={got[:16]}…")
        if adt.cid.duplicated().any():
            raise RuntimeError("冻结表 cid 不唯一")

        des = pd.read_csv(ROOT / "02_druglikeness/tables/descriptors.csv", low_memory=False)
        univ = pd.read_csv(ROOT / "03_gates/tables/candidate_universe.csv",
                           usecols=["cid", "CNN_VS", "qed", "herg_prob"], low_memory=False)
        pool = univ.merge(adt.drop(columns=[c for c in ("CNN_VS",) if c in adt.columns]),
                         on="cid", how="left")
        n_cov = int(pool.ADMET_hERG.notna().sum())
        log.info(f"池覆盖：{n_cov}/{len(pool)}；未覆盖 cid={sorted(pool.loc[pool.ADMET_hERG.isna(), 'cid'])}")
        if n_cov != int(spec["expected_pool_coverage"]):
            raise RuntimeError(f"池覆盖 {n_cov} != 规则期望 {spec['expected_pool_coverage']}")

        # ---------- 池内重标定切点（库内相对切点按本包池标定，不用大 run 的 54,932 口径）
        cp = spec["cutpoint_recalibration"]
        cut = {
            "herg": float(np.percentile(pool.ADMET_hERG.dropna(), cp["hERG_pool_percentile"])),
            "DILI": float(np.percentile(pool.ADMET_DILI.dropna(), cp["DILI_pool_percentile"])),
            "solubility": float(np.percentile(pool.ADMET_Solubility_AqSolDB.dropna(),
                                              cp["solubility_pool_percentile"])),
            "AMES": float(cp["absolute"]["AMES"]),
            "ClinTox": float(cp["absolute"]["ClinTox"]),
            "Carcinogens_Lagunin": float(cp["absolute"]["Carcinogens_Lagunin"]),
        }
        pool["admet_coverage"] = np.where(pool.ADMET_hERG.notna(), "covered", "not_covered")
        pool["admet_solubility_gate"] = (pool.ADMET_Solubility_AqSolDB >= cut["solubility"]).astype(int)
        pool["admet_flag_herg"] = (pool.ADMET_hERG >= cut["herg"]).astype(int)
        pool["admet_flag_DILI"] = (pool.ADMET_DILI >= cut["DILI"]).astype(int)
        pool["admet_flag_AMES"] = (pool.ADMET_AMES >= cut["AMES"]).astype(int)
        pool["admet_flag_ClinTox"] = (pool.ADMET_ClinTox >= cut["ClinTox"]).astype(int)
        pool["admet_flag_Carc"] = (pool.ADMET_Carcinogens_Lagunin >= cut["Carcinogens_Lagunin"]).astype(int)
        pool["admet_struct_alert"] = ((pool.ADMET_BRENK_alert >= 1) | (pool.ADMET_NIH_alert >= 1)
                                      | (pool.ADMET_PAINS_alert >= 1)).astype(int)
        pool["admet_liability_count"] = (pool[["admet_flag_herg", "admet_flag_DILI",
                                               "admet_flag_AMES", "admet_flag_ClinTox",
                                               "admet_flag_Carc"]].sum(axis=1)).astype(int)
        pool["admet_dev_ok"] = ((pool.ADMET_Lipinski >= 3) & (pool.admet_struct_alert == 0)
                                & (pool.admet_solubility_gate == 1)).astype(int)
        drop = ((pool.admet_liability_count >= 2) | (pool.ADMET_Lipinski <= 2)
                | ((pool.admet_struct_alert == 1) & (pool.admet_liability_count >= 1)))
        lead = ((pool.admet_liability_count == 0) & (pool.ADMET_hERG < cp["absolute"]["hERG_tier_pass"])
                & (pool.admet_dev_ok == 1))
        pool["admet_tier"] = np.where(drop, "建议淘汰", np.where(lead, "优先跟进", "需警惕"))
        soft, hard = float(thr["soft_line"]), float(thr["library_head_tradeoff_expected"])
        pool["admet_herg_band"] = np.where(pool.ADMET_hERG >= hard, "高",
                                           np.where(pool.ADMET_hERG >= soft, "中", "低"))
        pool.loc[pool.admet_coverage == "not_covered", "admet_herg_band"] = "not_covered"
        # 距门边距：把「一个概率值」变成可排序的分诊量。三档 band 在交付名单上恒为「低」
        # （软线 0.744 高于门 0.5），信息量为 0；边距则直接给出「离 0.5 有多远」，供 MD 前置
        # 复核按优先级排序（容差阈值取自 step3_gates.admet_layer.herg_margin_tight）。
        gate_t = float(rule["step3_gates"]["herg"]["threshold"])
        tight = float(cp.get("herg_margin_tight", 0.05))
        pool["admet_herg_margin"] = (gate_t - pool.ADMET_hERG).round(6)
        pool["admet_herg_margin_tight"] = np.where(
            pool.admet_coverage == "covered",
            (pool.admet_herg_margin < tight).astype(int), np.nan)
        pool["admet_md_precheck"] = ((pool.admet_liability_count >= 1) | (pool.admet_struct_alert == 1)
                                     | (pool.admet_dev_ok == 0)).astype(int)
        pool.loc[pool.admet_coverage == "not_covered", "admet_md_precheck"] = 0

        # ---------- 池注释全表（原始 105 列 + 派生列）
        keep_base = [c for c in adt.columns if c.startswith("ADMET_")]
        cols = ["cid", "CNN_VS", "rank"] + keep_base + DERIVED
        write_csv(OUT / "tables/admet_annotation_pool.csv", pool[cols])
        log.info(f"池注释全表：{len(pool)} 行 × {len(cols)} 列 → tables/admet_annotation_pool.csv")

        # ---------- 名单层精选注释（交付层追加列用；未覆盖者显式标 not_covered）
        cards = pd.DataFrame({"cid": pool.cid})
        for src, (dst, _) in RENAME.items():
            cards[dst] = pool[src]
        for c in DERIVED:
            cards[c] = pool[c]
        for c in cards.columns:
            if c != "admet_coverage":
                cards[c] = cards[c].where(cards.admet_coverage == "covered", "")
        write_csv(OUT / "tables/admet_annotation_cards.csv", cards)
        log.info(f"名单层精选注释：{len(cards)} 行 × {len(cards.columns)} 列 → "
                 f"tables/admet_annotation_cards.csv")

        # ---------- 描述符独立复算交叉核验（本包 RDKit 口径 vs ADMET-AI 内置实现）
        # 口径对齐：本包 05 步已按「最大片段」规范化（去盐）。ADMET-AI 侧的 MW/logP/TPSA
        # 等是按其**输入 SMILES**（含反离子）算的，故多片段分子上两口径**预期不等**，
        # 必须单列报告，不得混进严格容差统计（否则会把口径差异误读为回归）。
        x = des.merge(pool[["cid"] + list(RENAME)], on="cid", how="inner")
        _mf_col = ("n_fragments_standardized" if "n_fragments_standardized" in x.columns
                   else "n_fragments")
        x_single = x[x[_mf_col] <= 1]
        x_multi = x[x[_mf_col] > 1]
        pairs = [("MW", "ADMET_molecular_weight"), ("cLogP", "ADMET_logP"), ("TPSA", "ADMET_tpsa"),
                 ("HBD", "ADMET_hydrogen_bond_donors"), ("HBA", "ADMET_hydrogen_bond_acceptors"),
                 ("qed", "ADMET_QED")]
        xc = pd.DataFrame({"cid": x.cid})
        num_stats, num_stats_multi = {}, {}
        for ours, theirs in pairs:
            d = (x_single[ours].astype(float) - x_single[theirs].astype(float)).abs()
            num_stats[ours] = {"max_abs_diff": r6(float(d.max())) if len(d) else None,
                               "n_over_tol": int((d > tol).sum()), "n": int(len(d)),
                               "scope": "单片段分子（严格容差）"}
            if len(x_multi):
                dm = (x_multi[ours].astype(float) - x_multi[theirs].astype(float)).abs()
                num_stats_multi[ours] = {
                    "max_abs_diff": r6(float(dm.max())),
                    "n_over_tol": int((dm > tol).sum()), "n": int(len(dm)),
                    "scope": "多片段分子（本包已去盐，ADMET-AI 未去盐 → 预期不等）"}
            xc[f"d_{ours}"] = pd.concat(
                [d, (x_multi[ours].astype(float) - x_multi[theirs].astype(float)).abs()]
            ).reindex(x.index).round(9)
        xc["n_frag_std"] = x[_mf_col].to_numpy()
        rule_pairs = [("ro5", (x.ro5_viol == 0), (x.ADMET_Lipinski == 4)),
                      ("pains", (x.pains == 0), (x.ADMET_PAINS_alert == 0)),
                      ("brenk", (x.brenk == 1), (x.ADMET_BRENK_alert == 1))]
        rule_stats = {}
        for name, ours, theirs in rule_pairs:
            xc[f"agree_{name}"] = (ours == theirs).astype(int)
            rule_stats[name] = {
                "agreement_pct": r4(float((ours == theirs).mean()) * 100),
                "agreement_pct_single_fragment": (
                    r4(float((ours[x_single.index] == theirs[x_single.index]).mean()) * 100)
                    if len(x_single) else None),
                "pack_hits": int(ours.sum()), "admet_hits": int(theirs.sum()),
                "pack_only": int((ours & ~theirs).sum()), "admet_only": int((~ours & theirs).sum()),
                "n": int(len(x))}
        write_csv(OUT / "tables/admet_descriptor_crosscheck.csv", xc)
        xc_json = {
            "n_compared": int(len(x)),
            "tolerance": tol,
            "numeric": num_stats,
            "numeric_multi_fragment": num_stats_multi,
            "n_single_fragment": int(len(x_single)),
            "n_multi_fragment": int(len(x_multi)),
            "rule_agreements": rule_stats,
            "expected_differences": spec["crosscheck"]["expected_differences"],
            "all_numeric_within_tol": bool(all(v["n_over_tol"] == 0
                                              for v in num_stats.values())),
            "multi_fragment_note": "多片段分子的本包描述符按「最大片段」口径、ADMET-AI 按含反离子"
                                   "口径，逐列差值预期非零且方向已知（本包 MW 更小、logP/TPSA 相应更小）；"
                                   "该组不参与严格容差判定，其逐分子差值见 admet_descriptor_crosscheck.csv"
                                   "（列 n_frag_std）与 02_druglikeness/tables/standardization_impact.csv",
        }
        write_json(OUT / "reports/admet_descriptor_crosscheck.json", xc_json)
        log.info(f"描述符交叉核验：单片段 {len(x_single)} 个，数值列最大偏差 "
                 f"{max(v['max_abs_diff'] for v in num_stats.values())}（容差 {tol}）；"
                 f"多片段 {len(x_multi)} 个单列（去盐口径差异，预期不等）；"
                 f"规则一致率 Ro5 {rule_stats['ro5']['agreement_pct']}%、"
                 f"PAINS {rule_stats['pains']['agreement_pct']}%、"
                 f"BRENK {rule_stats['brenk']['agreement_pct']}%")

        # ---------- MD 前置复核清单（交付分子中带 ≥1 项警告）
        cards2 = cards.copy()
        cards2["CNN_VS"] = pool.CNN_VS.to_numpy()
        for c in DERIVED:
            cards2[c] = pd.to_numeric(cards2[c], errors="coerce")
        ml = pd.read_csv(ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv", low_memory=False)
        cand = ml.loc[ml.candidate_type == "candidate", "cid"].astype(str)
        pre = cards2[cards2.cid.isin(set(cand)) & (cards2.admet_md_precheck == 1)].copy()
        pre = pre.merge(des[["cid", "smiles"]], on="cid", how="left")
        pre_c = ["cid", "smiles", "CNN_VS", "admet_hERG", "admet_herg_margin",
                 "admet_herg_margin_tight", "admet_herg_drugbank_pct", "admet_herg_band",
                 "admet_DILI", "admet_AMES", "admet_ClinTox", "admet_Carcinogens",
                 "admet_Solubility", "admet_Lipinski", "admet_BRENK_alert", "admet_NIH_alert",
                 "admet_struct_alert", "admet_liability_count", "admet_tier"]
        pre = pre[pre_c].sort_values(["admet_liability_count", "CNN_VS"],
                                     ascending=[False, False]).reset_index(drop=True)
        write_csv(OUT / "tables/admet_md_precheck.csv", pre)
        log.info(f"MD 前置复核清单：{len(pre)}/{len(cand)} 个交付候选带 ≥1 项警告 → "
                 f"tables/admet_md_precheck.csv")

        # ---------- 门语义证据：0.5 这个数在分布上意味着什么
        # 与 08b（模型依赖度）互补：08b 回答「换模型会怎样」，本块回答「这个阈值本身处在
        # 什么位置」——库内/池内百分位、已上市药百分位映射、以及交付名单的贴线情况。
        _near = pool.ADMET_hERG.between(gate_t - 0.01, gate_t + 0.01)
        _dpc = pool.loc[_near, "ADMET_hERG_drugbank_approved_percentile"]
        _cand = pool[pool.cid.isin(set(cand))]
        _m = _cand.admet_herg_margin.dropna()
        sem = {
            "gate_threshold": gate_t,
            "library_pct_at_or_above_gate": r4(float((adt.ADMET_hERG >= gate_t).mean()) * 100),
            "pool_pct_at_or_above_gate": r4(float((pool.ADMET_hERG >= gate_t).mean()) * 100),
            "pool_empirical_percentile_of_gate": r4(float((pool.ADMET_hERG < gate_t).mean()) * 100),
            "approved_drug_percentile_at_gate_pool_near": (r4(float(_dpc.median()))
                                                           if len(_dpc) else None),
            "approved_drug_percentile_at_gate_library_near": r4(float(
                adt.loc[adt.ADMET_hERG.between(gate_t - 0.01, gate_t + 0.01),
                        "ADMET_hERG_drugbank_approved_percentile"].median())),
            "margin_tight_threshold": tight,
            "delivered_margin_min": (r6(float(_m.min())) if len(_m) else None),
            "delivered_margin_median": (r6(float(_m.median())) if len(_m) else None),
            "delivered_n_margin_lt_tight": int((_m < tight).sum()),
            "delivered_n_margin_lt_0.10": int((_m < 0.10).sum()),
            "delivered_n_in_half_band_0.40_0.50": int(
                _cand.ADMET_hERG.between(0.40, gate_t, inclusive="left").sum()),
            "note": "0.5 是未校准的模型分诊线：库内 66% 的分子在线上，并非化学上的安全界；"
                    "它大致对应已上市药在 ADMET-AI hERG 头下的中上段（≈ 62.8 百分位），"
                    "因此「过门」不得转述为「已判定安全」，只能转述为「ADMET-AI 口径下低风险标记」。",
        }
        log.info(f"门语义：库内 ≥{gate_t} = {sem['library_pct_at_or_above_gate']}%；"
                 f"已上市药百分位 ≈ {sem['approved_drug_percentile_at_gate_library_near']}；"
                 f"交付名单最小边距 {sem['delivered_margin_min']}（<{tight} 的 "
                 f"{sem['delivered_n_margin_lt_tight']} 个）")

        # ---------- 报告
        trig = {k: r4(float(pool[f"admet_flag_{k}"].mean()) * 100)
                for k in ("herg", "DILI", "AMES", "ClinTox", "Carc")}
        trig["struct_alert"] = r4(float(pool.admet_struct_alert.mean()) * 100)
        rep = {
            "role": spec["role"],
            "provenance": {"raw_table": spec["raw_table"], "raw_sha256": got,
                           "rows": int(adt.shape[0]), "columns": int(adt.shape[1]),
                           "admet_columns": int(n_admet), "producer": spec["producer"],
                           "regenerate_script": spec["regenerate_script"]},
            "coverage": {"pool_n": int(len(pool)), "covered": n_cov,
                         "not_covered": sorted(pool.loc[pool.admet_coverage == "not_covered", "cid"]),
                         "note": spec["coverage_note"]},
            "cutpoints_recalibrated_on_pool": {k: r6(v) for k, v in cut.items()},
            "cutpoints_reference_from_external_run": spec["reference_cutpoints_from_external_run"],
            "trigger_rates_pct": trig,
            "tier_counts": counts_dict(pool.admet_tier),
            "herg_band_counts_pool": counts_dict(pool.admet_herg_band),
            "herg_gate_semantics": sem,
            "liability_count_distribution": counts_dict(pool.admet_liability_count),
            "delivered_list": {"n_candidates": int(len(cand)),
                              "md_precheck_n": int(len(pre)),
                              "md_precheck_cids": pre.cid.astype(str).tolist(),
                              "tier_counts": counts_dict(
                                  pool.loc[pool.cid.isin(cand), "admet_tier"]),
                              "herg_band_counts": counts_dict(
                                  pool.loc[pool.cid.isin(cand), "admet_herg_band"]),
                              "herg_margin_min": sem["delivered_margin_min"],
                              "herg_margin_median": sem["delivered_margin_median"],
                              "n_margin_lt_tight": sem["delivered_n_margin_lt_tight"],
                              "n_margin_lt_0.10": sem["delivered_n_margin_lt_0.10"],
                              "margin_tight_cids": _cand.loc[_cand.admet_herg_margin < tight,
                                                             "cid"].astype(str).tolist(),
                              "herg_median": r6(float(pool.loc[pool.cid.isin(cand),
                                                             "ADMET_hERG"].median()))},
            "annotation_columns": spec["annotation_columns"],
            "non_interference": {
                "statement": "本步为分析步：不修改任何门（含 hERG<0.5）、候选池与名单成员；"
                             "注释列仅在交付层以追加列形式出现",
                "gate_threshold": rule["step3_gates"]["herg"]["threshold"],
                "list_size_main": rule["step4_selection"]["list_size"],
            },
        }
        write_json(OUT / "reports/admet_annotation.json", rep)

        md = [f"# ADMET-AI 注释层（步骤 3d · 分析步 10b）", "",
              f"- 冻结表：`{spec['raw_table']}`（{adt.shape[0]:,} 行 × {adt.shape[1]} 列，"
              f"{n_admet} 个 ADMET 列；sha256 `{got[:16]}…`）",
              f"- 生产者：{spec['producer']}",
              f"- 再生成：`{spec['regenerate_script']}`（手工执行、不入链）", "",
              "## 覆盖", "",
              f"- 本包池 4,990 分子全部有预测值；2 个阴性对照（run01 5,098 富集库来源）"
              f"标 `admet_coverage=not_covered`，不写空值。", "",
              "## 池内重标定切点（库内相对切点按本包池标定，避免照搬外部库口径）", "",
              "| 切点 | 本包池值 | 外部库原值 |", "|---|---|---|"]
        ref = spec["reference_cutpoints_from_external_run"]
        for k, zh in [("herg", "hERG 高危线"), ("DILI", "DILI 高危线"), ("solubility", "溶解度下限")]:
            md.append(f"| {zh} | {r6(cut[k])} | {ref[k]} |")
        for k in ("AMES", "ClinTox", "Carcinogens_Lagunin"):
            md.append(f"| {k} 绝对线 | {r6(cut[k])} | {ref[k]} |")
        md += ["", "## 池内触发率（%）", "",
               "| 端点 | 触发率 |", "|---|---|"]
        for k, v in trig.items():
            md.append(f"| {k} | {v} |")
        md += ["", "## 分诊层（池 4,990）", "",
               "| 层 | 分子数 |", "|---|---|"]
        for k, v in counts_dict(pool.admet_tier).items():
            md.append(f"| {k} | {v} |")
        md += ["", "## 描述符独立复算核验", "",
               f"- 数值列最大偏差："
               f"{max(v['max_abs_diff'] for v in num_stats.values())}（容差 {tol}）；"
               f"超容差项 {sum(v['n_over_tol'] for v in num_stats.values())} 个",
               f"- 判定一致率：Ro5 {rule_stats['ro5']['agreement_pct']}%、"
               f"PAINS {rule_stats['pains']['agreement_pct']}%、"
               f"BRENK {rule_stats['brenk']['agreement_pct']}%"
               f"（BRENK 差异为既有且已解释：{spec['crosscheck']['expected_differences']['BRENK']}）",
               "", "## MD 前置复核清单", "",
               f"- 交付 48 个候选中 **{len(pre)}** 个带 ≥1 项警告（详见 "
               f"`tables/admet_md_precheck.csv`）："
               + ("、".join(f"`{c}`" for c in pre.cid.astype(str).head(20))
                  + ("…" if len(pre) > 20 else "") if len(pre) else "无"),
               "", "## hERG 阈值门语义（0.5 在分布上意味着什么）", "",
               f"| 量 | 值 |", "|---|---|",
               f"| 库内 ≥ 门的分率 | {sem['library_pct_at_or_above_gate']}% |",
               f"| 本包池 ≥ 门的分率 | {sem['pool_pct_at_or_above_gate']}% |",
               f"| 门 ≈ 已上市药百分位 | {sem['approved_drug_percentile_at_gate_library_near']} |",
               f"| 交付名单最小距门边距 | {sem['delivered_margin_min']} |",
               f"| 边距 < {tight} 的候选数 | {sem['delivered_n_margin_lt_tight']} |",
               "",
               "> 三档 band（软线 0.744 / 硬线候选 0.80）在交付名单上恒为「低」；"
               "相比之下**距门边距**是可排序的分诊量，MD 前人工复核应按边距升序取用。",
               "",
               "> 本层非判据：不改门、不改池、不改名单；预测值非实验数据，"
               "采购/合成/立项前需体外确认。"]
        write_text(OUT / "ADMET_LAYER.md", "\n".join(md) + "\n")
        log.info(f"注释层完成：ADMET_LAYER.md + 2 个 JSON（非判据，门与名单未改动）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
