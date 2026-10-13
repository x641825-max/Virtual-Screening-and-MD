#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 2a：RDKit 全量重算描述符 + 类药性 + 注释列（PAINS 交叉验证）。

对去重后的池（≈4,990）逐分子重算：MW/cLogP/HBD/HBA/TPSA/RotB/环数/fsp3/形式电荷、
PAINS（FilterCatalog 重算）、Brenk（标注）、Ro5/Veber 违规数、QED、CNS MPO（6 分项）、
化学型（2 SMARTS 组合）、Bemis-Murcko 骨架、AM1-BCC 风险档、ECFP4 指纹。

结构规范化（frozen_rule.step2_druglikeness.standardization）：描述符与门所描述的物种 =
对接/模拟所处理的物种。成盐输入（盐 + 反离子）会使 MW/cLogP/TPSA/QED 与 Ro5/Veber/BRENK
判定被反离子污染（实测本池 61 个多片段分子，MW 虚高中位 +46 Da、最大 +192 Da）。
本步先取**最大片段**（保留原电荷态、不中和）再算描述符，并把反离子与输入口径列为独立列
（`smiles_input` / `mw_input` / `n_fragments` / `n_fragments_standardized` / `counter_ion_smiles`）。
单片段分子原样通过（描述符逐位不变）；规范化的影响面单列 `tables/standardization_impact.csv`。

交叉校验：与上游表 mw/logp/hbd/hba/tpsa/rotb 逐值比对（应在舍入精度内一致；
多片段分子因口径不同而预期不等，单列报告）。
输出：02_druglikeness/tables/descriptors.csv、fingerprints_ecfp4.npz（位矩阵）、
      tables/standardization_impact.csv、reports/descriptor_crosscheck.json、
      figures/fig_descriptor_hist.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import chem_core as cc
from vendor import plots
from vendor.cns_mpo import cns_mpo, cns_mpo_bounded, regression_check
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json

UPSTREAM = ROOT / "00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv"
CROSSCHECK_COLS = {"MW": "mw", "cLogP": "logp", "HBD": "hbd", "HBA": "hba",
                   "TPSA": "tpsa", "RotB": "rotb"}


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    up = pd.read_csv(UPSTREAM, low_memory=False)
    pool = pd.read_csv(ROOT / "01_pose_qc/tables/pose_pool.csv", low_memory=False)
    with StepLog(ROOT, "05_descriptors") as log:
        mpo_check = regression_check(verbose=True)
        log.info(f"CNS MPO 回归（原文 Table 3 三例）：最大偏差 {mpo_check['max_abs_diff']}")

        rows, fps, bad, std_impact = [], [], [], []
        for r in pool.itertuples():
            mol_raw = Chem.MolFromSmiles(r.smiles_upstream)
            if mol_raw is None:
                bad.append(r.cid)
                continue
            mol, n_frag_std, counter_ions = cc.standardize_main_fragment(mol_raw)
            d = {"cid": r.cid, "smiles": cc.canon_smiles(mol),
                 "smiles_input": r.smiles_upstream, "ik": r.ik_upstream,
                 "CNN_VS": float(r.CNN_VS), "n_fragments": int(r.n_fragments),
                 "n_fragments_standardized": int(n_frag_std),
                 "counter_ion_smiles": counter_ions,
                 "stereo_source": r.stereo_source, "qc_class": r.qc_class}
            if n_frag_std > 1:
                # 输入口径（含反离子）的描述符：只对多片段分子多算一次，用于影响面留痕
                d_raw = cc.druglikeness(mol_raw)
                d["mw_input"] = d_raw["MW"]
                d["qed_input"] = d_raw["qed"]
                d["pains_input"] = d_raw["pains"]
                d["brenk_input"] = d_raw["brenk"]
            d.update(cc.druglikeness(mol))
            d["smiles_canonical"] = cc.canon_smiles(mol)
            d["chemotype"] = cc.chemotype(mol)
            d["scaffold"] = cc.scaffold(mol)
            d["scaffold_generic"] = cc.scaffold_generic(cc.canon_smiles(mol))
            u = up.loc[up.cid == r.cid]
            pk = float(u["basic_pka"].iloc[0]) if len(u) and pd.notna(u["basic_pka"].iloc[0]) else None
            cd = float(u["clogd74"].iloc[0]) if len(u) and pd.notna(u["clogd74"].iloc[0]) else None
            st = u["pka_status"].iloc[0] if len(u) and "pka_status" in u else None
            cd_lo = (float(u["clogd74_lower_80"].iloc[0])
                     if len(u) and pd.notna(u["clogd74_lower_80"].iloc[0]) else None)
            cd_hi = (float(u["clogd74_upper_80"].iloc[0])
                     if len(u) and pd.notna(u["clogd74_upper_80"].iloc[0]) else None)
            d.update(cns_mpo_bounded(d["cLogP"], cd, cd_lo, cd_hi, d["MW"], d["TPSA"],
                                     d["HBD"], pk, st))
            d.update({"clogd74": cd, "basic_pka": pk, "pka_status": st,
                      "clogd74_lower_80": cd_lo, "clogd74_upper_80": cd_hi})
            d["cns_mpo_ge4"] = int(d["cns_mpo"] >= rule["step2_druglikeness"]["cns_mpo"]["threshold"])
            lo, hi = rule["step2_druglikeness"]["tpsa_band"]["range"]
            d["tpsa_in_band"] = int(lo <= d["TPSA"] <= hi)
            rows.append(d)
            fps.append(cc.fp_ecfp4(mol))
        des = pd.DataFrame(rows)
        if bad:
            log.info(f"[警告] SMILES 解析失败 {len(bad)} 条：{bad[:5]}")

        # ---- 结构规范化影响面（逐分子留痕：不过→过、新入池、主池候选）
        _multi = des[des.n_fragments_standardized > 1].copy()
        if len(_multi):
            _multi["mw_delta"] = (_multi.MW - _multi.mw_input).round(4)
            _multi["qed_delta"] = (_multi.qed - _multi.qed_input).round(6)
            _multi["brenk_delta_hit"] = _multi.brenk_input - _multi.brenk
            _multi["pains_delta_hit"] = _multi.pains_input - _multi.pains
            _q = rule["step2_druglikeness"]["qed"]["primary"]
            _multi["qed_ge_primary_after"] = (_multi.qed >= _q).astype(int)
            _multi["qed_ge_primary_before"] = (_multi.qed_input >= _q).astype(int)
            _impact = _multi[["cid", "smiles_input", "smiles", "n_fragments",
                              "n_fragments_standardized", "counter_ion_smiles",
                              "mw_input", "MW", "mw_delta", "qed_input", "qed", "qed_delta",
                              "pains_input", "pains", "brenk_input", "brenk",
                              "qed_ge_primary_before", "qed_ge_primary_after"]].copy()
            _impact = _impact.sort_values("mw_delta", ascending=False)
            write_csv(ROOT / "02_druglikeness/tables/standardization_impact.csv", _impact)
            log.info(f"结构规范化：多片段分子 {len(_impact)} 个（最大片段口径）→ "
                     f"tables/standardization_impact.csv；MW 虚高中位 "
                     f"{float((_multi.MW - _multi.mw_input).median()):+.1f} Da、最大 "
                     f"{float((_multi.MW - _multi.mw_input).max()):+.1f} Da；"
                     f"QED 最大偏移 {float((_multi.qed - _multi.qed_input).abs().max()):.3f}")
        else:
            log.info("结构规范化：池内无多片段分子（跳过影响面文件）")

        # ---- 阈值列（口径以 frozen_rule 为准）
        d2 = rule["step2_druglikeness"]
        des["pains_pass"] = (des.pains == 0).astype(int)
        des["ro5_pass"] = (des.ro5_viol == 0).astype(int)
        des["veber_pass"] = (des.veber_viol == 0).astype(int)
        des["druglike_pass"] = (des.pains_pass & des.ro5_pass & des.veber_pass).astype(int)
        for q in [d2["qed"]["primary"]] + list(d2["qed"]["sensitivity"]):
            des[f"qed_ge_{str(q).replace('.', '')}"] = (des.qed >= q).astype(int)
        write_csv(ROOT / "02_druglikeness/tables/descriptors.csv", des)

        arr = np.zeros((len(fps), 2048), dtype=np.uint8)
        for i, f in enumerate(fps):
            arr[i] = cc.fp_to_array(f)
        np.savez_compressed(ROOT / "02_druglikeness/tables/fingerprints_ecfp4.npz",
                            cid=des.cid.to_numpy(), bits=arr)

        # ---- 去盐后的身份冗余：上游去重键是**含盐** SMILES 的 InChIKey，故「盐 / 游离碱」
        # 同物种对会同时留在池内。本包**不重新去重**（会改变池规模口径，属另一项口径决定），
        # 但把冗余显式登记：选序的骨架 ≤1 与成对 Tc<0.65 约束使同一名单内不会同时出现两者
        # （实测 0 例），下游核对与后续去重决策都需要这个数。
        des["ik_main"] = [cc.inchi_key(Chem.MolFromSmiles(s)) for s in des.smiles]
        _vc = des.ik_main.value_counts()
        des["main_fragment_dup_group_n"] = des.ik_main.map(_vc).fillna(1).astype(int)
        _dup = des[des.main_fragment_dup_group_n > 1]
        if len(_dup):
            log.info(f"去盐后身份冗余：{_dup.ik_main.nunique()} 组 / {len(_dup)} 个分子"
                     f"（盐-游离碱同物种对；本包不重新去重，仅登记）")

        # ---- 交叉校验（与上游列）：单片段分子严格比对，多片段分子预期不等（口径不同）
        cc_rows = []
        merged = des.merge(up[["cid"] + list(CROSSCHECK_COLS.values())], on="cid", how="left")
        _single = merged[merged.n_fragments_standardized <= 1]
        _multi_m = merged[merged.n_fragments_standardized > 1]
        for ours, theirs in CROSSCHECK_COLS.items():
            diff = (_single[ours] - _single[theirs]).abs()
            cc_rows.append({"descriptor": ours, "upstream_col": theirs,
                            "max_abs_diff": r4(float(diff.max())),
                            "mean_abs_diff": r4(float(diff.mean())),
                            "n_exact": int((diff < 1e-6).sum()), "n": int(len(diff)),
                            "scope": "单片段分子（严格对比）"})
            if len(_multi_m):
                _dm = (_multi_m[ours] - _multi_m[theirs]).abs()
                cc_rows.append({"descriptor": ours, "upstream_col": theirs,
                                "max_abs_diff": r4(float(_dm.max())),
                                "mean_abs_diff": r4(float(_dm.mean())),
                                "n_exact": int((_dm < 1e-6).sum()), "n": int(len(_dm)),
                                "scope": "多片段分子（含反离子口径，预期不等）"})
            _m_txt = (f"；多片段（去盐口径，预期不等）max|Δ|={cc_rows[-2]['max_abs_diff']}"
                      if len(_multi_m) else "")
            log.info(f"交叉校验 {ours:6s} vs {theirs:5s}：单片段（严格）max|Δ|="
                     f"{cc_rows[-1]['max_abs_diff']}{_m_txt}")
        up_pains = int(up.loc[up.cid.isin(des.cid), "pains"].notna().sum())
        cross = {"descriptor_crosscheck": cc_rows, "mpo_regression_vs_paper": mpo_check,
                 "pains_recount_hits": int(des.pains.sum()),
                 "pains_upstream_nonnull": up_pains,
                 "pains_note": "上游 pains 列全为空值 → 本包以 RDKit FilterCatalog 重算为唯一口径",
                 "n_molecules": int(len(des)),
                 "main_fragment_redundancy": {
                     "n_groups": int(_dup.ik_main.nunique()) if len(_dup) else 0,
                     "n_molecules": int(len(_dup)),
                     "n_groups_involving_multi_fragment": int(
                         _dup[_dup.n_fragments_standardized > 1].ik_main.nunique())
                     if len(_dup) else 0,
                     "note": "去盐后同一物种在池内出现两个条目（盐 / 游离碱）。上游去重键是含盐 "
                             "SMILES 的 InChIKey，抓不到这类冗余。本包不重新去重（会改变池规模"
                             "口径，属独立的口径决定），只登记冗余；选序的骨架 ≤1 与成对 "
                             "Tc<0.65 约束已使其不会同时出现在同一名单内（实测 0 例）。",
                 },
                 "standardization": {
                     "method": "largest_fragment（保留原电荷态，不中和）",
                     "n_multi_fragment": int(len(_multi)),
                     "n_standardized_mol_changed": int(
                         (_multi.smiles != _multi.smiles_input).sum()) if len(_multi) else 0,
                     "mw_inflation_median": (round(float((_multi.MW - _multi.mw_input).median()), 4)
                                             if len(_multi) else None),
                     "mw_inflation_max": (round(float((_multi.MW - _multi.mw_input).max()), 4)
                                          if len(_multi) else None),
                     "qed_shift_max_abs": (round(float((_multi.qed - _multi.qed_input).abs().max()), 6)
                                           if len(_multi) else None),
                     "impact_table": "02_druglikeness/tables/standardization_impact.csv",
                     "note": "描述符与门所描述的物种 = 对接/模拟所处理的物种；含反离子的输入会把"
                             "MW/cLogP/TPSA/QED 与 Ro5/Veber/BRENK 判定污染（例如甲酸反离子会命中"
                             "BRENK aldehyde 而误报）",
                 },
                 "chemotype_counts": {str(k): int(v) for k, v in
                                      des.chemotype.value_counts().items()},
                 "scaffold_unique": int(des.scaffold.nunique()),
                 "brenk_hits": int(des.brenk.sum()),
                 "cns_mpo_ge4": int(des.cns_mpo_ge4.sum()),
                 "tpsa_in_band": int(des.tpsa_in_band.sum())}
        write_json(ROOT / "02_druglikeness/reports/descriptor_crosscheck.json", cross)

        plots.setup_cjk()
        plots.hist(des.qed, "QED 分布（去重后池）", ROOT / "02_druglikeness/figures/fig_qed_hist",
                   "QED", bins=40)
        plots.hist(des.cns_mpo, f"CNS MPO 分布（≥4 者 {cross['cns_mpo_ge4']} 个）",
                   ROOT / "02_druglikeness/figures/fig_cnsmpo_hist", "CNS MPO", bins=30,
                   color=plots.GREEN)
        log.info(f"描述符表 → 02_druglikeness/tables/descriptors.csv（{len(des)} × {des.shape[1]}）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
