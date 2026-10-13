#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 1c：16 残基接触指纹 + 4 锚点命中 + 盒子几何 → 姿态闸门（全 5,000 姿态）。

口径（与 run01 `02_enrichment_5098/scripts/contact_fp.py` 一致，frozen_rule.receptor）：
  16 残基 = V84 L85 F168 E169 M174 M177 N181 W246 L249 H250 N253 H264 L267 M270 Y271 I274
  接触判据 = 配体重原子与残基重原子最近距离 < 4.0 Å
  锚点 = N253 / E169 / F168 / W246，命中数 anchor_hits_4A ∈ {0..4}
  闸门 = reject（出盒 或 锚点 ≤1）/ warn（锚点 = 2）/ pass（锚点 ≥3）
另做三件校验：
  ⓐ 与 run01 `contact_fingerprint_top200.csv` 逐距离对齐（期望最大偏差 0.00 Å）；
  ⓑ 5,098 带标签集的 IFP（每分子取 CNN_VS 最大姿态）→ 供锚点门证据链复核；
  ⓒ **共晶配体正对照**：以 4EIY 共晶配体 ZMA（ZM241385）在同一受体副本上重算
     16 残基接触指纹与锚点命中，给出参照读数（锚点距离、真实接触集、盒子余量）；
     并对每个上传姿态算「与共晶接触集的 Jaccard」（`contact_jaccard_cocrystal`）。
     没有正对照时，「anchor=4」只表明姿态接近 4 个残基，无法判断它离晶体真值的距离；
     该项把判据从「接近性」补上「与已知结合模式的相似性」参照。
  ⓓ **方向性与质量档（只报告、不落门）**：对两个极性锚点（N253/E169）算重原子近似的
     接触方向角 `angle_deg_*`（无氢几何下的角度代理），并按
     `anchor_hits × Jaccard` 给 `pose_quality_tier`（A/B/C/D，阈值见规则
     step1_pose_qc.pose_quality）。档位**不改变** `pose_gate`（pass/warn/reject）语义，
     以免破坏 anchor≥3 敏感性臂的口径；它给 MD 优先级与人工复核提供排序量。
输出：01_pose_qc/tables/pose_gate.csv、contact_fp_16res.csv、labeled5098_contact_fp.csv、
      cocrystal_reference_contacts.csv
      01_pose_qc/reports/pose_gate_summary.json、01_pose_qc/figures/fig_pose_gate.png
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

from vendor import chem_core as cc
from vendor import plots
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json

RAW_SDF = ROOT / "00_inputs/raw/inputs/top5000_best_pose.sdf"   # 闸门在全 5,000 上传姿态上判定
RECEPTOR = ROOT / "00_inputs/extracted/run01/01_redocking/receptor/4EIY_protein.pdb"
COCRYSTAL = ROOT / "00_inputs/extracted/run01/01_redocking/receptor/4EIY_ligand.pdb"
REF_TOP200 = ROOT / "00_inputs/extracted/run01/02_enrichment_5098/triage/contact_fingerprint_top200.csv"
LAB_POSES = ROOT / "00_inputs/extracted/run01/02_enrichment_5098/poses_raw/poses.sdf"
LAB_RANK = ROOT / "00_inputs/extracted/run01/02_enrichment_5098/scores/gnina_ranking.csv"


def cocrystal_reference(log, res) -> tuple[dict, set]:
    """4EIY 共晶配体（ZMA）正对照：锚点距离、真实接触集、盒子余量。

    这是姿态判据的**参照读数**（不参与落门）：它回答「anchor=4 离晶体真值有多远」——
    共晶配体本身在同一受体副本上应命中全部 4 个锚点，且给出 16 残基中真正接触的那部分。
    """
    ref: dict = {"ligand": "ZMA (ZM241385)", "source": str(COCRYSTAL.relative_to(ROOT))}
    if not COCRYSTAL.is_file():
        ref["status"] = "missing"
        log.info("[提示] 未找到 4EIY 共晶配体 PDB，正对照参照读数不可得")
        return ref, set()
    xyz = cc.read_pdb_heavy_coords(COCRYSTAL, resname="ZMA")
    _xyz_el = cc.read_pdb_heavy(COCRYSTAL, resname="ZMA")
    if xyz.size == 0:
        ref["status"] = "empty"
        log.info("[提示] 共晶配体 PDB 未解出 ZMA 重原子，正对照跳过")
        return ref, set()
    row = cc.ifp_row(xyz, res)
    contacts = [n for n in (row["contacts_4A"].split(";") if row["contacts_4A"] else [])]
    anchor_d = {a: row[f"d_{a}"] for a in cc.ANCHORS}
    box = cc.box_geometry(xyz)
    half = cc.BOX_EDGE / 2.0
    import numpy as _np
    d_inf = float(_np.abs(xyz - _np.asarray(cc.BOX_CENTER)).max())
    ref.update({
        "status": "ok",
        "n_heavy": int(xyz.shape[0]),
        "anchor_distances_A": anchor_d,
        "anchor_hits_4A": int(row["anchor_hits_4A"]),
        "contact_residues_4A": contacts,
        "n_contact_residues_4A": len(contacts),
        "pocket16_not_contacted": [n for n in cc.POCKET16.values() if n not in contacts],
        "pocket16_extra_vs_cocrystal_n": len([n for n in cc.POCKET16.values()
                                              if n not in contacts]),
        "box_edge_A": cc.BOX_EDGE,
        "box_half_edge_A": round(half, 3),
        "ligand_max_abs_dev_from_center_A": round(d_inf, 3),
        "box_margin_A": round(half - d_inf, 3),
        "outside_box": int(box["outside_box"]),
        "design_note": "盒子余量 = 半边长 − 配体坐标相对盒中心的最大绝对偏差；规则设计为「配体 bbox + 4 Å"
                       "」量级，故余量 ≈ 4 Å 属预期",
        "polar_angle_deg": cc.anchor_directionality(_xyz_el[0], _xyz_el[1], res),
    })
    log.info("共晶正对照（ZMA）：锚点距离 " + ", ".join(f"{k}={v} Å" for k, v in anchor_d.items())
             + f"；命中 {row['anchor_hits_4A']}/4；接触残基 {len(contacts)} 个；盒子余量 "
             f"{ref['box_margin_A']} Å")
    log.info("POCKET16 中在共晶里不接触的残基：" + ", ".join(ref["pocket16_not_contacted"]))
    return ref, set(contacts)


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    res = cc.load_receptor(RECEPTOR)
    with StepLog(ROOT, "04_pose_gate") as log:
        log.info("受体口袋残基原子数：" + ", ".join(f"{cc.POCKET16[k]}={len(v)}"
                                                    for k, v in sorted(res.items())))
        # ---- 共晶配体正对照（参照读数，不参与落门）
        cocrystal, ref_contacts = cocrystal_reference(log, res)
        if ref_contacts:
            write_csv(ROOT / "01_pose_qc/tables/cocrystal_reference_contacts.csv",
                      pd.DataFrame([{"residue": n,
                                     "distance_A": cocrystal["anchor_distances_A"].get(n)
                                     if n in cocrystal["anchor_distances_A"] else "",
                                     "is_anchor": int(n in cc.ANCHORS),
                                     "role": "共晶 ZMA 在 4.0 Å 内的接触残基"}
                                    for n in sorted(ref_contacts, key=lambda x: (
                                        list(cc.POCKET16.values()).index(x)))]))
        # ---- 全 5,000 上传姿态（闸门口径与规则集一致）
        _tier = rule["step1_pose_qc"].get("pose_quality", {})
        _tier_thr = float(_tier.get("jaccard_threshold", 0.60))
        rows = []
        for block, props, title in cc.iter_sdf_blocks(RAW_SDF):
            cid = props.get("cid", title)
            hc = cc.block_heavy(block)
            coords = hc[0] if hc else None
            elems = hc[1] if hc else []
            if coords is None:
                continue
            r = {"cid": cid, "CNN_VS": float(props["CNN_VS"]),
                 "n_heavy": int(coords.shape[0])}
            r.update(cc.ifp_row(coords, res))
            if ref_contacts:
                _c = set(r["contacts_4A"].split(";")) if r["contacts_4A"] else set()
                _u = _c | ref_contacts
                r["contact_jaccard_cocrystal"] = round(len(_c & ref_contacts) / max(len(_u), 1), 4)
                r["n_cocrystal_contacts_hit"] = len(_c & ref_contacts)
            r.update(cc.anchor_directionality(coords, elems, res))
            r.update(cc.box_geometry(coords))
            r["pose_gate"] = cc.pose_gate(r["anchor_hits_4A"], r["outside_box"])
            r["pose_quality_tier"] = cc.pose_quality_tier(
                r["anchor_hits_4A"], r.get("contact_jaccard_cocrystal"), _tier_thr)
            r["warn_or_reject_reason"] = (
                "" if r["pose_gate"] == "pass" else
                ("outside_box" if r["outside_box"] else f"anchor_hits={r['anchor_hits_4A']}"))
            rows.append(r)
        gate = pd.DataFrame(rows).sort_values("cid").reset_index(drop=True)
        write_csv(ROOT / "01_pose_qc/tables/pose_gate.csv", gate)
        write_csv(ROOT / "01_pose_qc/tables/contact_fp_16res.csv",
                  gate[["cid", "n_contacts_4A"] + [f"d_{n}" for n in cc.POCKET16.values()]])
        log.info(f"姿态闸门：{gate.pose_gate.value_counts().to_dict()}")
        log.info(f"锚点命中分布：{gate.anchor_hits_4A.value_counts().sort_index().to_dict()}")

        # ---- 5,098 带标签集 IFP（锚点门证据链）
        # ---- 5,098 带标签集 IFP（锚点门证据链）
        lab_stats = {}
        if LAB_POSES.is_file():
            best = {}
            for block, props, title in cc.iter_sdf_blocks(LAB_POSES):
                mid = props.get("cid") or props.get("mol_id") or title
                cv = float(props.get("CNN_VS", "nan"))
                if mid not in best or cv > best[mid][0]:
                    best[mid] = (cv, cc.block_coords(block)[0])
            rank = pd.read_csv(LAB_RANK)
            lab = pd.DataFrame([{"mol_id": m, "cnn_vs": best[m][0],
                                 "activity_label": int(rank.loc[rank.mol_id == m,
                                                                "activity_label"].iloc[0]),
                                 **cc.ifp_row(best[m][1], res)}
                                for m in best if (rank.mol_id == m).any()])
            write_csv(ROOT / "01_pose_qc/tables/labeled5098_contact_fp.csv", lab)
            a = lab[lab.activity_label == 1]
            d = lab[lab.activity_label == 0]
            lab_stats = {
                "n_labeled": int(len(lab)), "n_actives": int(len(a)), "n_decoys": int(len(d)),
                "anchor4_active_rate": r4(float((a.anchor_hits_4A == 4).mean())),
                "anchor4_decoy_rate": r4(float((d.anchor_hits_4A == 4).mean())),
                "anchor_ge3_active_rate": r4(float((a.anchor_hits_4A >= 3).mean())),
                "anchor_ge3_decoy_rate": r4(float((d.anchor_hits_4A >= 3).mean())),
                "anchor4_enrichment": r4(float((a.anchor_hits_4A == 4).mean()
                                               / max((d.anchor_hits_4A == 4).mean(), 1e-9))),
                "anchor4_actives_n": int((a.anchor_hits_4A == 4).sum()),
                "prevalence": r4(float(lab.activity_label.mean())),
                "p_active_given_anchor4": r4(float(
                    a[a.anchor_hits_4A == 4].shape[0] / max((lab.anchor_hits_4A == 4).sum(), 1))),
                "p_active_given_anchor3": r4(float(
                    a[a.anchor_hits_4A == 3].shape[0] / max((lab.anchor_hits_4A == 3).sum(), 1))),
                "anchor4_n_in_labeled": int((lab.anchor_hits_4A == 4).sum()),
                "anchor3_n_in_labeled": int((lab.anchor_hits_4A == 3).sum()),
            }
            log.info(f"5,098 带标签集（{lab_stats['n_labeled']}）：锚点 4/4 活性率 "
                     f"{lab_stats['anchor4_active_rate']:.4f} vs 诱饵 "
                     f"{lab_stats['anchor4_decoy_rate']:.4f}（富集 "
                     f"{lab_stats['anchor4_enrichment']:.2f}×）")
        # ---- 与 run01 top200 参考表逐距离对齐（同一 5,098 富集库，mol_id 对齐）
        align = {"n": 0, "max_abs_dist_diff": None}
        if REF_TOP200.is_file() and lab_stats:
            ref = pd.read_csv(REF_TOP200)
            labtab = pd.read_csv(ROOT / "01_pose_qc/tables/labeled5098_contact_fp.csv")
            dcols = [f"d_{n}" for n in cc.POCKET16.values()]
            m = ref.merge(labtab[["mol_id"] + dcols], on="mol_id", suffixes=("_ref", "_new"))
            if len(m):
                diff = np.abs(m[[c + "_ref" for c in dcols]].to_numpy()
                              - m[[c + "_new" for c in dcols]].to_numpy())
                pred = m[[f"d_{a}_new" for a in cc.ANCHORS]].lt(4.0).sum(axis=1).to_numpy()
                align = {"n": int(len(m)),
                         "max_abs_dist_diff": r4(float(np.nanmax(diff))),
                         "anchor_hits_match": bool((pred == m["anchor_hits_4A"].to_numpy()).all())}
                log.info(f"与 run01 top200 对齐：n={align['n']}，距离最大偏差 "
                         f"{align['max_abs_dist_diff']} Å，锚点命中全一致="
                         f"{align['anchor_hits_match']}")
        else:
            log.info("[提示] 参考表或带标签集 IFP 缺失，跳过对齐校验")

        if not lab_stats:
            log.info("[跳过] 未解包 5,098 姿态集（extract_large=false）或参考表缺失")

        summary = {
            "n_poses": int(len(gate)),
            "gate_counts": {str(k): int(v) for k, v in gate.pose_gate.value_counts().items()},
            "anchor_hits_counts": {str(k): int(v) for k, v in
                                   gate.anchor_hits_4A.value_counts().sort_index().items()},
            "outside_box_n": int(gate.outside_box.sum()),
            "warn_cids": sorted(gate.loc[gate.pose_gate == "warn", "cid"].tolist()),
            "reject_cids": sorted(gate.loc[gate.pose_gate == "reject", "cid"].tolist()),
            "n_contact_res_median": int(gate.n_contacts_4A.median()),
            "align_with_run01_top200": align,
            "labeled5098_anchor_evidence": lab_stats,
            "cocrystal_reference": cocrystal,
        }
        if ref_contacts and "contact_jaccard_cocrystal" in gate.columns:
            _pass = gate[gate.pose_gate == "pass"]
            _a4 = _pass[_pass.anchor_hits_4A == 4]
            _a3 = _pass[_pass.anchor_hits_4A == 3]
            summary["cocrystal_reference"]["jaccard_distribution"] = {
                "anchor4_median": float(_a4.contact_jaccard_cocrystal.median()),
                "anchor3_median": float(_a3.contact_jaccard_cocrystal.median()),
                "n_anchor4": int(len(_a4)), "n_anchor3": int(len(_a3)),
                "note": "与共晶接触集的 Jaccard（主口径内）；anchor=4 应显著高于 anchor=3",
            }
            log.info(f"与共晶接触集 Jaccard 中位：anchor=4 → "
                     f"{summary['cocrystal_reference']['jaccard_distribution']['anchor4_median']:.4f}，"
                     f"anchor=3 → "
                     f"{summary['cocrystal_reference']['jaccard_distribution']['anchor3_median']:.4f}")

        # ---- 方向性 / 质量档（报告口径；不参与落门）
        if "pose_quality_tier" in gate.columns:
            _tc = gate.pose_quality_tier.value_counts().to_dict()
            _dir = {
                "tier_counts_all5000": {k: int(_tc.get(k, 0)) for k in ("A", "B", "C", "D")},
                "jaccard_threshold": _tier_thr,
                "polar_anchor_contacts_n": int((gate.polar_contact_n > 0).sum()),
                "angle_deg_min_median": (None if gate.angle_deg_min.dropna().empty
                                         else float(gate.angle_deg_min.median())),
                "angle_deg_min_by_tier_median": {
                    k: (None if gate.loc[gate.pose_quality_tier == k, "angle_deg_min"].dropna().empty
                        else float(gate.loc[gate.pose_quality_tier == k, "angle_deg_min"].median()))
                    for k in ("A", "B", "C", "D")},
                "cocrystal_reference_angle_deg_min": (cocrystal.get("polar_angle_deg", {})
                                                      .get("angle_deg_min")),
                "n_below_cocrystal_minus_25": int(
                    (gate.angle_deg_min < float(cocrystal.get("polar_angle_deg", {})
                                                .get("angle_deg_min") or 0) - 25).sum()),
                "note": "极性锚点方向角为重原子近似（无氢几何）：θ=∠(X–D–A)。参照系必须用共晶配体"
                        "自身的读数（ZMA 实测 ≈126–134°），不能用「接近 180°」的理想氢键角——"
                        "重原子代理天然低于理想值。**实测该量不能区分 anchor=4 与 anchor=3**"
                        "（两档中位数接近），真正有区分度的是 contact_jaccard_cocrystal；"
                        "本项只作 MD 前人工复核的诊断列，不落门、不改变 pose_gate 语义",
            }
            summary["pose_directionality"] = _dir
            log.info(f"姿态质量档（全 5,000）：{_dir['tier_counts_all5000']}；"
                     f"极性锚点方向角中位 {_dir['angle_deg_min_median']}°"
                     f"（共晶参照 {_dir['cocrystal_reference_angle_deg_min']}°；"
                     f"按档中位 {_dir['angle_deg_min_by_tier_median']}——该量不区分 A/C 档）")
        write_json(ROOT / "01_pose_qc/reports/pose_gate_summary.json", summary)

        fig_steps = {"姿态输入": len(gate),
                     "pass": int((gate.pose_gate == "pass").sum()),
                     "warn（进留档）": int((gate.pose_gate == "warn").sum()),
                     "reject": int((gate.pose_gate == "reject").sum()),
                     "pass 且 anchor=4": int(((gate.pose_gate == "pass")
                                              & (gate.anchor_hits_4A == 4)).sum())}
        plots.setup_cjk()
        plots.funnel_barh(fig_steps, "姿态闸门与锚点（5,000 上传最优姿态）",
                          ROOT / "01_pose_qc/figures/fig_pose_gate")
        log.info("图 → 01_pose_qc/figures/fig_pose_gate.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
