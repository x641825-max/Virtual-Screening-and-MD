#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 1a：姿态属性一致性 + 结构身份分级 + 坐标往返偏差。

对 5,000 个最优姿态逐条：
  ① 属性一致性：|CNNscore × CNNaffinity − CNN_VS| ≤ 5e-4；SDF 属性 vs 上游排名表逐值比对；
  ② 结构身份：姿态（3D）vs 上游 SMILES → identical / stereo_only / protonation_diff / graph_mismatch；
  ③ 坐标往返：molblock → RDKit → 重写 molblock → 再读入，逐原子最大偏差。
输出：01_pose_qc/tables/pose_identity_qc.csv、reports/pose_identity_summary.json
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
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json

RDLogger.DisableLog("rdApp.*")
PROPS = ["CNN_VS", "CNNscore", "CNNaffinity", "CNNaffinity_variance",
         "minimizedAffinity", "rank_score", "n_fragments"]


def coord_roundtrip(block: str):
    """返回 (最大偏差 Å, 重读原子数, 原原子数)"""
    lines = block.split("\n")
    ci = next((i for i, l in enumerate(lines) if l.rstrip().endswith("V2000")), None)
    if ci is None:
        return None, 0, 0
    na = int(lines[ci][0:3])
    orig = np.array([[float(x) for x in l.split()[:3]] for l in lines[ci + 1: ci + 1 + na]])
    m = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
    if m is None:
        return None, 0, na
    rt = Chem.MolToMolBlock(m).split("\n")
    ci2 = next((i for i, l in enumerate(rt) if l.rstrip().endswith("V2000")), None)
    na2 = int(rt[ci2][0:3])
    new = np.array([[float(x) for x in l.split()[:3]] for l in rt[ci2 + 1: ci2 + 1 + na2]])
    if na != na2:
        return None, na2, na
    return float(np.abs(orig - new).max()), na2, na


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    tol = rule["step1_pose_qc"]["prop_consistency"]["tolerance"]
    sdf = ROOT / "00_inputs/raw/inputs/top5000_best_pose.sdf"
    up = pd.read_csv(ROOT / "00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv",
                     low_memory=False).set_index("cid")
    with StepLog(ROOT, "02_pose_identity_qc") as log:
        rows, prop_diff = [], {p: [] for p in PROPS}
        n_blocks = 0
        for block, props, title in cc.iter_sdf_blocks(sdf):
            n_blocks += 1
            cid = props.get("cid", title)
            ref = up.loc[cid] if cid in up.index else None
            for p in PROPS:
                if ref is not None and p in ref.index and p in props:
                    prop_diff[p].append(abs(float(props[p]) - float(ref[p])))
            mol = cc.mol_from_block(block)
            rec = {"cid": cid, "sdf_title": title, "n_heavy_pose": None,
                   "cnn_vs_x_affinity_delta": None, "coord_max_dev": None,
                   "n_atoms_before": None, "n_atoms_after": None}
            if mol is not None:
                rec["n_heavy_pose"] = int(mol.GetNumHeavyAtoms())
            if "CNNscore" in props and "CNNaffinity" in props and "CNN_VS" in props:
                rec["cnn_vs_x_affinity_delta"] = r4(
                    float(props["CNNscore"]) * float(props["CNNaffinity"]) - float(props["CNN_VS"]))
            dev, n_after, n_before = coord_roundtrip(block)
            rec["coord_max_dev"] = None if dev is None else round(dev, 6)
            rec["n_atoms_before"], rec["n_atoms_after"] = n_before, n_after
            rec.update(cc.classify_pose(mol, ref["smiles"] if ref is not None else ""))
            rec["smiles_upstream"] = ref["smiles"] if ref is not None else None
            rec["ik_upstream"] = cc.inchi_key(Chem.MolFromSmiles(ref["smiles"])) if ref is not None else None
            rec["ik_pose"] = cc.inchi_key(mol) if mol is not None else None
            rec["CNN_VS"] = float(props["CNN_VS"]) if "CNN_VS" in props else None
            rec["n_fragments_prop"] = props.get("n_fragments")
            rows.append(rec)
            if n_blocks % 1000 == 0:
                log.info(f"...已核验 {n_blocks} 条")
        df = pd.DataFrame(rows)
        write_csv(ROOT / "01_pose_qc/tables/pose_identity_qc.csv", df)

        ident_bad = int((df.cnn_vs_x_affinity_delta.abs() > tol).sum())
        summary = {
            "n_poses": int(len(df)),
            "identity_class_counts": {str(k): int(v) for k, v in
                                      df.identity_class.value_counts().sort_index().items()},
            "qc_class_counts_technical": {str(k): int(v) for k, v in
                                          df.qc_class.value_counts().sort_index().items()},
            "stereo_source_counts": {str(k): int(v) for k, v in
                                     df.stereo_source.value_counts().sort_index().items()},
            "n_multifragment_template": int((df.n_fragments_template > 1).sum()),
            "n_fragments_template_dist": {str(k): int(v) for k, v in
                                          df.n_fragments_template.value_counts().sort_index().items()},
            "cnn_vs_x_affinity_violations": ident_bad,
            "cnn_vs_x_affinity_tolerance": tol,
            "coord_max_dev_max": r4(float(pd.to_numeric(df.coord_max_dev, errors="coerce").max())),
            "prop_max_abs_diff_vs_upstream": {p: r4(max(v)) if v else None
                                              for p, v in prop_diff.items()},
            "heavy_atom_delta_counts": {str(k): int(v) for k, v in
                                        df.heavy_atom_delta.value_counts().sort_index().items()},
            "protonation_diff_cids": sorted(df.loc[df.qc_class == "protonation_diff", "cid"].tolist()),
            "graph_mismatch_cids": sorted(df.loc[df.qc_class == "graph_mismatch", "cid"].tolist()),
            "pass": bool(ident_bad == 0 and (pd.to_numeric(df.coord_max_dev, errors="coerce")
                                             .fillna(0).max() < 1e-6)),
        }
        write_json(ROOT / "01_pose_qc/reports/pose_identity_summary.json", summary)
        log.info(f"身份分级（简化）：{summary['identity_class_counts']}")
        log.info(f"身份分级（技术）：{summary['qc_class_counts_technical']}")
        log.info(f"CNNscore×CNNaffinity vs CNN_VS 越界数 = {ident_bad}（容差 {tol}）")
        log.info(f"坐标往返最大偏差 = {summary['coord_max_dev_max']} Å")
        log.info(f"多片段（输入 SMILES 含 '.'）分子数 = {summary['n_multifragment_template']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
