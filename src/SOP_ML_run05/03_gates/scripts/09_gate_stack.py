#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3c：候选主表 + 门叠加统计 + 已知物核验。

产出：
  03_gates/tables/candidate_universe.csv  候选主表（描述符/姿态/hERG/maxTc/注释/上游 pKa-cLogD）
  03_gates/tables/gate_stack.csv          各臂门叠加计数
  03_gates/reports/gate_stack.json        漏斗 JSON
  03_gates/tables/known_ligand_similarity.csv  与 ChEMBL 已知 A2A 配体的最近邻（novelty 核验）
  03_gates/figures/fig_gate_stack.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import chem_core as cc
from vendor import plots
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json
from vendor.master import build_master

CHEMBL = ROOT / "00_inputs/extracted/screening_kit/07_independent_validation/independent_a2a_chembl.csv"


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    with StepLog(ROOT, "09_gate_stack") as log:
        m = build_master(ROOT, log=log)
        write_csv(ROOT / "03_gates/tables/candidate_universe.csv", m)   # novelty 列在下方核验后重写

        # ---- 门叠加（主口径链）
        chain = {
            "0_去重后池": int(len(m)),
            "1_类药性硬门（PAINS/Ro5/Veber）": int((m.druglike_pass == 1).sum()),
            "2_姿态闸门 pass": int(((m.druglike_pass == 1) & (m.pose_gate == "pass")).sum()),
            "3_anchor=4": int(((m.druglike_pass == 1) & (m.pose_gate == "pass")
                               & (m.anchor_hits_4A == 4)).sum()),
            "4_maxTc<0.9": int(((m.druglike_pass == 1) & (m.pose_gate == "pass")
                                & (m.anchor_hits_4A == 4) & (m.maxTc_actives < 0.9)).sum()),
            "5_QED>=0.50（主池）": int(((m.druglike_pass == 1) & (m.pose_gate == "pass")
                                       & (m.anchor_hits_4A == 4) & (m.maxTc_actives < 0.9)
                                       & (m.qed >= 0.50)).sum()),
        }
        if m.herg_status.iloc[0] == "applied":
            chain["6_hERG<0.5"] = int(((m.druglike_pass == 1) & (m.pose_gate == "pass")
                                       & (m.anchor_hits_4A == 4) & (m.maxTc_actives < 0.9)
                                       & (m.qed >= 0.50) & (m.herg_gate == 1)).sum())
        log.info(f"门叠加：{chain}")

        # ---- 各臂池（3 QED × 2 anchor）与四臂并集
        arm_rows = []
        pools = {}
        for q in (0.40, 0.50, 0.60, 0.70):
            for a in (4, 3):
                key = f"QED{q:.2f}_anchor{a}"
                sub = m[(m.druglike_pass == 1) & (m.pose_gate == "pass")
                        & (m.anchor_hits_4A >= a) & (m.maxTc_actives < 0.9) & (m.qed >= q)]
                pools[key] = set(sub.cid)
                arm_rows.append({"arm": key, "qed_min": q, "anchor_min": a,
                                 "pool_n": len(sub),
                                 "herg_status": str(m.herg_status.iloc[0]),
                                 "pool_n_herg_pass": (int((sub.herg_gate == 1).sum())
                                                      if m.herg_status.iloc[0] == "applied" else None)})
        four = ["QED0.60_anchor4", "QED0.50_anchor4", "QED0.40_anchor4", "QED0.60_anchor3"]
        union4 = set().union(*[pools[k] for k in four])
        write_csv(ROOT / "03_gates/tables/gate_stack.csv", pd.DataFrame(arm_rows))

        # ---- 已知物核验（novelty）
        kn = {"n_known": 0, "exact_inchikey_hits": 0, "n_tc_ge_0.9": 0, "max_tc": None,
              "top": []}
        if CHEMBL.is_file():
            kn_df = pd.read_csv(CHEMBL)
            known_smi = kn_df.smiles.dropna().unique().tolist()
            known_ik = set()
            kfps, kids = [], []
            for s in known_smi:
                mm = Chem.MolFromSmiles(s)
                if mm is None:
                    continue
                ik = cc.inchi_key(mm)
                if ik:
                    known_ik.add(ik)
                kfps.append(cc.fp_ecfp4(mm))
                kids.append(ik)
            rows = []
            for r in m.itertuples():
                mm = Chem.MolFromSmiles(r.smiles)
                fp = cc.fp_ecfp4(mm)
                sims = DataStructs.BulkTanimotoSimilarity(fp, kfps)
                j = int(np.argmax(sims))
                rows.append({"cid": r.cid, "ik": r.ik,
                             "known_exact_hit": int(r.ik in known_ik),
                             "max_tc_known": round(float(sims[j]), 4),
                             "nearest_known_ik": kids[j]})
            kk = pd.DataFrame(rows)
            write_csv(ROOT / "03_gates/tables/known_ligand_similarity.csv", kk)
            m = m.merge(kk[["cid", "known_exact_hit", "max_tc_known", "nearest_known_ik"]],
                        on="cid", how="left")
            write_csv(ROOT / "03_gates/tables/candidate_universe.csv", m)
            kn = {"n_known_smiles": len(kfps),
                  "exact_inchikey_hits": int(kk.known_exact_hit.sum()),
                  "n_tc_ge_0.9": int((kk.max_tc_known >= 0.9).sum()),
                  "max_tc": r4(float(kk.max_tc_known.max())),
                  "top": kk.nlargest(5, "max_tc_known")[["cid", "max_tc_known"]].to_dict("records")}
            log.info(f"已知物核验：精确重合 {kn['exact_inchikey_hits']}；"
                     f"Tc≥0.9 相似 {kn['n_tc_ge_0.9']}；最大 Tc {kn['max_tc']}")
        else:
            log.info("[跳过] 未找到 ChEMBL 独立已知配体集")

        out = {"gate_chain": chain, "arm_pools": arm_rows,
               "novelty_gate": {
                   "rule": "选序层剔除 InChIKey 与 ChEMBL 已知 A2A 配体精确重合者",
                   "n_pool_exact_known": int(m.known_exact_hit.sum()) if "known_exact_hit" in m else None,
                   "pool_exact_known_cids": (sorted(m.loc[m.known_exact_hit == 1, "cid"].tolist())
                                             if "known_exact_hit" in m else [])},
               "union_four_arms": len(union4),
               "hERG": {"status": str(m.herg_status.iloc[0]),
                        "n_pending": int((m.herg_status == "pending").sum())},
               "known_ligand_check": kn,
               "chemotype_pool_distribution": {str(k): int(v) for k, v in
                                               m.chemotype.value_counts().items()},
               "n_risk_am1bcc_pool": int(m.risk_am1bcc.sum()),
               "n_pose_warn_kept_out": int((m.pose_gate != "pass").sum())}
        write_json(ROOT / "03_gates/reports/gate_stack.json", out)

        plots.setup_cjk()
        plots.funnel_barh(chain, "门叠加漏斗（主口径链）",
                          ROOT / "03_gates/figures/fig_gate_stack")
        log.info(f"候选主表 → 03_gates/tables/candidate_universe.csv（{m.shape[0]} × {m.shape[1]}）")
        log.info(f"四臂并集 = {len(union4)}（hERG 补算范围）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
