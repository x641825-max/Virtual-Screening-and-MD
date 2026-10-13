#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3a：活性相似度门 maxTc（ECFP4 最大 Tanimoto vs 479 个已知活性）。

参考集 = run01 5,098 带标签富集库 activity_label==1（479 个）。
主门 maxTc < 0.90；敏感性档 maxTc < 0.70。
同时输出最相似活性（nearest active）供决策卡交叉核对。
输出：03_gates/tables/maxTc_actives.csv、reports/maxTc_summary.json
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import chem_core as cc
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json

LIB = ROOT / "00_inputs/extracted/run01/02_enrichment_5098/library/quick_5100_library.csv"


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    g = rule["step3_gates"]["maxTc_actives"]
    primary, sens = g["primary"], g["sensitivity"]
    des = pd.read_csv(ROOT / "02_druglikeness/tables/descriptors.csv", low_memory=False)
    lib = pd.read_csv(LIB)
    with StepLog(ROOT, "07_maxTc") as log:
        act = lib[lib.activity_label == 1]
        log.info(f"参考活性集：{len(act)} 个（库 {len(lib)} 条，活性占比 "
                 f"{len(act) / len(lib):.4f}）")
        act_fps, act_ids = [], []
        for cid, s in zip(act.mol_id, act.smiles):
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            act_fps.append(cc.fp_ecfp4(m))
            act_ids.append(str(cid))
        log.info(f"可用活性指纹：{len(act_fps)}")

        rows = []
        for r in des.itertuples():
            m = Chem.MolFromSmiles(r.smiles)
            fp = cc.fp_ecfp4(m)
            sims = DataStructs.BulkTanimotoSimilarity(fp, act_fps)
            j = int(np.argmax(sims))
            rows.append({"cid": r.cid, "maxTc_actives": round(float(sims[j]), 4),
                         "nearest_active_id": act_ids[j],
                         "n_actives_tc_ge_07": int(sum(1 for s in sims if s >= 0.7)),
                         "n_actives_tc_ge_09": int(sum(1 for s in sims if s >= 0.9))})
        mt = pd.DataFrame(rows)
        for thr in [primary] + list(sens):
            mt[f"pass_maxTc_lt_{str(thr).replace('.', '')}"] = (mt.maxTc_actives < thr).astype(int)
        write_csv(ROOT / "03_gates/tables/maxTc_actives.csv", mt)

        hist = Counter(np.digitize(mt.maxTc_actives, np.arange(0.1, 1.0, 0.1)).tolist())
        summary = {
            "reference_actives": int(len(act_fps)),
            "n_molecules": int(len(mt)),
            "maxTc_percentiles": {f"p{int(p)}": r4(float(np.percentile(mt.maxTc_actives, p)))
                                  for p in (0, 5, 25, 50, 75, 95, 100)},
            "pass_primary_lt_0.9": int(mt["pass_maxTc_lt_09"].sum()),
            "pass_sensitivity_lt_0.7": int(mt["pass_maxTc_lt_07"].sum()),
            "n_with_any_active_ge_0.9": int((mt.n_actives_tc_ge_09 > 0).sum()),
            "histogram_0.1_bins": {str(k): int(v) for k, v in sorted(hist.items())},
        }
        write_json(ROOT / "03_gates/reports/maxTc_summary.json", summary)
        log.info(f"maxTc（vs 479 活性）分位：{summary['maxTc_percentiles']}")
        log.info(f"maxTc<0.90 通过 {summary['pass_primary_lt_0.9']}；"
                 f"<0.70 通过 {summary['pass_sensitivity_lt_0.7']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
