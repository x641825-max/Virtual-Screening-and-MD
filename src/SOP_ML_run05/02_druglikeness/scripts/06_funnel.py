#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 2b：类药性漏斗（硬门 + 三档 QED + 标注列统计）。

口径：PAINS == 0 且 Ro5 违规 == 0 且 Veber 违规 == 0 为硬门；
QED 主口径 ≥0.50，敏感性 ≥0.60 / ≥0.40；CNS MPO ≥4 与 TPSA∈[90,140] 为标注列。
输出：02_druglikeness/reports/druglikeness_funnel.json、tables/druglikeness_columns.csv
      figures/fig_druglikeness_funnel.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import plots
from vendor.io_utils import StepLog, read_json, r4, write_csv, write_json


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    d2 = rule["step2_druglikeness"]
    des = pd.read_csv(ROOT / "02_druglikeness/tables/descriptors.csv", low_memory=False)
    with StepLog(ROOT, "06_funnel") as log:
        n0 = len(des)
        keep = des[des.druglike_pass == 1].copy()
        funnel = {"0_去重后池": n0,
                  "1_硬门剔除": n0 - len(keep),
                  "2_PAINS0_Ro5Veber零违规": len(keep)}
        attribution = {
            "pains_fail": int((des.pains_pass == 0).sum()),
            "ro5_fail": int((des.ro5_pass == 0).sum()),
            "veber_fail": int((des.veber_pass == 0).sum()),
            "ro5_fail_only": int(((des.ro5_pass == 0) & (des.veber_pass == 1)
                                  & (des.pains_pass == 1)).sum()),
            "veber_fail_only": int(((des.veber_pass == 0) & (des.ro5_pass == 1)
                                    & (des.pains_pass == 1)).sum()),
            "pains_fail_only": int(((des.pains_pass == 0) & (des.ro5_pass == 1)
                                    & (des.veber_pass == 1)).sum()),
        }
        qed_tiers = {}
        for q in [d2["qed"]["primary"]] + list(d2["qed"]["sensitivity"]):
            col = f"qed_ge_{str(q).replace('.', '')}"
            qed_tiers[f"QED>={q:.2f}"] = int(keep[col].sum())
            funnel[f"3_QED>={q:.2f}"] = int(keep[col].sum())
        # 违规项归因细节
        viol_detail = {}
        for col, name in (("ro5_viol", "Ro5"), ("veber_viol", "Veber")):
            vc = des.loc[des[col] > 0, col].value_counts().sort_index()
            viol_detail[name] = {str(k): int(v) for k, v in vc.items()}
        ann = {
            "CNS MPO>=4": int(keep.cns_mpo_ge4.sum()),
            "TPSA in [90,140]": int(keep.tpsa_in_band.sum()),
            "CNS MPO>=4 且 TPSA∈[90,140]": int((keep.cns_mpo_ge4 & keep.tpsa_in_band).sum()),
            "Brenk 命中": int(keep.brenk.sum()),
            # cLogD 80% 预测区间界定下的 CNS MPO（标注用）：悲观界 = cLogD 取区间上界
            "CNS MPO>=4（点估计）": int(keep.cns_mpo_ge4.sum()),
            "CNS MPO>=4（80% 区间悲观界）": int((keep.cns_mpo_lo80 >= 4).sum()),
            "CNS MPO>=4（80% 区间乐观界）": int((keep.cns_mpo_hi80 >= 4).sum()),
            "cLogD 区间宽度中位": r4(float((keep.clogd74_upper_80
                                            - keep.clogd74_lower_80).median())),
            "cLogD 区间使 CNS MPO 变动>=0.1 的分子": int(
                (keep.cns_mpo_pi80_span >= 0.1).sum()),
        }
        out = {"funnel": funnel, "attribution": attribution,
               "violation_count_distribution": viol_detail,
               "annotation_columns_on_druglike_pool": ann,
               "primary_qed": d2["qed"]["primary"],
               "coverage": {k: r4(v / n0) for k, v in funnel.items() if k != "0_去重后池"}}
        write_json(ROOT / "02_druglikeness/reports/druglikeness_funnel.json", out)
        write_csv(ROOT / "02_druglikeness/tables/druglikeness_columns.csv",
                  keep[["cid", "smiles", "CNN_VS", "MW", "cLogP", "HBD", "HBA", "TPSA",
                        "RotB", "qed", "pains", "ro5_viol", "veber_viol", "brenk",
                        "brenk_names", "cns_mpo", "cns_mpo_ge4", "cns_mpo_lo80",
                        "cns_mpo_hi80", "cns_mpo_pi80_span", "basic_pka", "clogd74",
                        "clogd74_lower_80", "clogd74_upper_80", "pka_status", "tpsa_in_band",
                        "chemotype", "scaffold", "risk_am1bcc", "n_aromN", "n_N"]])

        plots.setup_cjk()
        plots.grouped_barh(
            ["PAINS 命中", "Ro5 违规", "Veber 违规", "Ro5 或 Veber", "任一违规"],
            {"分子数": [attribution["pains_fail"], attribution["ro5_fail"],
                       attribution["veber_fail"],
                       int(((des.ro5_pass == 0) | (des.veber_pass == 0)).sum()),
                       n0 - len(keep)]},
            "类药性硬门剔除归因（去重后池）",
            ROOT / "02_druglikeness/figures/fig_druglikeness_funnel")
        log.info(f"漏斗：{funnel}")
        log.info(f"归因：{attribution}")
        log.info(f"标注列（药物样池）：{ann}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
