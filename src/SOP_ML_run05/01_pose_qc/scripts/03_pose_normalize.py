#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 1b：InChIKey 去重 + 姿态规范化 + 标注列。

规则（frozen_rule.step1_pose_qc.dedup）：同 InChIKey(上游 SMILES) 保留 CNN_VS 最大者，
并列取 cid 字典序最小；被剔除者逐条登记。
“规范化” 指：统一属性块（补 cid/smiles/ik/stereo_source/n_fragments/protonation_note/
qc_class）、按 cid 字典序重排、几何坐标保持原样（不重新嵌入、不做能量最小化）。
输出：01_pose_qc/outputs/SOP-05_poses_normalized.sdf
      01_pose_qc/tables/pose_pool.csv、pose_dedup_map.csv
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

from vendor import chem_core as cc
from vendor.io_utils import StepLog, read_json, write_csv, write_json

ANNOT_COLS = ["cid", "smiles_upstream", "ik_upstream", "stereo_source", "n_fragments_template",
              "protonation_note", "qc_class", "identity_class", "CNN_VS",
              "n_fragments_prop", "n_heavy_pose"]


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    exp_kept = rule["step1_pose_qc"]["dedup"]["expected_kept"]
    qc = pd.read_csv(ROOT / "01_pose_qc/tables/pose_identity_qc.csv", low_memory=False)
    with StepLog(ROOT, "03_pose_normalize") as log:
        # ---- 去重（先按 CNN_VS 降序、cid 升序定序，取每个 InChIKey 的首条）
        qc = qc.sort_values(["CNN_VS", "cid"], ascending=[False, True]).reset_index(drop=True)
        qc["dup_rank"] = qc.groupby("ik_upstream").cumcount()
        keep = qc[qc.dup_rank == 0].copy()
        dropped = qc[qc.dup_rank > 0].copy()
        keep["dedup_group_size"] = keep.ik_upstream.map(qc.ik_upstream.value_counts())
        log.info(f"姿态 {len(qc)} → 去重后 {len(keep)}（剔除 {len(dropped)}；"
                 f"目标 {exp_kept}）")
        if len(keep) != exp_kept:
            log.info(f"[警告] 去重后数量与规则期望（{exp_kept}）不一致，请检查上游表版本")

        dedup = qc[["cid", "ik_upstream", "CNN_VS", "dup_rank"]].copy()
        dedup["group_size"] = dedup.ik_upstream.map(qc.ik_upstream.value_counts())
        dedup["kept"] = (dedup.dup_rank == 0).astype(int)
        write_csv(ROOT / "01_pose_qc/tables/pose_dedup_map.csv", dedup)

        keep = keep.sort_values("cid").reset_index(drop=True)
        annot = keep[ANNOT_COLS].copy()
        annot["n_fragments"] = annot["n_fragments_prop"].fillna(annot["n_fragments_template"])
        annot["warn_or_reject_reason"] = ""       # 由 04 步填姿态闸门原因
        annot["pose_gate"] = "pending"
        write_csv(ROOT / "01_pose_qc/tables/pose_pool.csv", annot)

        # ---- 规范化 SDF（几何不变，仅统一属性块）
        wdir = ROOT / "01_pose_qc/outputs"
        wdir.mkdir(parents=True, exist_ok=True)
        lookup = annot.set_index("cid")
        out = open(wdir / "SOP-05_poses_normalized.sdf", "w", encoding="utf-8")
        written, missing = 0, []
        for block, props, title in cc.iter_sdf_blocks(ROOT / "00_inputs/raw/inputs/top5000_best_pose.sdf"):
            cid = props.get("cid", title)
            if cid not in lookup.index:
                continue
            r = lookup.loc[cid]
            # 去掉原有属性块，写统一属性块（数据部分以 M  END 收尾，保证块合法）
            core = cc.molblock_data_lines(block)
            out.write("\n".join(core) + "\n")
            for key, val in [("cid", cid), ("smiles_upstream", str(r["smiles_upstream"])),
                             ("ik_upstream", str(r["ik_upstream"])),
                             ("CNN_VS", f"{float(r['CNN_VS']):.6f}"),
                             ("n_fragments", str(r["n_fragments"])),
                             ("stereo_source", str(r["stereo_source"])),
                             ("protonation_note", str(r["protonation_note"]) or "none"),
                             ("qc_class", str(r["qc_class"])),
                             ("identity_class", str(r["identity_class"])),
                             ("pose_source", "top5000_best_pose.sdf（平台回传最优姿态）")]:
                out.write(f"> <{key}>\n{val}\n\n")
            out.write("$$$$\n")
            written += 1
        out.close()
        if written != len(lookup):
            log.info(f"[警告] 规范化 SDF 写入 {written} 条，池内 {len(lookup)} 条，差值需排查")

        summary = {
            "n_in": int(len(qc)), "n_kept": int(len(keep)), "n_dropped": int(len(dropped)),
            "dropped_cids": sorted(dropped.cid.tolist()),
            "dropped_detail": dropped[["cid", "ik_upstream", "CNN_VS"]].to_dict("records"),
            "n_written_sdf": int(written),
            "identity_class_counts_kept": {str(k): int(v) for k, v in
                                           keep.identity_class.value_counts().sort_index().items()},
            "n_multifragment_in_pool": int((keep.n_fragments_template > 1).sum()),
            "positive_expectation_met": bool(len(keep) == exp_kept),
        }
        write_json(ROOT / "01_pose_qc/reports/pose_normalize_summary.json", summary)
        log.info(f"规范化 SDF → 01_pose_qc/outputs/SOP-05_poses_normalized.sdf（{written} 条）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
