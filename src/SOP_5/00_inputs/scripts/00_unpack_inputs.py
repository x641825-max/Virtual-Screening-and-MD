#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 00：按白名单解包归档 → 00_inputs/extracted/

只解出本工程真正需要的成员（GF 白名单），不解压整个归档；
大体积文件（5,098 集姿态 poses.sdf）受 run_order.options.extract_large 控制。
输出：00_inputs/manifests/EXTRACTED_MANIFEST.md5、extracted_manifest.json
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.io_utils import StepLog, md5_file, read_json, sha256_file, write_json, write_text

ARCH = ROOT / "00_inputs/raw/archives"
DEST = ROOT / "00_inputs/extracted"
MAN = ROOT / "00_inputs/manifests"

WHITELIST = {
    "run01_redock_enrichment.zip": {
        "prefix": "VS_run01_archive_redock_enrichment/",
        "dest": "run01",
        "members": [
            "README.md", "run_config.json", "MANIFEST.md5",
            "01_redocking/receptor/4EIY_protein.pdb",
            "01_redocking/receptor/4EIY_ligand.pdb",
            "01_redocking/receptor/ZMA.cif",
            "01_redocking/receptor/receptor_meta.json",
            "01_redocking/receptor/docking_box_rule_23.042.json",
            "01_redocking/results/redock_summary.json",
            "01_redocking/results/redock_validation.json",
            "02_enrichment_5098/scores/gnina_ranking.csv",
            "02_enrichment_5098/scores/gnina_all_poses.csv",
            "02_enrichment_5098/library/quick_5100_library.csv",
            "02_enrichment_5098/library/library_meta.json",
            "02_enrichment_5098/triage/contact_fingerprint_top200.csv",
            "02_enrichment_5098/triage/top_hits.csv",
            "02_enrichment_5098/triage/triage_summary.json",
            "02_enrichment_5098/triage/scaffold_clusters.csv",
            "02_enrichment_5098/enrichment/enrichment_multi.json",
            "02_enrichment_5098/enrichment/funnel_backtest.json",
            "02_enrichment_5098/scripts/contact_fp.py",
        ],
        "large": ["02_enrichment_5098/poses_raw/poses.sdf"],
    },
    "screening_kit.zip": {
        "prefix": "A2A_screening_kit/",
        "dest": "screening_kit",
        "members": [
            "README.md", "SOP_report.md",
            "01_external_screen/README_step1.md",
            "01_external_screen/requirements.txt",
            "01_external_screen/example_input.csv",
            "02_return_to_platform/README_step2.md",
            "02_return_to_platform/docking_protocol.md",
            "02_return_to_platform/receptor/4EIY_protein.pdb",
            "03_model_card/model_card.md",
            "03_model_card/training_summary.json",
            "03_model_card/dataset_summary.json",
            "03_model_card/portable_check.json",
            "05_screen_test/README_test.md",
            "07_independent_validation/README_independent.md",
            "07_independent_validation/independent_a2a_chembl.csv",
            "08_split_for_retrain/README_split.md",
            "10_full_pipeline_test/README_full_test.md",
            "12_docking_validation/README_docking_validation.md",
        ],
        "large": [],
    },
    "redocking.zip": {
        "prefix": "A2A_redocking/",
        "dest": "redocking",
        "members": [
            "00_3EML_structure/3EML_crystal_zma.pdb",
            "00_3EML_structure/redocking_report_summary.md",
            "00_3EML_structure/3EML_redock.log",
            "SOP-02/README.md",
            "SOP-02/SOP-02_final_summary_report.md",
            "SOP-02/01_raw_data/actives_final.ism",
            "SOP-02/01_raw_data/decoys_final.ism",
            "SOP-02/02_cleaned_benchmark/actives_cleaned.csv",
            "SOP-02/02_cleaned_benchmark/decoys_cleaned.csv",
            "SOP-02/02_cleaned_benchmark/benchmark_sample_5100_eval.csv",
            "SOP-02/02_cleaned_benchmark/property_matching_validation.csv",
            "SOP-02/04_docking_results/docked_molecules.csv",
            "SOP-02/04_docking_results/undocked_molecules.csv",
            "SOP-02/05_evaluation_metrics/docking_scores_ranked.csv",
            "SOP-02/05_evaluation_metrics/qc_gate2_results.csv",
        ],
        "large": [],
    },
}


def main() -> int:
    opts = read_json(ROOT / "config/run_order.json")["options"]
    extract_large = bool(opts.get("extract_large", True))
    with StepLog(ROOT, "00_unpack") as log:
        log.info(f"extract_large={extract_large}")
        recs, missing = [], []
        for zname, spec in WHITELIST.items():
            zpath = ARCH / zname
            if not zpath.is_file():
                log.info(f"[缺失] 归档不存在：{zpath.relative_to(ROOT)}")
                return 2
            names = list(spec["members"]) + (list(spec["large"]) if extract_large else [])
            with zipfile.ZipFile(zpath) as zf:
                have = set(zf.namelist())
                for rel in names:
                    member = spec["prefix"] + rel
                    if member not in have:
                        missing.append(f"{zname}::{member}")
                        continue
                    out = DEST / spec["dest"] / rel
                    out.parent.mkdir(parents=True, exist_ok=True)
                    with zf.open(member) as src, open(out, "wb") as dst:
                        dst.write(src.read())
                    recs.append({"archive": zname, "member": member,
                                 "extracted": str(out.relative_to(ROOT)),
                                 "size": out.stat().st_size,
                                 "md5": md5_file(out), "sha256": sha256_file(out)})
            log.info(f"[解包] {zname} → extracted/{spec['dest']}/（{len(names)} 个成员）")
        if missing:
            log.info(f"[致命] 白名单成员缺失 {len(missing)} 个：{missing[:5]}")
            return 3

        write_json(MAN / "extracted_manifest.json",
                   {"extract_large": extract_large, "n_files": len(recs),
                    "total_bytes": sum(r["size"] for r in recs),
                    "files": sorted(recs, key=lambda r: r["extracted"])})
        write_text(MAN / "EXTRACTED_MANIFEST.md5",
                   "\n".join(f"{r['md5']}  {r['extracted']}"
                             for r in sorted(recs, key=lambda r: r["extracted"]))
                   + "\n\n# 源归档成员 sha256\n"
                   + "\n".join(f"{r['sha256']}  {r['archive']}::{r['member']}"
                               for r in sorted(recs, key=lambda r: r["extracted"])))
        log.info(f"解包完成：{len(recs)} 文件，"
                 f"{sum(r['size'] for r in recs) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
