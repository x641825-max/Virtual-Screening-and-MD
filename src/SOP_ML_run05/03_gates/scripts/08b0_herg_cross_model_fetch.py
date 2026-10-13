#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3b0（**手工执行、不入链**）：抓取第二模型交叉验证所需的 ChEMBL hERG 活性数据。

为什么单独出来：全链必须离线可复现，因此**联网抓取不放进 run_all**。本脚本只在需要
刷新数据时手工运行；产物（原始活性表 + 抓取元信息，含 md5/sha256 与抓取时间）随包固化，
链上的 `08b_herg_cross_model.py` 只读这份固化数据。

数据来源：ChEMBL REST API（EBI）
  靶点 CHEMBL240 = KCNH2 / hERG（Homo sapiens, SINGLE PROTEIN）
  过滤：standard_type ∈ {IC50, Ki}、standard_relation = '='、pchembl_value 非空

输出：
  03_gates/herg/cross_model/raw/chembl_herg_kcnh2_activities.csv
  03_gates/herg/cross_model/raw/fetch_meta.json
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
FIELDS = ("molecule_chembl_id,canonical_smiles,standard_type,standard_relation,"
          "standard_value,standard_units,pchembl_value,assay_chembl_id,potential_duplicate")
KEEP = ("molecule_chembl_id", "canonical_smiles", "standard_type", "standard_relation",
        "standard_value", "standard_units", "pchembl_value", "assay_chembl_id",
        "potential_duplicate")
OUTDIR = ROOT / "03_gates/herg/cross_model/raw"


def fetch(standard_type: str) -> tuple[list[dict], int, list[str]]:
    rows, offset, total, urls = [], 0, None, []
    while True:
        params = ["target_chembl_id=CHEMBL240", f"standard_type={standard_type}",
                  "standard_relation=%3D", "pchembl_value__isnull=false",
                  f"only={FIELDS}", "limit=1000", f"offset={offset}"]
        url = BASE + "?" + "&".join(params)
        urls.append(url)
        out = subprocess.run(["curl", "-s", "-m", "120", url],
                             capture_output=True, check=True).stdout.decode("utf-8")
        d = json.loads(out)
        total = d["page_meta"]["total_count"]
        acts = d["activities"]
        rows.extend({k: a.get(k) for k in KEEP} for a in acts)
        offset += len(acts)
        if not acts or offset >= total:
            break
    return rows, int(total or 0), urls


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows, totals, urls = [], {}, []
    for st in ("IC50", "Ki"):
        r, n, u = fetch(st)
        rows.extend(r)
        totals[st] = n
        urls.extend(u)
        print(f"  {st}：{len(r)} 行（API 报告 {n}）")
    if not rows:
        print("[ERROR] 未抓到任何记录，未写文件")
        return 1

    out = OUTDIR / "chembl_herg_kcnh2_activities.csv"
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(KEEP))
        w.writeheader()
        w.writerows(rows)
    raw = out.read_bytes()
    meta = {
        "source": "ChEMBL REST API（EBI）",
        "url_base": BASE,
        "query": {"target_chembl_id": "CHEMBL240 (KCNH2 / hERG, Homo sapiens)",
                  "standard_type": ["IC50", "Ki"], "standard_relation": "=",
                  "pchembl_value": "not null"},
        "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_rows": len(rows), "api_total_count": totals, "n_requests": len(urls),
        "first_url": urls[0],
        "retrieved_by": "ATLAS 手工执行（联网步骤不入链）",
        "file": out.name,
        "md5": hashlib.md5(raw).hexdigest(),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }
    (OUTDIR / "fetch_meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"写出 {out.relative_to(ROOT)}（{len(rows)} 行，md5 {meta['md5']}，"
          f"耗时 {time.time() - t0:.1f}s）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
