#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 5G：确定性重跑（全量产物 sha256 比对）。

哈希范围：包内全部「数值产物」= 各步骤 tables/、outputs/、reports/、figures/、herg/、
07_md_batch/ 与 06_deliverables/ 下的文件（排除 00_inputs/raw、00_inputs/extracted、logs、
determinism 自身与 *.pyc）。首次运行写入 baseline，之后每次比对；
不一致即视为不确定（抖动源只有显式随机源，应完全一致）。

写入语义（防歧义）：`--freeze` / `--refreeze` 写基线时 mode=baseline_written，
`determinism_verified=false`、`identical=null` —— 只写基线的运行**不能**算作确定性通过
（它与基线同源，不是独立重跑）；只有 mode=compared 且 n_diff=0 才是判据。
输出：05_robustness/determinism/sha256_baseline.json、sha256_current.json
      05_robustness/reports/determinism_compare.json
"""
from __future__ import annotations

import sys
from pathlib import Path

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.io_utils import StepLog, read_json, sha256_file, write_json

INCLUDE_DIRS = ["01_pose_qc", "02_druglikeness", "03_gates", "04_selection",
                "05_robustness", "06_deliverables", "07_md_batch", "00_inputs/manifests"]
EXCLUDE_PARTS = ("logs", "determinism", "__pycache__")
EXCLUDE_SUFFIX = (".pyc",)
# 排除两类「自指」文件，否则判据永远到不了不动点（会 True/False 振荡）：
#  ① 比对报告本身（记录上一轮差异数）；
#  ② 交付层的「状态载体」：报告 / 索引 / 验收对账 / MANIFEST 会把 G 臂结论写进自身内容，
#     而该结论又由包含它们的哈希范围决定 —— 属纯自指。
# 其余全部解析型产物（各步骤 tables/outputs/figures/reports + 输入清单）都在范围内。
# 状态载体由被比对产物唯一决定：范围一旦封闭，稳态下它们也逐字节稳定（实测：连续两轮
# 载体内容一致，第二轮起 n_diff=0）。
EXCLUDE_REL = {"05_robustness/reports/determinism_compare.json",
               "06_deliverables/SOP-05_determinism_compare.json",
               "06_deliverables/SOP-05_report.md",
               "06_deliverables/INDEX.md",
               "06_deliverables/MANIFEST.sha256",
               "06_deliverables/acceptance_report.json"}


def collect(root: Path) -> dict:
    files = {}
    for d in INCLUDE_DIRS:
        base = root / d
        if not base.exists():
            continue
        for f in sorted(base.rglob("*")):
            if not f.is_file() or f.suffix in EXCLUDE_SUFFIX:
                continue
            if any(part in EXCLUDE_PARTS for part in f.parts):
                continue
            rel = str(f.relative_to(root))
            if rel in EXCLUDE_REL:
                continue
            files[rel] = sha256_file(f)
    return files


def main() -> int:
    freeze = "--freeze" in sys.argv
    refreeze = "--refreeze" in sys.argv
    det = ROOT / "05_robustness/determinism"
    det.mkdir(parents=True, exist_ok=True)
    with StepLog(ROOT, "12_determinism") as log:
        cur = collect(ROOT)
        write_json(det / "sha256_current.json", {"n_files": len(cur), "files": cur})
        log.info(f"哈希范围：{len(cur)} 个解析型产物文件"
                 f"（排除 6 个自指文件：2 个比对报告 + 4 个交付层状态载体）")
        base_path = det / "sha256_baseline.json"
        if refreeze or freeze or not base_path.is_file():
            overwriting = base_path.is_file()
            write_json(base_path, {"n_files": len(cur), "files": cur})
            why = ("--refreeze 显式重建" if refreeze else
                   "--freeze 强制覆盖" if freeze else "首次运行")
            log.info(f"已写入确定性基线（{len(cur)} 个文件，{why}"
                     f"{'；原有基线被覆盖' if overwriting else ''}）——**本轮未做比对**，"
                     f"不得记为「确定性通过」")
            compare = {"mode": "baseline_written", "n_files": len(cur),
                       "baseline_sha256": sha256_file(base_path),
                       "determinism_verified": False, "identical": None, "n_diff": None,
                       "note": "本轮只写基线（确定性「未检验」而非「通过」）：写入基线的运行与基线"
                               "自身同源，不构成独立重跑证据；需再跑一轮 mode=compared 才能给出判据"}
        else:
            base = read_json(base_path)["files"]
            diff = sorted(k for k in set(base) | set(cur) if base.get(k) != cur.get(k))
            compare = {"mode": "compared", "n_files_current": len(cur),
                       "n_files_baseline": len(base), "n_diff": len(diff),
                       "baseline_sha256": sha256_file(base_path),
                       "determinism_verified": True,
                       "identical": bool(len(diff) == 0),
                       "diff": diff[:50],
                       "missing_in_current": sorted(set(base) - set(cur))[:20],
                       "new_in_current": sorted(set(cur) - set(base))[:20]}
            log.info(f"确定性比对：{len(cur)} 个文件，差异 {len(diff)} 个"
                     f"{'（完全一致）' if not diff else ''}")
            if diff:
                log.info(f"差异文件（前 10）：{diff[:10]}")
        write_json(ROOT / "05_robustness/reports/determinism_compare.json", compare)
    if compare["mode"] == "baseline_written":
        return 0
    return 0 if compare["identical"] else 5


if __name__ == "__main__":
    raise SystemExit(main())
