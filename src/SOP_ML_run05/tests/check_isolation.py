#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""静态隔离检查：确认包内脚本不引用包外路径。

规则（本工程的硬约束）：步骤脚本、vendor、tests 中不得出现
  - 绝对路径字面量（以 / 开头且不属于 sys.executable 之类的运行期拼接）
  - 对工作区其他目录（如 ../、Other/、A2A_work）的引用
允许：PurePosixPath 之外的相对文件名常量（由 pack_root 拼接）。
退出码：0 = 干净；6 = 发现越界引用。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCAN_DIRS = ["00_inputs", "01_pose_qc", "02_druglikeness", "03_gates", "04_selection",
             "05_robustness", "06_deliverables", "vendor", "tests"]
BANNED = [
    (re.compile(r"""["']/workspace"""), "引用 /workspace 绝对路径"),
    (re.compile(r"""["']/(home|Users|mnt|data|public|skills)/"""), "引用包外绝对路径"),
    (re.compile(r"""\.\./(?!\.\./)*[A-Za-z]"""), "使用 ../ 上溯引用"),
    (re.compile(r"\bOther/"), "引用 Other/ 目录"),
    (re.compile(r"A2A_work|SOP03_run04_v6|A2AR_Screening_Workflow"), "引用既有工程目录"),
]
ALLOW_LINE = re.compile(r"sys\.path\.insert|# noqa: isolation")


def main() -> int:
    problems = []
    files = []
    for d in SCAN_DIRS:
        for p in (ROOT / d).rglob("*.py"):
            rel = p.relative_to(ROOT)
            # 解包出来的上游归档内容（00_inputs/raw、00_inputs/extracted）不是本包代码：
            # 它们逐字节受输入清单哈希约束，不可修改，故不参与「包内脚本禁止包外引用」扫描，
            # 但其存在性/哈希由 01_verify_inputs 与 tests/verify_pack.py 负责校验。
            if rel.parts[:2] in (("00_inputs", "raw"), ("00_inputs", "extracted")):
                continue
            files.append(p)
    files += [ROOT / "run_all.py"]
    self_path = Path(__file__).resolve()
    for f in sorted(set(files)):
        if f.resolve() == self_path:
            continue      # 检查器自身含违例模式字面量，不参与扫描
        for i, line in enumerate(f.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if ALLOW_LINE.search(line):
                continue
            for pat, why in BANNED:
                if pat.search(line):
                    problems.append((f.relative_to(ROOT), i, why, line.strip()[:110]))
    print(f"扫描 {len(set(files))} 个包内 Python 文件（解包的上游归档脚本已排除，见注释）；"
          f"发现越界引用 {len(problems)} 处")
    for p in problems:
        print(f"  [{p[0]}:{p[1]}] {p[2]} → {p[3]}")
    if problems:
        print("\n结论：**不通过**（包外引用被禁止）")
        return 6
    print("结论：通过 —— 全部路径均在包根（REPRO_INFO.json 所在目录）之下解析")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
