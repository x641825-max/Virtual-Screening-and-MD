#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 00b（换数据版专用，非判据步）：登记「全库打分/排名表」的换版。

它做四件事，顺序固定：
  ① 契约校验：调用 tests/check_swap_csv.py 对新表做结构与取值校验（FAIL 则中止）；
  ② 落位：把新表复制到固定路径 00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv
     （文件名不可改：11 个脚本按字面量读取它）；
  ③ 重登记：更新 config/frozen_rule.json → inputs.<文件名> 的 md5/sha256，并把
     rules_sha256 置为 "PENDING"，交由下一次 01_verify_inputs 重新冻结
     （阈值与门一字不改；变的只是「输入登记」）；
  ④ 重建 00_inputs/INPUTS_MD5.txt（从 rules 登记项生成，便于 md5sum 快速核对）。

默认 dry-run，只有 --confirm 才落盘；只写 ① 指定的三处，不动其它文件。
退出码：0 = 完成（或 dry-run 完成）；7 = 契约校验 FAIL；8 = 文件缺失；2 = 参数/前置问题。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

NAME = "CNN_VS_ranking_49963_with_pKa_cLogD74.csv"
CANON = ROOT / "00_inputs/raw/inputs" / NAME
RULE = ROOT / "config/frozen_rule.json"
CHECKER = ROOT / "tests/check_swap_csv.py"
MD5_LIST = ROOT / "00_inputs/INPUTS_MD5.txt"
REPORT = ROOT / "00_inputs/manifests/swap_check.json"


def md5_file(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description="登记排名表换版（默认 dry-run）")
    ap.add_argument("--file", required=True, help="新表路径（可不在包内，会复制到固定位置）")
    ap.add_argument("--role", default=None, help="覆盖登记里的 role 说明（默认沿用原值）")
    ap.add_argument("--note", default="", help="换表原因/出处，写入登记项 swap_note")
    ap.add_argument("--skip-check", action="store_true", help="跳过契约校验（不推荐）")
    ap.add_argument("--force", action="store_true", help="契约校验 FAIL 也继续（须配合 --note 说明）")
    ap.add_argument("--confirm", action="store_true", help="真正落盘（默认只预览）")
    a = ap.parse_args()

    src = Path(a.file)
    if not src.is_file():
        print(f"[ERR] 找不到新表：{src}")
        return 8
    if a.force and not a.note:
        print("[ERR] --force 必须配合 --note 说明放行原因")
        return 2
    if not RULE.is_file():
        print(f"[ERR] 规则文件缺失：{RULE}")
        return 2
    rule = json.loads(RULE.read_text(encoding="utf-8"))
    if NAME not in rule.get("inputs", {}):
        print(f"[ERR] 规则里的 inputs 登记项不含 {NAME}")
        return 2

    print(f"[1/4] 契约校验：{src.name}（{src.stat().st_size/1048576:.2f} MB）")
    if a.skip_check:
        print("      --skip-check：跳过")
    else:
        cmd = [sys.executable, str(CHECKER), "--file", str(src)]
        # 2026-09-24：一并带上耦合件做覆盖校验（E1 行数下限取 SDF 记录数，
        # E12/E13 姿态覆盖与属性同源、E14 hERG 落门表、E15 ADMET 冻结表）
        for flag, rel in (("--sdf", "00_inputs/raw/inputs/top5000_best_pose.sdf"),
                          ("--herg", "03_gates/herg/herg_predictions.csv"),
                          ("--admet", "03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv")):
            if (ROOT / rel).is_file():
                cmd += [flag, str(ROOT / rel)]
        if a.confirm:
            cmd += ["--json", str(REPORT)]
        proc = subprocess.run(cmd, text=True)
        if proc.returncode != 0 and not a.force:
            print(f"      [中止] 校验未通过（exit={proc.returncode}）。修掉 FAIL 后重跑，"
                  f"确需放行请加 --force --note <原因>")
            return 7
        if proc.returncode != 0:
            print(f"      [放行] 校验 exit={proc.returncode}，按 --force 继续")

    new_md5, new_sha = md5_file(src), sha256_file(src)
    old = rule["inputs"][NAME]
    print(f"[2/4] 落位：{CANON.relative_to(ROOT)}")
    if src.resolve() == CANON.resolve():
        print("      源文件已在固定位置，无需复制")
    else:
        if CANON.is_file():
            print(f"      将覆盖现有同名文件（现 md5 {md5_file(CANON)}）——如需留档请先自行备份到包外")
        print(f"      复制 {src} → {CANON}")

    print(f"[3/4] 重登记：inputs.{NAME}")
    print(f"      md5    {old.get('md5')} → {new_md5}")
    print(f"      sha256 {str(old.get('sha256'))[:16]}… → {new_sha[:16]}…")
    print('      rules_sha256 → "PENDING"（下次 01_verify_inputs 重新冻结；阈值与门不动）')
    if a.note:
        print(f"      swap_note: {a.note}")

    print("[4/4] 重建 00_inputs/INPUTS_MD5.txt（从登记项生成）")

    if not a.confirm:
        print("\n[dry-run] 未落盘。确认无误后加 --confirm 执行。")
        return 0

    if src.resolve() != CANON.resolve():
        CANON.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, CANON)
    entry = {"md5": new_md5, "sha256": new_sha,
             "role": a.role or old.get("role", "")}
    if a.note:
        entry["swap_note"] = a.note
        entry["swapped_on"] = __import__("datetime").date.today().isoformat()
    rule["inputs"][NAME] = entry
    rule["rules_sha256"] = "PENDING"
    RULE.write_text(json.dumps(rule, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = []
    for name, spec in sorted(rule["inputs"].items()):
        sub = "raw/archives" if name.endswith(".zip") else "raw/inputs"
        lines.append(f"{spec['md5']}  {sub}/{name}")
    MD5_LIST.write_text("\n".join(lines) + "\n", encoding="utf-8")

    missing = [n for n in rule["inputs"] if not (ROOT / "00_inputs/raw/inputs" / n).is_file()
               and not (ROOT / "00_inputs/raw/archives" / n).is_file()]
    print("\n[完成] 已登记。")
    if missing:
        print(f"[提示] 以下登记输入当前不在包内，跑 01_verify_inputs 会显式失败（exit 4），"
              f"需先按 RESTORE_LIST.md 恢复：{missing}")
    print("下一步：python3 tests/make_index.py && python3 run_all.py（首轮顺带重新冻结 rules_sha256）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
