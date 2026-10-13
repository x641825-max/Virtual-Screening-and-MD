#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SOP_5 一键复现入口。

用法：
  python3 run_all.py --list            列出步序（不执行）
  python3 run_all.py                   全链执行（00→14）+ 末尾自检
  python3 run_all.py --steps 00_unpack,01_verify_inputs
  python3 run_all.py --from 10_select_lists
  python3 run_all.py --freeze          把本次产物固化为确定性基线（重建判据基准）
  python3 run_all.py --refreeze       同 --freeze，但语义显式：**重建基线**，本轮不构成确定性通过
  python3 run_all.py --skip-tests      跳过末尾自检（verify_pack / check_isolation / make_index）

步序末步为 12_determinism（G 臂）：它覆盖 13/14 交付层产物做 sha256 比对，
首轮写基线（mode=baseline_written，determinism_verified=false），第二次起给出判据。
起链前按 config/run_order.json 的 outputs 预建产物父目录：即使整目录被删
（而非只删文件），重跑也能照常写出，不必先手工 mkdir。
单实例锁：同一包目录同时只允许一条链运行（防止两条链交叉写同一批产物，
使「确定性」结论无效）；锁文件在 logs/.run_all.lock，退出时自动释放。

约定：所有步骤脚本以子进程运行，环境变量 PYTHONHASHSEED=0（与规则集一致），
工作目录固定为包根，任何失败立即中止并回传非零退出码。
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ORDER = json.loads((ROOT / "config/run_order.json").read_text(encoding="utf-8"))


def run_step(step: dict, env: dict, logf, extra_argv: list | None = None) -> int:
    script = ROOT / step["script"]
    if not script.is_file():
        print(f"[FAIL] 缺少脚本 {step['script']}")
        return 2
    t0 = time.time()
    line = f"\n>>> [{step['id']}] {step['desc']}\n    {step['script']}"
    print(line, flush=True)
    logf.write(line + "\n")
    logf.flush()
    p = subprocess.run([sys.executable, str(script)] + (extra_argv or []),
                       cwd=str(ROOT), env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    dt = time.time() - t0
    tail = p.stdout.strip().splitlines()[-6:] if p.stdout else []
    for t in tail:
        print("    " + t)
    print(f"    exit={p.returncode}  用时 {dt:.1f}s", flush=True)
    logf.write(p.stdout + f"\n--- exit={p.returncode} 用时 {dt:.1f}s ---\n")
    logf.flush()
    return p.returncode


def ensure_output_dirs(steps: list) -> list:
    """按 run_order.json 的 outputs 预建产物父目录。

    产物目录若被整体删除（而非只删文件），步骤脚本用 Path.write_text 写文件时
    不会自建父目录，会在中段以 FileNotFoundError 中止。这里在起链前把本轮
    所有声明的产物目录建好，使「删目录 → 重跑」与「删文件 → 重跑」行为一致。
    只建目录、不建文件，因此不影响任何产物哈希与确定性判据。
    """
    created = []
    for step in steps:
        for out in step.get("outputs", []):
            d = (ROOT / out).parent
            if not d.is_dir():
                d.mkdir(parents=True, exist_ok=True)
                created.append(str(d.relative_to(ROOT)) + "/")
    return created


def acquire_lock():
    """单实例锁：同一包目录不允许两条链并行（并行链会让确定性判据失效）。"""
    lock = ROOT / "logs/.run_all.lock"
    lock.parent.mkdir(exist_ok=True)
    if lock.is_file():
        pid = lock.read_text(encoding="utf-8").strip()
        alive = False
        if pid.isdigit():
            try:
                os.kill(int(pid), 0)
                alive = True
            except OSError:
                alive = False
        if alive:
            print(f"[ABORT] 已有链在运行（pid {pid}，锁 {lock}）；"
                  f"同一包内只允许单实例，以免两条链交叉写产物。")
            return None
        print(f"[提示] 发现残留锁（pid {pid} 已不存在），接管。")
    fd = os.open(lock, os.O_CREAT | os.O_TRUNC | os.O_WRONLY)
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    return lock


def main() -> int:
    args = sys.argv[1:]
    steps = ORDER["steps"]
    if "--list" in args:
        print(f"{ORDER['package']} · {ORDER['title']}\n")
        for i, s in enumerate(steps, 1):
            print(f"{i:2d}. {s['id']:20s} {s['script']:48s} {s['desc']}")
        print("\n重复上述执行即完成全链复现；末尾自动运行 tests/verify_pack.py。")
        return 0

    if "--steps" in args:
        want = set(args[args.index("--steps") + 1].split(","))
        steps = [s for s in steps if s["id"] in want]
    if "--from" in args:
        start = args[args.index("--from") + 1]
        ids = [s["id"] for s in steps]
        if start not in ids:
            print(f"[FAIL] 未知步骤 {start}；用 --list 查看")
            return 2
        steps = steps[ids.index(start):]

    env = dict(os.environ)
    env["PYTHONHASHSEED"] = "0"
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "1")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    (ROOT / "logs").mkdir(exist_ok=True)
    lock = acquire_lock()
    if lock is None:
        return 3
    try:
        return _run(args, steps, env)
    finally:
        try:
            lock.unlink()
        except OSError:
            pass


def _run(args: list, steps: list, env: dict) -> int:
    tmp = ROOT / "logs/run_all.log"
    _freezing = ("--freeze" in args) or ("--refreeze" in args)
    opts = ["--freeze"] if ("--freeze" in args) else (["--refreeze"] if "--refreeze" in args else [])

    # 每轮滚动：上一轮整体移入 logs/history/，本轮从空文件开始写（避免多轮数字混装）
    hist = ROOT / "logs/history"
    hist.mkdir(exist_ok=True)
    if tmp.is_file() and tmp.stat().st_size > 0:
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        tmp.rename(hist / f"run_all_{stamp}.log")

    print(f"{ORDER['package']} · {ORDER['title']}")
    print(f"包根：{ROOT}")
    print(f"步数：{len(steps)}；PYTHONHASHSEED={env['PYTHONHASHSEED']}")
    t0 = time.time()
    with open(tmp, "w", encoding="utf-8") as logf:
        rule_file = ROOT / "config/frozen_rule.json"
        try:
            rule_sha = json.loads(rule_file.read_text(encoding="utf-8")).get("rules_sha256", "NA")
        except Exception:                                     # noqa: BLE001
            rule_sha = "NA"
        logf.write(f"# 运行开始 {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}"
                   f"｜步骤 {','.join(s['id'] for s in steps)}"
                   f"｜rules_sha256={rule_sha}"
                   f"｜argv={' '.join(args) if args else '(默认)'}\n")
        # 预建本轮声明产物的父目录（容忍「整目录被删」的重跑场景）
        created = ensure_output_dirs(steps)
        if created:
            note = f"[预建] 本轮补齐 {len(created)} 个产物目录：{', '.join(sorted(set(created)))}\n"
            print("    " + note.strip())
            logf.write(note)
            logf.flush()
        for s in steps:
            extra_argv = (["--freeze"] if "--freeze" in args else
                          ["--refreeze"] if "--refreeze" in args else []) \
                if s["id"] == "12_determinism" else []
            rc = run_step(s, env, logf, extra_argv)
            if rc != 0:
                print(f"\n[ABORT] 步骤 {s['id']} 失败（exit={rc}），已中止。"
                      f"日志：logs/{s['id']}.log、logs/run_all.log")
                return rc
        if "--skip-tests" not in args:
            # 12_determinism 已是步序末步（覆盖 13/14 产物），此处不重复执行
            # 顺序：先重建包级清单（INDEX_PACK/PACK_MANIFEST 排除自身，供 [12] 断言），
            # 再验收、最后隔离性扫描
            for extra in [["tests/make_index.py"],
                          ["tests/verify_pack.py", *opts],
                          ["tests/check_isolation.py"]]:
                cmd = [sys.executable, str(ROOT / extra[0]), *extra[1:]]
                print(f"\n>>> 自检 {' '.join(extra)}", flush=True)
                p = subprocess.run(cmd, cwd=str(ROOT), env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
                for t in (p.stdout.strip().splitlines()[-8:] if p.stdout else []):
                    print("    " + t)
                logf.write(p.stdout + "\n")
                logf.flush()
                if p.returncode != 0:
                    print(f"\n[ABORT] 自检失败（exit={p.returncode}）")
                    return p.returncode
    print(f"\nALL DONE（{len(steps)} 步 + 自检，用时 {time.time() - t0:.1f}s）")
    print("产物：06_deliverables/（名单矩阵、决策卡、验收对账、报告、MANIFEST.sha256）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
