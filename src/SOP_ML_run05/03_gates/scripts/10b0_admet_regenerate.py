#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ADMET-AI 预测表**再生成**脚本（手工执行、不入链，与 08b0 同一处理方式）。

用途：当需要重算 `03_gates/admet/raw/` 下的 ADMET-AI 预测表时（例如库更新、模型升级），
在包外环境调用 HPC「admet-ai」应用重跑，并把结果落到本脚本指定的输出路径；
随后更新 `config/frozen_rule.json → step3_gates.admet_layer` 中的 md5/sha256，
再整链重跑 `run_all.py`（注释层与阈值锚定会随新表重算）。

为什么不在链内跑：本包不打包 ADMET-AI / chemprop / torch 运行时与权重，且包内脚本
禁止引用包外路径与联网（见 tests/check_isolation.py）。因此预测层是**登记输入**，
包内复现的是「冻结预测值 → 标注列 / 阈值锚定 / 保留曲线 / 交叉核验」的确定性计算。

参考运行（本次冻结表的来源）：
  - 应用：HPC 目录应用 admet-ai；job_id f610251d9cdb473aa346a94c0e7c4331（2026-09-22）
  - 输入：去重后的 54,932 分子 SMILES 表（md5 0f84b5b5aaf7de43b62b85eac3f2ae31）
  - 原始输出：predictions.csv（sha256 a8c952d642445926bb4b1b323a1ba591eb4158d7a9af4cfe16063ed79d6a85e3）
  - 回填后冻结表：CNN_VS_ranking_49963_ADMET.csv
    （md5 e19ab1f7c5aafc3267d95f1aeb588c52 / sha256 5d693ea7d9d06e29375f7b73fc496911ce514444ef2e243e3b2feb4703a67476）

用法（在包外执行，路径自定）：
  python 10b0_admet_regenerate.py --check 03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv
  python 10b0_admet_regenerate.py --fill-input ranking.csv --input-out mol.csv
  python 10b0_admet_regenerate.py --record-meta out/meta.json --input mol.csv --pred pred.csv --table table.csv
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def _hash(path: Path, algo: str, chunk: int = 1 << 20) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def _desc(path: Path) -> dict:
    return {"file": path.name, "bytes": path.stat().st_size,
            "md5": _hash(path, "md5"), "sha256": _hash(path, "sha256")}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="ADMET-AI 再生成助手（手工执行、不入链）")
    ap.add_argument("--check", metavar="TABLE",
                    help="校验现有冻结表：行数/列数/哈希，并打印登记用信息")
    ap.add_argument("--fill-input", metavar="RANKING_CSV",
                    help="从筛选排名表抽 cid+smiles 生成 admet-ai 应用输入（smiles,cid）")
    ap.add_argument("--input-out", metavar="CSV", default="molecules_admet.csv")
    ap.add_argument("--record-meta", metavar="JSON",
                    help="记录输入/原始输出/回填表的哈希与规模，供登记进 frozen_rule")
    ap.add_argument("--input", metavar="CSV", help="admet-ai 应用的输入表")
    ap.add_argument("--pred", metavar="CSV", help="admet-ai 应用的原始输出 predictions.csv")
    ap.add_argument("--table", metavar="CSV", help="回填后的库级 ADMET 表")
    a = ap.parse_args(argv)

    if not any([a.check, a.fill_input, a.record_meta]):
        ap.print_help()
        return 2

    if a.fill_input:
        import pandas as pd
        src = Path(a.fill_input)
        df = pd.read_csv(src, usecols=lambda c: c in ("cid", "smiles"), low_memory=False)
        miss = [c for c in ("cid", "smiles") if c not in df.columns]
        if miss:
            print(f"[FAIL] 排名表缺少列 {miss}")
            return 3
        if df.smiles.duplicated().any():
            print(f"[WARN] 存在重复 SMILES {int(df.smiles.duplicated().sum())} 条；"
                  f"admet-ai 应用按唯一分子推理，重复行可在回填阶段展开")
        df.to_csv(a.input_out, index=False)
        print(f"已写出应用输入 {a.input_out}（{len(df)} 行）→ 交给 HPC 应用 admet-ai")
        print(f"  {_desc(Path(a.input_out))}")

    if a.check:
        import pandas as pd
        p = Path(a.check)
        head = pd.read_csv(p, nrows=5, low_memory=False)
        n_cols = len(head.columns)
        n_admet = sum(1 for c in head.columns if c.startswith("ADMET_"))
        n_rows = sum(1 for _ in open(p, encoding="utf-8", errors="replace")) - 1
        print(json.dumps({"rows": n_rows, "columns": n_cols, "admet_columns": n_admet,
                          **_desc(p)}, ensure_ascii=False, indent=2))
        print("把上面的 md5/sha256 与 rows/columns 登记到 "
              "config/frozen_rule.json → step3_gates.admet_layer")

    if a.record_meta:
        need = {"input": a.input, "pred": a.pred, "table": a.table}
        miss = [k for k, v in need.items() if not v]
        if miss:
            print(f"[FAIL] 记录元数据需要 --input/--pred/--table（缺 {miss}）")
            return 3
        meta = {k: _desc(Path(v)) for k, v in need.items()}
        out = Path(a.record_meta)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"已写出元数据 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
