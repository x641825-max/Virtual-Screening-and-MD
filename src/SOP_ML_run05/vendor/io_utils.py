"""IO / 哈希 / 清单 / 日志 工具。所有写出均为确定性（键排序、固定缩进、无时间戳进正文）。"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------- 哈希
def _hash_file(path: str | Path, algo: str, chunk: int = 1 << 20) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def md5_file(path) -> str:
    return _hash_file(path, "md5")


def sha256_file(path) -> str:
    return _hash_file(path, "sha256")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def canonical_json(obj) -> str:
    """规范 JSON 文本：键排序、无多余空白、非 ASCII 原样保留。"""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def canonical_hash(obj) -> str:
    return sha256_text(canonical_json(obj))


# ---------------------------------------------------------------- 读写
def read_json(path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path, obj, indent: int = 2) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=indent, ensure_ascii=False, sort_keys=False)
        fh.write("\n")
    return p


def write_csv(path, df: pd.DataFrame, float_format: str | None = "%.6g") -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False, float_format=float_format)
    return p


def write_text(path, text: str) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    return p


def read_table(path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


# ---------------------------------------------------------------- 日志
class StepLog:
    """同时打印到 stdout 并落盘到 logs/<step>.log（带步骤横幅）。"""

    def __init__(self, root: Path, step: str):
        self.step = step
        self.path = Path(root) / "logs" / f"{step}.log"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(self.path, "w", encoding="utf-8")
        self.banner(f"=== {step} ===")

    def banner(self, msg: str) -> None:
        line = f"\n{'=' * 4} {msg} {'=' * 4}"
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()

    def info(self, msg: str) -> None:
        print(f"[{self.step}] {msg}", flush=True)
        self.fh.write(f"[{self.step}] {msg}\n")
        self.fh.flush()

    def close(self) -> None:
        self.fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


# ---------------------------------------------------------------- 清单
def manifest_lines(root: Path, files: list[Path], algo: str = "sha256") -> list[str]:
    lines = []
    for f in sorted(files, key=lambda x: str(x)):
        if Path(f).is_file():
            lines.append(f"{_hash_file(f, algo)}  {Path(f).resolve().relative_to(Path(root).resolve())}")
    return lines


def write_manifest(path, root: Path, files: list[Path], algo: str = "sha256") -> Path:
    return write_text(path, "\n".join(manifest_lines(root, files, algo)))


# ---------------------------------------------------------------- 数值
def as_bool_series(df: pd.DataFrame, col: str) -> pd.Series:
    return df[col].astype(bool)


def r6(x) -> float:
    """统一四舍五入到 6 位，用于 JSON 报告的确定性表述。"""
    if x is None:
        return None
    if isinstance(x, (int, np.integer)):
        return int(x)
    try:
        f = float(x)
    except (TypeError, ValueError):
        return x
    if np.isnan(f) or np.isinf(f):
        return None
    return round(f, 6)


def r4(x):
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return x
    if np.isnan(f) or np.isinf(f):
        return None
    return round(f, 4)


def counts_dict(series) -> dict:
    return {str(k): int(v) for k, v in series.value_counts().sort_index().items()}


def eprint(*a) -> None:
    print(*a, file=sys.stderr, flush=True)
