"""包内路径解析与隔离守卫。

设计原则：**任何脚本都不得引用包根之外的路径**。
所有路径由 `pack_root()` 向上查找 `REPRO_INFO.json` 得到，
再由 `paths()` 拼出标准子目录；`guard_inside()` 对最终路径做越界断言。
"""
from __future__ import annotations

import os
from pathlib import Path

MARKER = "REPRO_INFO.json"

# 标准子目录（相对包根）
STD_DIRS = [
    "00_inputs/raw/inputs", "00_inputs/raw/archives", "00_inputs/extracted",
    "00_inputs/manifests",
    "01_pose_qc/scripts", "01_pose_qc/tables", "01_pose_qc/outputs",
    "01_pose_qc/figures", "01_pose_qc/reports",
    "02_druglikeness/scripts", "02_druglikeness/tables", "02_druglikeness/figures",
    "02_druglikeness/reports",
    "03_gates/scripts", "03_gates/tables", "03_gates/herg", "03_gates/figures",
    "03_gates/reports",
    "04_selection/scripts", "04_selection/tables", "04_selection/figures",
    "04_selection/reports",
    "05_robustness/scripts", "05_robustness/tables", "05_robustness/figures",
    "05_robustness/reports", "05_robustness/determinism",
    "06_deliverables", "06_deliverables/figures",
    "config", "docs", "environment", "vendor", "tests", "logs",
]


def pack_root(start: str | os.PathLike | None = None) -> Path:
    """从 start（默认本文件）向上查找含 REPRO_INFO.json 的目录。"""
    p = Path(start or __file__).resolve()
    for cand in [p, *p.parents]:
        if (cand / MARKER).is_file():
            return cand
    raise RuntimeError(f"未找到 {MARKER}：从 {p} 向上无包根，脚本不可在包外运行")


def paths(root: Path | None = None, ensure: bool = False) -> dict[str, Path]:
    """返回标准子目录字典 {相对路径: 绝对 Path}；ensure=True 时创建。"""
    r = Path(root) if root is not None else pack_root()
    out = {d: r / d for d in STD_DIRS}
    if ensure:
        for d in out.values():
            d.mkdir(parents=True, exist_ok=True)
    return out


def pth(root: Path, rel: str, ensure_parent: bool = True) -> Path:
    """包内文件路径：越界即报错。"""
    f = (Path(root) / rel).resolve()
    guard_inside(root, f)
    if ensure_parent:
        f.parent.mkdir(parents=True, exist_ok=True)
    return f


def guard_inside(root: Path, path: Path) -> None:
    r = Path(root).resolve()
    p = Path(path).resolve()
    if r != p and r not in p.parents:
        raise RuntimeError(f"路径越界（包外引用被禁止）：{p} 不在 {r} 之下")


def rel(root: Path, path: Path) -> str:
    """包内相对路径（用于清单与报告）。"""
    return str(Path(path).resolve().relative_to(Path(root).resolve()))
