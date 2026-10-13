#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成包索引 INDEX_PACK.md 与总哈希清单 PACK_MANIFEST.sha256。

INDEX_PACK.md：目录树 + 每个子目录的用途 + 关键文件说明（自动生成，不手写）。
PACK_MANIFEST.sha256：包内全部文件（排除大体积输入/解包产物、日志、__pycache__）
的 sha256，用于跨机器校验一致性。
退出码：0。
"""
from __future__ import annotations

import sys
from pathlib import Path

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.io_utils import manifest_lines, read_json, sha256_file, write_text

PURPOSE = {
    "00_inputs": "数据层：5 个原始输入、3 个归档、白名单解包产物、哈希清单",
    "00_inputs/raw/inputs": "原始输入（姿态 SDF、全库排名+pKa/cLogD 表）",
    "00_inputs/raw/archives": "三个归档原件（run01 重对接+富集、筛选套件、重对接基准）",
    "00_inputs/extracted": "白名单解包产物（受体、5,098 带标签库、ChEMBL 独立集、3EML 基准等）",
    "00_inputs/manifests": "输入固化清单、解包清单、规则冻结与脚本哈希登记",
    "01_pose_qc": "步骤 1：姿态属性一致性、结构身份分级、去重、标注、姿态闸门与接触指纹",
    "02_druglikeness": "步骤 2：RDKit 重算描述符/PAINS/Brenk/QED/CNS MPO/化学型/骨架/指纹 + 漏斗",
    "03_gates": "步骤 3：anchor、maxTc、hERG 接口、候选主表、门叠加、已知物核验",
    "03_gates/herg": "hERG 外部接口：补算清单 + 接口说明书（本包不内置 hERG 预测）",
    "03_gates/herg/cross_model": "hERG 第二模型交叉验证（分析步 08b，非判据）",
    "03_gates/admet": "ADMET-AI 注释层（分析步 10b/10c，非判据）：冻结 f3 库表、池级注释、"
                      "hERG 阈值锚定与保留曲线、MD 前置复核清单",
    "04_selection": "步骤 4：6 组名单矩阵、决策卡、约束触发率、重叠热图",
    "05_robustness": "步骤 5：A–F 稳健性臂、Butina 聚类、G 臂确定性哈希基线",
    "06_deliverables": "交付层：名单矩阵、对照与决策卡、姿态子集、验收对账、报告、MANIFEST",
    "config": "frozen_rule.json（阈值唯一来源）、run_order.json（步序与选项）",
    "vendor": "共享模块：路径守卫、IO/哈希、化学核心、CNS MPO、选序、主表、绘图",
    "tests": "verify_pack.py（自检）、check_isolation.py（包外引用静态扫描）、make_index.py",
    "docs": "方法与来源、路径映射、核验说明",
    "environment": "环境快照与依赖版本",
    "logs": "运行日志（每步一份 + run_all.log）",
}
EXCLUDE_DIRS = ("00_inputs/raw", "00_inputs/extracted", "logs", "__pycache__")


def main() -> int:
    # 排除清单自身：否则第 N 轮记录的是第 N-1 轮的自哈希，每次都变、永不自洽
    self_name = "PACK_MANIFEST.sha256"
    files = [p for p in ROOT.rglob("*")
             if p.is_file() and not any(str(p.relative_to(ROOT)).startswith(d)
                                        for d in EXCLUDE_DIRS)
             and p.name != self_name
             and p.suffix != ".pyc"]
    write_text(ROOT / "PACK_MANIFEST.sha256", "\n".join(manifest_lines(ROOT, files, "sha256")))

    pkg = read_json(ROOT / "REPRO_INFO.json").get("package", "pack")
    lines = [f"# {pkg} 工程包索引（自动生成）", "",
             f"- 包根：`{pkg}/`（标记文件 `REPRO_INFO.json`）",
             f"- 清单文件数（排除 raw/extracted/logs）：{len(files)}",
             f"- 总数（含解包产物）见 `00_inputs/manifests/EXTRACTED_MANIFEST.md5`", "",
             "## 目录用途", "", "| 目录 | 用途 |", "|---|---|"]
    for d, why in PURPOSE.items():
        lines.append(f"| `{d}/` | {why} |")
    lines += ["", "## 顶层文件", "",
              "| 文件 | 用途 |", "|---|---|",
              "| `run_all.py` | 一键复现入口（`--list` 看步序；`--freeze` 固化确定性基线） |",
              "| `README.md` | 使用与复现说明、口径与局限 |",
              "| `REPRO_INFO.json` | 包标记与元信息（版本、日期、环境、验收口径） |",
              "| `PACK_MANIFEST.sha256` | 本清单（包内文件哈希） |",
              "| `config/frozen_rule.json` | 全部阈值与算法口径（唯一来源） |",
              "| `config/run_order.json` | 步序、脚本、产物、选项 |", "",
              "## 产物速查", "",
              "| 需求 | 文件 |", "|---|---|",
              "| 最终 50 名单 | `06_deliverables/SOP-05_MD_shortlist_main.csv` |",
              "| 逐分子证据卡 | `06_deliverables/SOP-05_decision_cards_main.csv` |",
              "| 起始姿态（MD 用） | `06_deliverables/SOP-05_poses_main_list.sdf` |",
              "| 验收对账 | `06_deliverables/acceptance_report.json` |",
              "| 最终报告 | `06_deliverables/SOP-05_report.md` |",
              "| 下游 MD 计划 | `06_deliverables/SOP-05_downstream_MD_plan.md` |",
              "| hERG 待补算清单 | `03_gates/herg/SOP-05_herg_input_union.csv` |",
              "| 稳健性结果 | `05_robustness/reports/robustness.json` |",
              "| ADMET 注释（主名单） | `06_deliverables/SOP-05_admet_annotation_main.csv` |",
              "| MD 前置复核清单 | `06_deliverables/SOP-05_admet_md_precheck.csv` |",
              "| hERG 阈值锚定 SOP | `06_deliverables/SOP-05_HERG_THRESHOLD_SOP.md` |",
              "| hERG 保留曲线 | `06_deliverables/SOP-05_herg_retention_curve.csv` |", ""]
    write_text(ROOT / "INDEX_PACK.md", "\n".join(lines))
    print(f"生成 INDEX_PACK.md 与 PACK_MANIFEST.sha256（{len(files)} 个文件）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
