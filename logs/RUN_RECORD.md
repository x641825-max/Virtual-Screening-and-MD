# RUN_RECORD · 交付前复跑与核验记录

日期：**2026-10-01**（UTC）｜环境：Python 3.12.2 / numpy 2.4.6 / pandas 2.3.3 / scipy 1.15.3 /
scikit-learn 1.9.1 / matplotlib 3.11.1 / rdkit 2025.09.5 / openpyxl 3.1.5｜纯 CPU，全程离线

---

## 1. 复跑范围

默认一键入口只复跑**「MD 候选选择 → 联立合并 → 出结果文件」**这一段：
该段的全部输入（名单表、资格池表、ECFP4 指纹 npz、可合成性表、冻结规则）随包发布，
不需要外网与上游大体积文件（上游原始数据见 `data/DATA_SOURCES.md`）。

| # | 命令 | 产物 |
|---|---|---|
| 1 | `python3 src/SOP_5/07_md_batch/scripts/15_md_batch_select.py` | 目录链 20 席批次 |
| 2 | `python3 src/SOP_ML_run05/07_md_batch/scripts/15_md_batch_select.py` | 生成链 19 行批次 |
| 3 | `python3 src/SOP_5/10_final_md/scripts/17_final_md_merge.py` | 合并层 22 行名册 + 最优性证书 |
| 4 | `python3 src/SOP_5/10_final_md/scripts/18_final_md_figures.py` | 22 张结构卡图 + 批次构成图 |
| 5 | `python3 predict.py --source repro --out results/reproduced` | 复跑清单（61 行 × 31 列） |
| 6 | `python3 predict.py --verify results/reproduced/results.csv results/results.csv` | 一致性核验 |
| 7 | `python3 tests/check_isolation.py`（两个包各一次） | 越界路径静态扫描 |
| 8 | `python3 notebooks/md_candidates_walkthrough.ipynb` 的全部代码单元 | Notebook 端到端执行 |

上述 1–7 即 `./run.sh` 的内容（第 8 项另行执行）。

---

## 2. 核验结果

| 核验项 | 结果 |
|---|---|
| 目录链选批复现 | `t* = 4`、跨度 `0.7337`、Σ CNN_VS `136.5031`；四阶段 `status=0`、`mip_gap=0`、重解一致 —— 与冻结值一致 |
| 生成链选批复现 | `t* = 3`、跨度 `3.5015`、Σ CNN_VS `119.7563`；四阶段 `status=0`、`mip_gap=0` —— 与冻结值一致 |
| 合并层复现 | 合并池 35、资格池 33、名册 22 行（目录 13 + 生成 7）、`t* = 5`、跨度 `3.5015`、Σ CNN_VS `150.4726`、敏感性重叠 13/20 —— 与冻结值一致 |
| 逐字节一致 | `md_batch20.csv`（两条链）与 `final_md_20.csv`、`final_md_excluded.csv` 的 sha256 与源工程**完全相同** |
| 候选清单核验 | `[PASS] 两份清单逐行逐字段一致（61 行 × 31 列）` |
| 越界路径扫描 | SOP_5：扫描 32 个 Python 文件，**越界引用 0 处**；SOP_ML_run05：扫描 34 个，**越界引用 0 处** |
| Notebook 执行 | 8 个代码单元串行执行 `exit 0`，末尾同样得到 `[PASS]` |

### 2.1 逐字节一致的判定方式

```bash
sha256sum <源工程>/<文件> <包内>/src/<同路径文件>
```

比对结果（源工程与包内）：

| 文件 | 结果 |
|---|---|
| `SOP_5/07_md_batch/tables/md_batch20.csv` | IDENTICAL |
| `SOP_ML_run05/07_md_batch/tables/md_batch20.csv` | IDENTICAL |
| `SOP_5/10_final_md/tables/final_md_20.csv` | IDENTICAL |
| `SOP_5/10_final_md/tables/final_md_excluded.csv` | IDENTICAL |
| `SOP_5/10_final_md/reports/final_md.json` | 内容不同（**仅因输入路径字段**：源工程记录绝对路径，本包已改为 `<包根>/` 相对占位符；所有 sha256 与读数一致） |

### 2.2 复跑日志

完整 stdout 见 `logs/run.log`（`./run.sh` 每次运行覆盖写入）。

---

## 3. 打包期对代码所做的改动（如实登记）

| 改动 | 原因 |
|---|---|
| `src/SOP_5/10_final_md/scripts/17_final_md_merge.py`：`ML_ROOT` 由硬编码的绝对路径（`<源工作区>/SOP_ML_run05`）改为 `SOPML_ROOT` 环境变量或包根同级目录 | 附件 5 §六要求「避免硬编码绝对路径」；同时使源工程布局与本包布局都能直接跑 |
| `src/SOP_5/10_final_md/scripts/18_final_md_figures.py`：同上一处硬编码路径改为 `ML_ROOT` | 同上 |
| 两处文档里的 `cd <源工作区>/...` 复现命令改为 `cd <包根>/src/...` | 同上 |
| `results/merged_final_md/final_md.json`、`src/SOP_5/10_final_md/reports/final_md.json`：输入绝对路径改为 `<包根>/` 占位符 | 避免暴露构建期宿主路径；输入 sha256 未改动 |
| `src/SOP_ML_run05/README.md` 重写、`REPRO_INFO.json` 的 `title`/`version` 更新 | 原文件是换表骨架期的陈旧 README（自称「骨架 v0.1.0、未跑链」），与本包实际产物不符 |
| 三套候选清单归一为 `results/results.csv` / `results.xlsx`（新增 `predict.py`） | 附件 5 §四要求标准化最终结果文件与必填字段 |

**未改动**：任何阈值、门、名单成员、规则集 `config/frozen_rule.json`、`PACK_MANIFEST.sha256`、
`REPRO_INFO.json` 的 `package` 字段与 `rule_set`、以及各阶段的中间产物内容。

> **2026-10-03 例外**：`data/manifests/PACK_MANIFEST.sha256` 在补充权重与结构文件后已按原格式
> **重生成**（601 → 621 条，格式与排序规则不变，未改动文件的哈希逐条一致，差异清单见 §5）。
> `src/SOP_5/PACK_MANIFEST.sha256` 与 `src/SOP_ML_run05/PACK_MANIFEST.sha256` 覆盖的 src 子包内容
> 未改动，保持原样。

---

## 4. 交付前检查清单

- [x] 三套 MD 候选清单齐全（20 / 19 / 22 行）并归一为 `results/results.csv` + `results.xlsx`
- [x] 一键入口 `./run.sh` 实测通过，产出与冻结清单逐行逐字段一致
- [x] 主运行入口 `python3 predict.py` 可独立生成最终结果文件
- [x] `requirements.txt` + `environment.yml` 列明版本；`README` 写明运行命令
- [x] Notebook 可执行且端到端跑通
- [x] 数据来源、版本、许可与哈希登记（`data/DATA_SOURCES.md` + `data/manifests/`）
- [x] 模型卡标注来源、版本、适用范围与已知局限（`models/MODEL_CARD.md`）
- [x] 无校名 / 标识 / 指导教师与学生姓名；无平台作业号与队列信息
- [x] 无硬编码绝对路径（`tests/check_isolation.py` 越界引用 0 处）
- [x] 三维结构文件（附件 5 §四）—— **已包含** `results/structures_3d/`（8 个文件，2026-10-03 补充，见 §5）
- [x] M2 导出权重 + 训练入口 + 训练/验证数据 + 训练日志随包（2026-10-03 补充，见 §5）

---

## 5. 2026-10-03 补充：随包权重 / 训练入口 / 三维结构文件

本次补充只**新增**文件与说明，未改动任何阈值、门、名单成员、冻结规则、随机种子与既有中间产物。

| 补充项 | 内容 | 来源与校验 |
|---|---|---|
| 模型权重 | `models/a2a_rf_v2_portable.npz`（4,595,821 B）、`models/a2a_rf_trained_portable.npz`（1,015,644 B）、`models/a2a_rf_portable.npz`（1,006,194 B） | 本项目 2026-09 训练运行的导出权重；入包时逐文件比对大小与 sha256（源副本与本包副本哈希一致）。三者实测结构：`n_trees=400`、`n_features=2048`、`n_classes=2` |
| 训练入口 | `train.py`（依赖 `src/features.py`、`src/metrics.py`、`src/__init__.py`） | 与 2026-09-18 训练运行同源。**本机未复跑**（打包环境无 rdkit/scikit-learn 运行时）；评审环境按 `requirements.txt` 安装后命令为 `python3 train.py` |
| 训练 / 验证数据 | `data/dataset_train.csv`（1,307 行）、`data/dataset_external.csv`（2,167 行） | 随包自有数据，字段与边界见 `data/DATA_SOURCES.md` §5 |
| 训练日志 | `logs/train.log` | 2026-09-18 训练与外部集评测记录（ROC-AUC 0.9872 / PR-AUC 0.9461 / BEDROC 0.8590 / EF@1% 11.60） |
| 模型卡 | `models/MODEL_CARD_A2AR-RF-v2.md`；`models/MODEL_CARD.md` 中指向未随包文件的引用已改为指向该卡 | 修正原「随包」指向不存在路径的问题 |
| 三维结构文件 | `results/structures_3d/`（8 个文件）：`4EIY_protein.pdb`、`4EIY_crystal_ligand_ZMA.pdb`、`ZMA.cif`、`4EIY_ZMA_redocked_poses_gnina.sdf`、`A2AR-05_MD05_docking_pose.sdf`、`A2AR-06_MD06_docking_pose.sdf`、`A2AR-05_MD05_A2A_complex_MD.pdb`、`A2AR-06_MD06_A2A_complex_MD.pdb` | 取自本项目 2026-09 SOP-05 结构产物；编号对应关系与坐标约定见该目录 `README.md` |
| 包内哈希清单 | `data/manifests/PACK_MANIFEST.sha256` 重生成 | 新增 20 条、内容变化 4 条、其余 597 条哈希逐条一致；清单本身不收录自己 |

> 说明：`models/MODEL_CARD_A2AR-RF-v2.md` 与 `logs/train.log` 记录的是本包随包的 `a2a_rf_v2_portable.npz`
> 的**训练运行**（4,975 阳性样本 + 平衡阴性）；它与 `models/MODEL_CARD.md` §2 表格中所述 **1,307 样本口径**
> 的 M2 描述并非同一次训练运行，两者均随包保留、口径以前者（含性能读数）为准做交叉参照。
