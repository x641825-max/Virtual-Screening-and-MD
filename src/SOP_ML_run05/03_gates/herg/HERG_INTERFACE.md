# hERG 门外部接口（SOP_5）

## 为什么需要外部接口

本工程包内**没有** hERG 预测实现：本地环境未安装 `admet_ai` / `chemprop` / `torch`
等运行时，且 ADMET-AI 权重不在包内。为保持"包内自包含、可复现"，本步只生成
**待补算清单**并预留落门位置，不伪造预测值。

## 需要外部产出的文件

路径（包内相对）：`03_gates/herg/herg_predictions.csv`

| 列 | 类型 | 说明 |
|---|---|---|
| `cid` | str | 与清单 `SOP-05_herg_input_union.csv` 的 `cid` 一致 |
| `hERG` | float | hERG 阻断概率，0–1，越大越危险（ADMET-AI 口径） |
| `source`（可选） | str | 模型名与版本，例如 `ADMET-AI v1.x / chemprop-2.x` |
| `weight_sha256`（可选） | str | 权重文件哈希，用于固化 |

参考实现（外部执行，命令示例）：

```bash
# ADMET-AI（Chemprop 后端）；需联网下载权重，建议在包外独立环境执行
pip install admet-ai
admet_predict --input-file SOP-05_herg_input_union.csv --output-file herg_predictions.csv
# 该命令输出的 hERG 列名为 "hERG"（阻断概率）
```

## 落门规则

`hERG < 0.5` 通过（阈值见 `config/frozen_rule.json → step3_gates.herg.threshold`）。
缺值**不得**当作「判为高风险」：状态记为 `not_evaluated`（未评估），并在 `herg_status` 与
报告中同时报出**作用域覆盖率**（`scope_coverage = 并集内分子数 / 池规模`）。
把缺值写成「保守不过门」会把「没测」冒名为「不合格」，属口径混淆。

## 外部产出溯源（登记于 `config/frozen_rule.json → step3_gates.herg`）

- 落门表：`03_gates/herg/herg_predictions.csv`（md5 `b22514dfbfcdbb98aa12fec9073a4e3e`、sha256 `81e30ec27ce09550f0707cc05a5b668ee0b32d02507489297ce0cfc526204766`）
- 来源：ADMET-AI (HPC app 'admet-ai') / run 2026-09-24
- 外部队列：HPC 目录应用 admet-ai；job_id eade0be46cb14cbb8b33b0ca46f43bf9；输入 03_gates/herg/hpc_herg/molecules.csv（4,979 行 cid,smiles）；原始输出 03_gates/herg/hpc_herg/predictions.csv（4,979 行 × 106 列）；落门表 = external_table（hERG 概率保留 6 位小数；本批覆盖全 4,979，与 scope=pool 对齐）

## 当前状态

- 清单规模：**3294** 个分子（四臂并集；规则期望 3236）
- 预测表是否存在：**True**
- 状态：**applied**
- 作用域覆盖率：**100.00%**（3294/4979；池内未评估 0 个）
- 为什么作用域小于池：hERG 外部补算的范围是四臂并集（代价控制），而非整个池。
  池内并集外的分子**没有预测值**，其门状态是 `not_evaluated`，不是「高风险」。
  当前并集外的分子均不在任何交付名单中，但下游若扩大候选池必须同步扩大补算范围，
  否则会把「未评估」当作「不合格」使用。

## 模型依赖度（第二模型交叉验证 · 分析步 08b）

本接口的 `hERG` 列来自**单一模型**（ADMET-AI）。为量化该判定的模型依赖度，链上分析步
`08b_herg_cross_model` 用第二个**独立训练**的 hERG 模型（ChEMBL `CHEMBL240` 活性数据自建，
特征仅 ECFP4）复算同一批四臂并集分子，给出阈值一致性、与阈值无关的排序一致性，
以及对池/名单的实测影响。产物见 `03_gates/herg/cross_model/`（`HERG_CROSS_MODEL.md`、
`reports/herg_cross_model.json`、`tables/`、`figures/`）。

要点（数字见上述产物与最终报告 §4b）：两模型在 0.5 阈值上的门决策一致率仅三分之一到一半，
且**与阈值无关**的排序交集同样低，说明差异不只是标定问题；因此
「hERG 过门」应转述为「在 ADMET-AI 下过门」，不宜当作可迁移的硬判据。
独立性的诚实边界：ADMET-AI 的 hERG 头同样基于公开 hERG 数据（含 TDC 集，多源自 ChEMBL）训练，
故该步量化的是**算法/训练管道依赖度**，不是数据独立的外部验证。

## 降级声明（表缺失时）

`hERG` 门记为 `pending`：名单**不因 hERG 被删减**，所有名单与决策卡带
`herg_gate=pending` 标注。此时：

1. 名单在"结构与活性相似度"两轴上仍是严格达标的（anchor/姿态/maxTc 均实测）；
2. hERG 风险需在 MD 前置的人工复核环节补齐，不得视为已通过；
3. 报告中 `hERG 区间` 一栏显示 `pending（未补算）` 而非区间值。
