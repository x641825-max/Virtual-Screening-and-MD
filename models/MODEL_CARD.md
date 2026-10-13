# Model Card · 本项目使用的模型与打分器

按赛事《代码提交要求》附件 5 §三（第三方模型/平台须说明名称、版本、调用时间、输入与关键参数，
并标注来源与许可）编写。**默认复跑链路上不训练、也不重跑任何模型**：所有预测值在上游一次性生成后
冻结为输入表，包内复现的是「冻结预测值 → 门/标注/分层 → 候选清单」的确定性计算。

> **2026-10-03 补充（附件 5 §一(2)、§二(3)）**：本包现已随包提供 M2 的**导出权重**（`models/a2a_rf_*.npz`）、
> **训练入口** `train.py`（含其依赖模块 `src/features.py`、`src/metrics.py`）、**训练与外部验证数据**
> （`data/dataset_train.csv`、`data/dataset_external.csv`）、**训练日志** `logs/train.log`，以及 M2 的
> **完整模型卡** `models/MODEL_CARD_A2AR-RF-v2.md`。该训练入口不在默认复跑链路上（`run.sh` / `predict.py`
> 不调用它），仅供评审独立复核权重与训练过程。

---

## 1. 一览

| # | 模型 / 打分器 | 类型 | 版本 / 时间 | 输入 | 输出（在本项目中的用法） | 是否本包自训练 |
|---|---|---|---|---|---|---|
| M1 | **gnina CNN** | 3D 卷积打分（对接） | gnina 1.3.2；`semi_flexible`，`cnn_scoring=rescore`，`exhaustiveness=8`，`num_modes=9` | 受体 4EIY（A 链）+ 配体构象 | `CNN_VS = CNNscore × CNNaffinity`（**关键预测指标**，越大越好） | 否（上游对接，本包不重跑） |
| M2 | **A2A 便携初筛模型** | ECFP4（2048 bit, r=2, useChirality=False）+ 随机森林（主力）/ HistGradientBoosting / MLP | 训练日期 2026-09-14 | SMILES | 「A2A 配体样」概率，**仅用于排序富集** | **是**（上游自建；**已随包**：权重 `models/a2a_rf_*.npz`、训练入口 `train.py`、训练日志 `logs/train.log`） |
| M3 | **ADMET-AI hERG 头** | chemprop（D-MPNN）多任务 ADMET | 冻结预测表，2026-09 | SMILES | hERG 阻断概率 → **分类门 < 0.5** | 否（第三方开源模型） |
| M4 | **第二 hERG 模型（本项目自建集成）** | ECFP4 + RF / LR / kNN，取集成口径 | 2026-09-23 训练；ChEMBL 抓取 2026-09-23T06:47:45Z | SMILES | hERG 阻断概率 → **仅作交叉验证标注，不作判据** | 是（训练数据可由公开 ChEMBL API 复取） |
| M5 | **MolGpKa** | 图神经网络 pKa 预测 | provider version `4dc835219baad18c2e705274ac54d49597703966` | SMILES | 最强碱性中心 pKa → CNS MPO 的 `basic_pka` | 否（上游第三方） |
| M6 | **cpLogD v2.0** | logD 预测 | 随上游排名表冻结 | SMILES | cLogD7.4 点估计 + 80% 区间 → CNS MPO 的 `clogd74` | 否（上游第三方） |
| M7 | **生成模型（SCALE- 分子来源）** | 上游 ML 分子生成 | 本包不含生成过程与权重 | — | 生成 4,979 个候选分子，交由 M1 打分 | 否（上游，**本包不代为描述其架构**） |
| M8 | **选批求解器** | MILP（分支定界） | `scipy.optimize.milp`，HiGHS 后端；scipy 1.15.3 | 候选池 + 约束 + 目标权重 | 20 席名册与最优性证书（`mip_gap=0`） | 不适用（通用优化器） |

> M7 的说明只到「上游 ML 生成、编号 SCALE-」这一层——本工作区没有该生成模型的架构、训练数据与权重
> 记录，故**不作任何架构层面的断言**（未经检索核实）。

---

## 2. M2 · A2A 便携初筛模型（本项目自建，唯一有完整模型卡的自训练模型）

| 项 | 内容 |
|---|---|
| 用途 | 从大规模可购买库中**排序并富集** A2A 腺苷受体配体，供后续对接精筛 |
| 靶点 | 人 A2A 腺苷受体 ADORA2A（UniProt **P29274**；PDB **4EIY**） |
| 输出 | 属于「A2A 配体样」的概率（0–1），**仅用于排序**，不用于二元判定 |
| 训练集 | 本地 DUD-E 风格库 479 个 A2A 活性（69 个兼为双拮抗剂）+ 文献 349 个真实实验非活性（Ki > 1000 nM）+ 479 个 1:1 下采样 DUD-E 诱饵 = **1,307** |
| 内部测试 | 按 Bemis-Murcko 骨架切分（跨骨架泛化）135 条（53 活性） |
| 外部独立测试 | 文献 167 个训练中未见过的双拮抗剂新化学型 + 2,000 held-out 诱饵 = **2,167** |
| 数据来源 | 本地库 `data/dataset_train.csv`（**随包**）；文献数据 Wang M, et al. *PLOS Comput Biol* 2021, DOI [10.1371/journal.pcbi.1008821](https://doi.org/10.1371/journal.pcbi.1008821)（仓库 [github.com/Houshujing/Adenosine](https://github.com/Houshujing/Adenosine)）；清洗时剔除文献内部标签自相矛盾的 3 个结构并做 canonical 去重 |
| 已知局限 | 训练集含 DUD-E 诱饵（存在已知诱饵偏置）；只在排序口径下使用，不作富集验收判据 |
| 完整模型卡 | `models/MODEL_CARD_A2AR-RF-v2.md`（**随包**，A2AR-RF-v2 档案卡，含外部独立集全部性能读数） |
| 随包权重 | `models/a2a_rf_v2_portable.npz`（4.59 MB，主权重）、`models/a2a_rf_trained_portable.npz`、`models/a2a_rf_portable.npz`；三者均为纯 numpy 便携包，实测结构 `n_trees=400`、`n_features=2048`、`n_classes=2`（ECFP4-2048 + 400 棵树的随机森林） |
| 训练入口 | `python3 train.py`：读 `data/dataset_train.csv` 训练 400 棵树，导出 `models/a2a_rf_trained_portable.npz`，并在 `data/dataset_external.csv` 上打印 ROC-AUC / PR-AUC / BEDROC(α=20) / EF@1%（实现见 `src/features.py`、`src/metrics.py`） |
| 训练日志 | `logs/train.log`（2026-09-18 运行记录：4,975 阳性样本 + 平衡阴性，400 棵树，外部集 ROC-AUC 0.9872 / PR-AUC 0.9461 / BEDROC 0.8590 / EF@1% 11.60，导出 4.59 MB 便携包） |

---

## 3. M3 / M4 · hERG 双模型与它们的诚实边界（**重要**）

| 项 | ADMET-AI（M3，主口径门） | 第二模型（M4，集成口径） |
|---|---|---|
| 训练数据 | TDC ADMET / 公开 hERG 数据集（多源自 ChEMBL） | ChEMBL REST API，`CHEMBL240`（KCNH2 / hERG，人） |
| 标签规则 | 模型内置 | blocker：pChEMBL ≥ 6.0；non-blocker：≤ 5.0；(5.0, 6.0) 丢弃；InChIKey 取 pChEMBL 中位数 |
| 样本量 | — | 原始记录 13,521 → 可用标签 6,057（2,411 blocker / 3,646 non-blocker） |
| 特征 | chemprop D-MPNN | ECFP4 2048 bit（RDKit） |
| 性能 | 见上游模型卡 | RF 随机 5 折 ROC-AUC **0.9489** / PR-AUC 0.9323；骨架 5 折 ROC-AUC 0.8968 |
| 在本项目中的角色 | **冻结门**（hERG < 0.5） | **交叉验证标注，不作判据** |

**必须转述的两模型一致性**：门决策一致率仅 **0.33–0.45**，Cohen's κ **0.017–0.073**
（接近随机一致）；与阈值无关的排序一致性同样低（高风险前 10% Jaccard 0.09–0.13）。
因此**「hERG 过门」只能表述为「在 ADMET-AI 口径下过门」**；第二模型不过门的候选
**保留入池并逐分子标注「模型依赖」**，正交出口是膜片钳实验——MD 不回答安全性问题。

原始一致性读数见 `src/SOP_5/03_gates/herg/cross_model/reports/herg_cross_model.json`。

---

## 4. 数据泄漏防控与去重

- 结构身份统一用 **InChIKey** 核对：跨来源（目录 vs 生成）结构零重叠，重复仅出现在 2 个共用阴性对照。
- 名单层设 **novelty 门**：剔除与 ChEMBL 已知 A2A 配体精确重合者（目录链 −11 个、生成链 −11 个）。
- 相似度门 **maxTc(对已知 A2A 活性物) < 0.90**，降低「复现已知活性物」的假新颖。
- 骨架层：Bemis-Murcko 通用骨架去重（保留 CNN_VS 高者）+ 名单内骨架 ≤ 1 + 成对 ECFP4 Tc < 0.65。
- M2 的内部评测按**骨架切分**（而非随机切分），外部测试用训练未见过的化学型，避免同类结构泄漏高估。

---

## 5. 复现方式（环境与入口）

```bash
python3 -m pip install -r requirements.txt   # 或 conda env create -f environment.yml
./run.sh                                     # 复跑「选批 → 合并 → 结果文件」并自检
python3 predict.py                           # 生成 results/results.csv 与 results.xlsx
python3 train.py                             # 可选：按 data/dataset_train.csv 复训 M2 并导出便携权重
```

M1/M3/M4/M5/M6/M7 的**预测层**不在包内重算（无相应运行时与权重）；
如需重生成，接口与再生成路径登记在各包文档中（如 `src/SOP_5/03_gates/herg/HERG_INTERFACE.md`、
`src/SOP_5/03_gates/scripts/10b0_admet_regenerate.py`，后者为手工步、不入链）。
