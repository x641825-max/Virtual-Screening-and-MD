# SOP_5 可视化目录（08_visuals）

> 数据快照：规则 1.7.0｜rules_sha256 `a983c4985ebf1fde…`｜计划 19 步｜读数生成 2026-09-23 19:14 UTC

本目录把 SOP_5 全链**可可视化部分**集中成一册图，供阅读、汇报与复核。所有数字都从链上产物现读（不重算口径、不改任何门/池/名单），由 `scripts/make_visuals.py` 一键重建。

## 1. 图目录

| 图 | 内容 | 说明 |
|---|---|---|
| [`figures/fig01_md_structure_cards.png`](figures/fig01_md_structure_cards.png) | MD 候选分子结构卡（20 = 18 候选 + 2 对照） | 每格一个分子的 RDKit 2D 结构与标签：SOP5-MD-序号 | CNS（CNS MPO ≥ 4）/ Peripheral / 对照；第二行为库 ID；△ 标记第二 hERG 模型分歧者。顺序 = CNN_VS 降序，对照固定末两位。 |
| [`figures/fig02_pipeline_flow.png`](figures/fig02_pipeline_flow.png) | 筛选流程图：gnina 对接表 → MD 候选分子 | 七段：固定输入 → 姿态核验与标准化 → 类药性硬门 → 姿态闸门与锚点 → 相似度与安全门 → 名单选序 → MD 批次。每段含规模、判据与「人工分析检查点」。 |
| [`figures/fig03_gate_cascade.png`](figures/fig03_gate_cascade.png) | 逐门级联 4,990 → 20 | 上：每道门的剩余量与剔除量（hERG 门收缩最强）；下：8 档臂池（QED 档 × 锚点档）在 hERG 门前后规模。 |
| [`figures/fig04_pose_identity.png`](figures/fig04_pose_identity.png) | 姿态与身份核验面板 | 姿态闸门 pass/warn/reject、锚点命中分布（含 5,098 实体库富集）、结构身份分级、maxTc 对数分布。 |
| [`figures/fig05_druglikeness.png`](figures/fig05_druglikeness.png) | 类药性面板 | QED 分布与 0.40/0.50/0.60 三档、CNS MPO 分布（含 cLogD 80% 区间悲观/乐观界）、Ro5/Veber/PAINS 触发、cLogP–TPSA 空间。 |
| [`figures/fig06_herg_safety.png`](figures/fig06_herg_safety.png) | hERG 安全面板 | 全池概率分布与 0.5 门、第二模型（ChEMBL 自建集成）判定对照、阈值保留曲线、0.5 在池/库/已上市药中的百分位。 |
| [`figures/fig07_robustness.png`](figures/fig07_robustness.png) | 稳健性面板 A/B/C/F | A 打分扰动 σ 扫描（含主判据线与通过与否）、B 阈值矩阵、C 单轴开关、F Butina 截断扫描。 |
| [`figures/fig08_list_and_batch.png`](figures/fig08_list_and_batch.png) | 交付名单 ↔ MD 批次画像 | 批次在 48 名单中的位置与中位线、QED/CNS MPO 覆盖、化学型配额、批次自检数字。 |
| [`figures/fig09_status_board.png`](figures/fig09_status_board.png) | 验收与确定性状态板 | 硬门 pass/pending/fail 逐条读数、G 臂 n_diff、关键质量数字、产物地图与转述注意。 |

网页版索引（自包含单文件、图片内嵌 base64，可离线分享）：[`gallery.html`](gallery.html)

## 2. 读图顺序

1. `fig02` 看全局：从 49,963 行排名表走到 20 个 MD 分子，每道关的判据与规模；
2. `fig03` 看每道门各裁掉多少 —— hERG 门（全池作用域）收缩最强；
3. `fig01` 看最终 20 个分子长什么样（结构卡，CNS/Peripheral/对照一目了然）；
4. 细节按步骤：`fig04` 姿态与身份 → `fig05` 类药性 → `fig06` 安全 → `fig07` 稳健性；
5. `fig08`/`fig09` 收口：名单与批次画像、验收与确定性状态。

## 3. 当前包状态（与图同一快照，逐条读数）

- 硬门验收：**6 通过 / 2 待判定（步序读数时差）/ 1 未通过**（`06_deliverables/acceptance_report.json`，共 9 条硬门）
  - 未通过：**σ=0.10 名单平均重叠 ≥ 0.8（主判据）** —— 实测 `0.7888`，期望 `0.8`
- 确定性 G 臂（链末文件 `05_robustness/reports/determinism_compare.json`）：`mode = compared`，`n_diff = 0`，文件数 187 / 187，`determinism_verified = True`
- 交付名单（主口径）：48 候选 + 2 对照；MD 批次 20 = 18 + 2

## 4. 本轮核对到的异常（需在链上处理，本目录不改链）

| # | 现象 | 证据 | 影响 |
|---|---|---|---|
| 1 | σ=0.10 名单平均重叠未达主判据 | 实测 `0.7888` < 0.80（`05_robustness/reports/robustness.json` → `pass_primary=false`） | 验收硬门 1 条未通过；对外须如实报告「分数扰动 σ=0.10 时名单不稳定」 |

读数时差说明：验收报告里那 2 条「待判定」不是缺陷，而是**报告与链末步的先后顺序**（`14_report` 早于 `12_determinism`，报告只登记上一轮读数）；链末文件本轮的实测读数为 `mode=compared` / `n_diff=0`，确定性判据可判为通过。唯一真正未达标的是 σ 主判据那一条，两者不要合并转述。

## 5. 图之外的既有图（各步骤自出，此处只做索引）

| 图 | 内容 |
|---|---|
| [`fig_pose_gate.png`](../01_pose_qc/figures/fig_pose_gate.png) | 步骤 1：姿态闸门与接触指纹（原文件存在：是）|
| [`fig_druglikeness_funnel.png`](../02_druglikeness/figures/fig_druglikeness_funnel.png) | 步骤 2：类药性漏斗（原文件存在：是）|
| [`fig_gate_stack.png`](../03_gates/figures/fig_gate_stack.png) | 步骤 3：门叠加（原文件存在：是）|
| [`fig_herg_cross_model.png`](../03_gates/herg/cross_model/figures/fig_herg_cross_model.png) | 步骤 3b：hERG 第二模型交叉验证（原文件存在：是）|
| [`fig_list_overlap.png`](../04_selection/figures/fig_list_overlap.png) | 步骤 4：名单重叠（原文件存在：是）|
| [`fig_gate_impact.png`](../04_selection/figures/fig_gate_impact.png) | 步骤 4：门对名单的影响（原文件存在：是）|
| [`fig_robA_sigma.png`](../05_robustness/figures/fig_robA_sigma.png) | 步骤 5A：σ 扰动（原文件存在：是）|
| [`fig_robB_matrix.png`](../05_robustness/figures/fig_robB_matrix.png) | 步骤 5B：阈值矩阵（原文件存在：是）|
| [`fig_robE_band.png`](../05_robustness/figures/fig_robE_band.png) | 步骤 5E：并列带宽（原文件存在：是）|
| [`fig_robF_clusters.png`](../05_robustness/figures/fig_robF_clusters.png) | 步骤 5F：Butina 聚类（原文件存在：是）|
| [`fig_md_batch.png`](../07_md_batch/figures/fig_md_batch.png) | 步骤 6：MD 批次分布与构成（原文件存在：是）|

## 6. 重建与边界

```bash
cd SOP_5
python3 08_visuals/scripts/make_visuals.py   # 重建本目录全部图 + README.md + gallery.html
python3 run_all.py                          # 复现整链（19 步 + 自检），约 5 分钟
python3 tests/verify_pack.py                # 包完整性自检
```

- 本目录**不参与** G 臂（`05_robustness/scripts/12_determinism.py` 的 `INCLUDE_DIRS` 不含 `08_visuals`），也不改任何门、池、名单与交付件；
- 图的内容随链上产物变化：**换快照（重冻/改规则）后必须重跑本脚本**，否则 `figures/` 会与包状态不一致；
- 图内全部为中文标注，字体由 `vendor/plots.py` 的 CJK 配置提供；PNG 200 dpi + SVG 双出，SVG 已固定 hashsalt 与元数据以保证复现。

## 7. 数据源对照

| 图 | 主要数据源 |
|---|---|
| fig01 | `07_md_batch/tables/md_batch20.csv`（SMILES、CNS MPO、批次序、`herg_cross_model_flag`）|
| fig02 | `03_gates/reports/gate_stack.json`、`01_pose_qc/reports/pose_gate_summary.json`、`01_pose_qc/reports/pose_identity_summary.json`、`02_druglikeness/reports/druglikeness_funnel.json`、`03_gates/reports/herg_status.json`、`03_gates/reports/maxTc_summary.json`、`04_selection/reports/selection_summary.json`、`07_md_batch/reports/md_batch.json` |
| fig03 | `03_gates/reports/gate_stack.json`（`gate_chain`、`arm_pools`）+ `04_selection/reports/selection_summary.json` + `07_md_batch/reports/md_batch.json` |
| fig04 | `01_pose_qc/reports/{pose_gate_summary,pose_identity_summary}.json`、`03_gates/reports/maxTc_summary.json` |
| fig05 | `02_druglikeness/tables/descriptors.csv`、`02_druglikeness/reports/druglikeness_funnel.json` |
| fig06 | `03_gates/admet/tables/admet_annotation_pool.csv`、`03_gates/reports/herg_status.json`、`03_gates/herg/cross_model/tables/herg_cross_model_predictions.csv` 与 `reports/…json`、`03_gates/admet/tables/herg_retention_curve.csv`、`03_gates/admet/reports/{admet_annotation,herg_threshold_analysis}.json` |
| fig07 | `05_robustness/reports/robustness.json`（armA/B/C/D/E/F + 全池 Butina 扫描）|
| fig08 | `04_selection/tables/list_main_QED0.50_anchor4.csv`、`07_md_batch/reports/md_batch.json` |
| fig09 | `06_deliverables/acceptance_report.json`、`05_robustness/reports/determinism_compare.json`、`REPRO_INFO.json`、`config/run_order.json`、`config/frozen_rule.json` |

## 8. 转述注意（写进对外材料前请先看这几条）

- **hERG**：单一外部模型（ADMET-AI）的分数，门 0.5 未校准（库内 66.4% 在线上）；第二模型（ChEMBL `CHEMBL240` 自建）与其 κ ≈ 0.02–0.07、门决策一致率 ≈ 33%–45%（近乎随机一致）。引用时必须写作「在 ADMET-AI 口径下过门」。
- **CNS MPO**：依赖上游 MolGpKa pKa 与 cpLogD7.4（第三方预测）；cLogD 带 80% 区间，`fig05` 已把悲观/乐观界一并给出（区间宽度中位 0.375 log）。
- **MD 批次**：机时分配层，不是「比其余候选更优」的结论；换人须改规则并重跑整链。
- **稳健性**：σ=0.10 名单平均重叠是主判据（阈值 0.80），现值见 `fig07` 与第 3 节；阈值矩阵、指纹更换（MACCS 0.5）说明名单对**相似度约束实现**同样敏感。
- **筛库形状**：全池 PAINS 命中 0 个，说明上游已做 PAINS 清洗；本包的 PAINS 判据在效度上属于「空转」而非「通过验证」。
