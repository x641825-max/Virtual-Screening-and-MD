# SOP_5 · A2A 拮抗剂 MD 候选名单工程包

自包含、结构化、可复现。**所有输入与产物都在本包内**；脚本不在包根（`REPRO_INFO.json`
所在目录）之外解析任何路径（`tests/check_isolation.py` 静态扫描保证）。

## 1. 一句话说明

从 5,000 个平台回传的 gnina 最优姿态出发，逐层施加「结构身份 → 类药性 → 结构/活性相似度
→ 毒性（hERG，外部接口）→ 多样性约束」的确定性门禁与约束贪心选序，产出 **50 个 MD 起始
分子**（48 候选 + 2 阴性对照）及其逐分子证据链、6 组名单矩阵、稳健性检验与验收对账。

## 2. 快速开始

```bash
cd SOP_5
python3 run_all.py --list          # 查看步序（19 步；其中 08b/10b/10c/15 为分析步，不改门/名单）
#   步序依赖：08b 必须晚于 10_select_lists（它要读主名单算名单级影响）；
#   单跑某一步可能因缺少上游产物而显式失败（不静默降级），用 --from 或全链运行
python3 run_all.py                 # 全链复现（末尾自动自检 + 确定性比对）
python3 run_all.py --from 10_select_lists    # 从某步继续
python3 run_all.py --freeze        # 重新固化确定性基线（改过规则/脚本后才需要）
```

产物目录：`06_deliverables/`（详见该目录 `INDEX.md`）；中期产物在 `01_pose_qc/`…`05_robustness/`；
运行日志在 `logs/`。

**重跑与清空的边界**（起链前 `run_all.py` 会按 `config/run_order.json` 的 `outputs` 预建所有产物父目录，
并把补齐清单打进 `logs/run_all.log`，故「整目录被删」与「只删文件」的重跑行为一致）：

| 可以删（会按需重建） | 不可删（属输入或代码层） |
|---|---|
| 各步骤的产物子目录与文件：`tables/`、`reports/`、`figures/`、`outputs/`、`00_inputs/manifests/` 的清单、`06_deliverables/` 的交付件 | `00_inputs/{raw,extracted}/`、`config/`、`vendor/`、`tests/`、`docs/`、各步 `scripts/`；以及**冻结输入层** `03_gates/admet/raw/`、`03_gates/herg/cross_model/raw/`、`03_gates/herg/herg_predictions.csv`（外部落门表）、`03_gates/herg/hpc_herg/`（惰性审计证据） |

冻结输入层被删时步骤会**显式失败并给出恢复指引**（`08b` 与 `10_select_lists` 均有前置检查），
不会静默重抓或降级。

## 3. 输入（仅 5 个，全部在包内固化哈希）

| 文件 | 角色 | md5 |
|---|---|---|
| `00_inputs/raw/inputs/top5000_best_pose.sdf` | 5,000 个 gnina 最优姿态 | `0ed679e1d61d7f562b27bfcdd1c7e956` |
| `00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv` | 全库 CNN_VS 排名 + 上游描述符 + MolGpKa pKa + cpLogD7.4 | `444f07fbe92cd15572ea459e6e955c66` |
| `00_inputs/raw/archives/run01_redock_enrichment.zip` | 受体 4EIY、5,098 带标签富集库（479 活性）、重对接与富集基线 | `82a64a41b06f49dc4327b6bd7fcb135e` |
| `00_inputs/raw/archives/screening_kit.zip` | ChEMBL 独立已知 A2A 配体集、模型卡、受体副本 | `ff5902239a51536d6bc34e42599ceac3` |
| `00_inputs/raw/archives/redocking.zip` | 3EML / SOP-02 重对接与诱饵基准（协议溯源） | `cac648447c54f688562ddafed90dcf65` |

归档按**白名单**解包到 `00_inputs/extracted/`（见 `00_inputs/scripts/00_unpack_inputs.py`
与 `00_inputs/manifests/EXTRACTED_MANIFEST.md5`）；大体积成员（5,098 集姿态 130 MB）
由 `config/run_order.json → options.extract_large` 控制，关闭时相关步骤**显式报错**而非静默降级。

另有 **3 个外部模型产出/抓取的数据表**在本包内固化并登记哈希（清单见 `REPRO_INFO.json → extra_inputs`）：
`03_gates/herg/herg_predictions.csv`（hERG 落门表，四臂并集 3,236）、
`03_gates/herg/cross_model/raw/chembl_herg_kcnh2_activities.csv`（第二模型训练/标注集）、
`03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv`（ADMET-AI f3 库全量表，49,963 × 124，
sha256 登记于 `config/frozen_rule.json → step3_gates.admet_layer`，81.1 MB）。

## 4. 口径（唯一来源：`config/frozen_rule.json`）

| 环节 | 口径 |
|---|---|
| 姿态核验 | SDF 属性与上游表逐值一致；`CNNscore×CNNaffinity == CNN_VS`（±5e-4）；坐标往返偏差 0 |
| 身份分级 | identical / stereo_only / protonation_diff / graph_mismatch（姿态 3D 指派 vs 上游 SMILES） |
| 去重 | InChIKey（上游 SMILES），保留 CNN_VS 最大者 |
| 接触/锚点 | 重原子最近距离 < 4.0 Å；口袋 16 残基；锚点 N253/E169/F168/W246 → `anchor_hits_4A` |
| 姿态闸门 | reject＝出盒或锚点 ≤1；warn＝锚点 =2；pass＝锚点 ≥3（warn 保留在池、不进主名单） |
| 类药性硬门 | PAINS（RDKit FilterCatalog 重算）=0 **且** Ro5 违规 =0 **且** Veber 违规 =0 |
| QED | 主口径 ≥0.50；敏感性 ≥0.60 / ≥0.40 |
| 标注列 | CNS MPO（Wager 2010；输入 `basic_pka`＝MolGpKa 点值、`clogd74`＝cpLogD v2.0 点估计 + 80% 区间）、`cns_mpo_lo80`/`hi80`/`pi80_span`（区间界定）、TPSA∈[90,140]、Brenk |
| 结构/相似度门 | 姿态 pass ∧（anchor=4 主 / ≥3 敏感）∧ maxTc(479 已知活性) <0.90（敏感 0.70） |
| hERG | <0.5，**外部接口**：预测表 `03_gates/herg/herg_predictions.csv`（ADMET-AI，已登记为冻结输入并**已落门**；缺表则 pending）。补算清单 `03_gates/herg/SOP-05_herg_input_union.csv` |
| 已知物（novelty） | 选序层剔除 InChIKey 与 ChEMBL 独立 A2A 配体集精确重合的分子（验收要求名单内已知物 = 0） |
| 选序 | CNN_VS 降序贪心 + 成对 ECFP4 Tc<0.65 + 一骨架一分子（敏感 ≤3/≤5）+ 化学型 ≤ round(8N/20) + AM1-BCC 风险档 ≤ round(3N/20) |
| ADMET 注释层（分析步 10b，**非判据**） | 读包内冻结的 ADMET-AI f3 全量表 → 池级标注列（毒性/ADME/CYP/描述符/已上市药百分位，43 列）+ 池内重标定切点 + RDKit 独立复算交叉核验 + MD 前置复核清单；**不改门、不改池、不改名单** |
| hERG 阈值锚定（分析步 10c，**非判据**） | 三锚定（① 库内相对仅排序 ② 已上市药参照定软线 ③ 取舍驱动定硬线候选）+ 0.01 步长保留曲线 + 名单档位标注 + 门作用域缺口披露；主口径仍 `hERG < 0.5` |
| 对照 | 2 个阴性对照（规则冻结，包内从 run01 归档重算 NEG1/NEG2） |

## 5. 交付清单（`06_deliverables/`）

| 文件 | 内容 |
|---|---|
| `SOP-05_MD_shortlist_main.csv` | 主名单 50（48 候选 + 2 对照），含逐分子门与注释 |
| `SOP-05_MD_shortlist_sens_*.csv`、`..._scaffoldcap{3,5}.csv` | 5 个敏感性名单 + 2 个骨架封顶档 |
| `SOP-05_decision_cards_main.csv` | 决策卡：分数、门、注释、最近已知活性、pKa/cLogD |
| `SOP-05_negative_controls.csv` | 阴性对照卡（规则 + 实测） |
| `SOP-05_poses_main_list.sdf` | 主名单候选起始姿态（几何与上传一致） |
| `acceptance_report.json` | 验收对账（逐条判据 + 实测 + 上下文） |
| `SOP-05_report.md` | 最终报告 |
| `SOP-05_downstream_MD_plan.md` | 下游 MD：体系、两级判据、机时预算 |
| `SOP-05_herg_cross_model_consistency.csv`、`..._rank_consistency.csv`、`..._list_flips.csv` | hERG **第二模型**（ChEMBL 自建 RF/LR/kNN）交叉验证：与 ADMET-AI 的阈值一致性、与阈值无关的排序一致性、交付名单中偏离的分子 |
| `SOP-05_admet_annotation_main.csv` | 交付 48 候选的 ADMET 注释全表（冻结表原精度，阴性对照标 `not_covered`） |
| `SOP-05_MD_batch20.csv` | MD 一级批次 20（18 候选 + 2 阴性对照），含入选理由、第二模型 hERG 与 CNS MPO |
| `SOP-05_MD_batch_excluded.csv` | 未入批分子及理由（hERG 第二模型不过门 / 通用骨架重复占位） |
| `SOP-05_MD_batch20.md`、`SOP-05_md_batch.json` | 批次选择说明（人读）与机器可读留痕 |
| `SOP-05_admet_md_precheck.csv` | MD 前置复核清单（带 ≥1 项 ADMET 警告者，13 个） |
| `SOP-05_admet_descriptor_crosscheck.csv`、`SOP-05_admet_descriptor_crosscheck.json` | 描述符/Ro5/PAINS/Brenk 独立复算交叉核验（数值 max\|Δ\|=0；BRENK 口径差异已登记） |
| `SOP-05_admet_annotation.json`、`SOP-05_admet_layer.md` | 注释层证据与说明（覆盖、重标定切点、触发率、分诊层） |
| `SOP-05_herg_retention_curve.csv`、`SOP-05_herg_threshold_anchors.csv`、`SOP-05_herg_threshold_analysis.json`、`SOP-05_HERG_THRESHOLD_SOP.md` | hERG 阈值三锚定与保留曲线（含输入 sha256 与复现核对） |
| `MANIFEST.sha256`、`INDEX.md` | 交付哈希清单与索引 |

## 6. 复现与核验

```bash
python3 tests/verify_pack.py       # 输入哈希、规则一致性、产物规模、名单硬判据、CNS MPO 回归
python3 tests/check_isolation.py   # 确认无包外路径引用
python3 tests/make_index.py        # 重新生成 INDEX_PACK.md 与 PACK_MANIFEST.sha256
python3 05_robustness/scripts/12_determinism.py   # G 臂：全量产物 sha256 比对
```

确定性：主链不使用随机数；稳健性臂使用 `PROJECT_SEED=20260905` 的显式随机源；
`run_all.py` 强制 `PYTHONHASHSEED=0`。

**G 臂（确定性重跑）口径**：比对范围 = 全部解析型产物（各步骤 `tables/`、`outputs/`、`figures/`、
`reports/` 与 `00_inputs/manifests/`）。6 个**自指文件**不在范围内——2 个比对报告本身，以及
交付层 4 个「状态载体」（最终报告、交付索引、验收对账、MANIFEST），因为它们把 G 臂结论写进
自身内容，纳入范围会造成 True/False 振荡、永无不动点。

引导流程：`python3 run_all.py --freeze`（写基线）→ `python3 run_all.py`（应报 n_diff=0）。
**从零重建的实证**（2026-09-23，沙箱副本已删）：删掉 **155 个产物文件并把 14 个产物目录整体 `rmtree`**
（仅保留脚本层、输入层与 G 臂基线）后执行 `python3 run_all.py` → ALL DONE（19 步 + 自检，约 308 s），
链末 G 臂 `mode=compared`、`n_diff=0`（187/187 个解析型产物与冻结基线逐字节一致），
`verify_pack` 全过（该轮 179 项；加入 MD 最优性证书组 [21] 后为 193 项 / 21 组）；与主包逐文件比对 213 个文件一致，差异只在 7 个自指状态载体
（见 `tests/verify_pack.py` 第 [20] 组与 `../SOP_5_review/SOP-05_运行逻辑审查_2026-09-23.md` §5.3）。
**v1.7.0（MD 批次改为可证明最优）的收敛证据**（2026-09-23，规则 v1.7.0 / `rules_sha256 = a983c4985ebf1fde…`）：`run_all.py --refreeze` 写基线 → 连续四轮普通跑，每轮均 `ALL DONE（19 步 + 自检，302–307 s）`、G 臂 `mode=compared` / `n_diff=0`（187/187 逐字节一致）、`determinism_verified=true` → **载体不动点 Z3 ≡ Z4（7/7 逐字节一致）**，证据 `logs/history/carrier_round{Z0…Z4}.md5`；随后链外重建图册（`08_visuals` 不进 G 臂，但其哈希进根清单）并重生成根清单，再跑一轮 → **W0 ≡ W（7/7）**，证据 `logs/history/carrier_round{W0,W}.md5`。终态：`verify_pack` 193 项 / 21 组全过、`check_isolation` 越界 0 处、硬门 9 条 = 6 通过 / 2 待 G 臂时点 / 1 未通过（σ=0.10 主判据，既有诚实项）。

**手工图册不在链上**：`python3 08_visuals/scripts/make_visuals.py` 重建 `08_visuals/`（9 张图 + README + 单文件
`gallery.html`，约 15 s）；该目录不在 G 臂范围，换快照后需手工重跑一次以保持与包状态一致。
**状态载体的稳定规则**：交付层模板（13/14）或脚本改动后，载体需「冻结后再跑 **2–3 轮**」
才收敛——因为 `SOP-05_determinism_compare.json` 在 13 之后才写、而 `MANIFEST.sha256` 覆盖它，
`MANIFEST` 因此恒慢一轮；实测 A≠B≠C、C==D（见 `logs/history/cc3_{A..D}.md5`）。判据取
「连续两轮载体逐字节一致 + G 臂 n_diff=0」。最近一次固化证据（19 步版见 `logs/history/`）：
`logs/history/carrier_round{B,C,D}.md5` —— 7 个载体/自检文件中，C 与 D 两轮
**全部逐字节一致**（B→C 仍有 `MANIFEST.sha256` 与 `PACK_MANIFEST.sha256` 变化，
即上述滞后一轮：`MANIFEST` 覆盖在其之后写出的 `determinism_compare.json`，
`PACK_MANIFEST` 又覆盖 `MANIFEST`）。

**包级清单（`PACK_MANIFEST.sha256`）**：由 `tests/make_index.py` 生成，**排除自身**
（自指会让清单每轮都变，永远无法自校验——`verify_pack` 的 `[12]` 块对此加了断言）。
`logs/run_all.log` 每轮滚动：上一轮整体移入 `logs/history/run_all_<UTC>.log`，文件头写入
运行时间、步骤清单与规则哈希，因此该文件始终只记录**当前这一轮**。
另注：报告只嵌入稳定字段（文件数 / 差异数 / 结论），刻意不嵌入 `mode`
（`baseline_written` 与 `compared` 的取值差会让载体在冻结后多漂移一轮）。
若中途改过规则或脚本，重复一次该两轮流程即可（脚本哈希登记在 `00_inputs/manifests/FREEZE.json`，
规则哈希登记在 `config/frozen_rule.json` 的 `rules_sha256`；改规则须先把它置为 `"PENDING"`）。

## 7. 已知局限（必须在引用本名单时一并转述）

1. **hERG 门已闭合（残余不确定性需转述）**：本地无 ADMET-AI/chemprop 运行时与权重，
   hERG 门做成外部接口；本轮外部 ADMET-AI 预测表（HPC 应用 `admet-ai`，3,236 个四臂并集分子）
   已落入接口位置并按 <0.5 硬落门（并集内覆盖 100.00%、通过 1,005 / 不通过 2,231，名单与各臂池
   均已剔除不通过者）。残余不确定性：① 并集外 1,754 个池内分子无预测值（按规则保守不过门，
   但它们本就不在任何臂/名单内）。**该保守口径的代价现已量化**（分析步 10c，非判据）：包内冻结的
   ADMET-AI 全量 f3 表（49,963）显示这 1,754 个中 **480 个（27.37%）的 ADMET-AI hERG 实际 < 0.5**；
   落门表与冻结表在同源交集 3,236 个内逐值一致（max\|Δ\|=1e-6）。若要改作用域，属口径变更
   （需重跑选名单与验收基线），本包不改门、不改名单；② 外部模型版本自身的适用域限制；③ 应用接口不回传权重哈希，
   故只登记输出表哈希与输入 `molecules.csv` md5。缺表时名单不删分子、全部标注 `herg_gate=pending`。
   ④ **hERG 判定强依赖于所选模型**：用第二个独立训练的 hERG 模型（ChEMBL `CHEMBL240` 自建，
   13,521 条活性 → 10,928 分子 → 6,057 标注，仅 ECFP4 特征）对同一 3,236 个并集分子交叉验证：
   门决策一致率仅 33.3%–45.0%（Cohen κ 0.017–0.073），Spearman ρ 0.096–0.182；
   **与阈值无关**的排序一致性同样低（高风险前 10% 集合 Jaccard 0.09–0.13；匹配过门率后
   0.21–0.27，仅覆盖 ADMET-AI 过门分子的 34.5%–57.1%），说明差异不只是阈值/标定问题。
   若改用「两模型一致过门」，hERG 前主池由 609 变为 551、交付 48 名单中 3 个会掉出
   （详见 §7-1④ 与 `03_gates/herg/cross_model/HERG_CROSS_MODEL.md`）。
   因此**不要用任一单模型 0.5 阈值永久剔除分子**，建议把「两模型一致判高风险」作为 MD 前置
   复核的优先级依据；本包主口径仍保持 ADMET-AI 单模型门不变（改动需重新选名单并重设验收口径）。
   ⑤ 独立性的诚实边界：ADMET-AI 的 hERG 头同样基于公开 hERG 数据（含 TDC hERG 集，多源自
   ChEMBL）训练，故上述量化的是**算法/训练管道依赖度**，**不是数据独立的外部验证**。
2. **富集指标不用于验收**：富集评价参考集为 DUD-E 诱饵（存在诱饵偏置），
   本包只把锚点命中率作为结构门证据。
3. **CNS MPO 依赖上游第三方模型，仅作标注列**：`basic_pka`＝MolGpKa 最强碱性中心 pKa
   （只有点值、无预测区间；无碱性中心者上游记 0，全表 192 个、本包池内 4 个、交付名单 0 个，
   按 pKa≤8 取满分）；`clogd74`＝cpLogD v2.0/ChEMBL33 点估计 + **80% 预测区间**。
   按该区间界定重算 cLogD 轴：池内 2,778/4,990 个 CNS MPO 变动 ≥0.1，交付 48 名单中 22 个变动
   ≥0.1（最大 0.3985，`ZINC000092920989`），点估计 ≥4 的 40 个里有 **8 个在悲观界回落至 4 以下**（最低 3.453）。
   决策卡给出 `cns_mpo_lo80`（悲观，cLogD 取区间上界）/ `cns_mpo_hi80`（乐观，取区间下界）
   与 `cns_mpo_pi80_span`；引用「CNS MPO≥4」时须一并转述。pKa 轴无区间，不确定性无法量化。
4. **不重跑对接**：姿态由上游提供；本包核验其一致性与几何，不重新对接/打分。
5. **自由度声明**：`化学型`（氨基烟腈母核 × 碱性胺）是工程化的化学系列近似，
   与 `Bemis-Murcko 骨架` 约束互为补充；两者定义均冻结在规则集中。
6. **打分扰动的头部敏感**：σ=0.10 时名单平均重叠 0.8046（主判据 ≥0.80，**贴线通过**；
   min 0.6875、P10 0.750），头部 20 平均重叠 0.86（次要指标，不设门）。σ=0.15 时平均重叠
   降至 0.7202 —— 头部名次对分数噪声更敏感，引用头部顺序时需一并转述。
7. **单分子边界差**：QED≥0.40 档实测 3,987，与既有会话数字 3,988 差 1 个分子；
   根因为包内用 RDKit 2025.09.5 重算 QED 与上游 `qed` 列在阈值边界的极小浮点差。
   主口径（QED≥0.50）与硬门（4,577）不受影响。
8. **Butina 簇数与旧口径不可比**：本包 ECFP4 距离阈值 0.35（Tanimoto≥0.65）在全池得
   2,322 簇 / 最大 144 / 单例 1,639；既有会话的 642/20/446 在显式阈值下未能复现，
   已给出阈值扫描（0.20/0.30/0.35/0.50/0.65）供对账，不据旧数字做验收。
9. **MD 批次是机时分配层，「最优」有明确限定口径**：步骤 6（`15_md_batch_select`）从 50 名单里挑 20 个进一级 MD，只做取舍与留痕，**不改名单与门**。v1.7.0 起选批为**可证明最优**：在配额可行域内用 MILP（`scipy.optimize.milp` / HiGHS）按字典序求解「① 最大化分层覆盖下限 t（CNN_VS / cLogP / CNS MPO 三轴各 4 分位，资格池自分布）→ ② 最大化 CNN_VS 跨度 → ③ 最大化 Σ CNN_VS（质量兜底）→ ④ 规范并列打破」，并出具 `mip_gap=0` 的最优性证书（`07_md_batch/reports/md_batch.json → optimality`，交付副本 `06_deliverables/SOP-05_md_batch.json`）；v1.6.0 的贪心填充降为对照实现（`heuristic_reference`）。实测最优解 t\*=4（12 格铺满）、跨度 0.7337（池内全量程）、下半区 8 席，代价是 Σ CNN_VS 比质量口径低 1.16%（名册与贪心对照 6 换 6）——**「最优」限定在该目标函数口径内，不是「这批分子比其余候选更好」**；低分位（CNN_VS Q1）分子入批后失败概率更高，这正是要**测量**的量（分数-稳定性关系需要量程与分位两端的样本），不得事后解读为「选坏了」。两个第二 hERG 模型分歧分子**保留入池并逐分子标注「模型依赖」**（机时代价 +0.2 卡日量级），其安全性需另用正交方法解决（MD 不回答安全性）。CNS MPO 轴依赖上游 MolGpKa/cpLogD7.4 预测（见第 3 条）。
10. **确定性判据的时点**：验收表中的 G 臂结果取自**最近一次**确定性比对（G 臂位于整链末步），
   故首轮为 pending、第二轮起为实测。报告内的「数据截止」为规则/报告冻结日 2026-09-23。

11. **ADMET 注释层是标注、不是判据；预测层不在包内重算**：包内无 ADMET-AI/chemprop 运行时与
   权重，复现的是「冻结预测值 → 43 个标注列 / 池内重标定切点 / 阈值锚定 / 保留曲线 / 描述符
   交叉核验」的确定性计算；再生成路径见 `03_gates/scripts/10b0_admet_regenerate.py`（手工执行、
   不入链）。注释层不删分子——交付 48 候选中 2 个被注释层标为「建议淘汰」仍保留在名单内
   （现状由既有门决定）；预测值非实验数据，采购/合成/MD 立项前需体外确认。
   另需注意：本包 hERG 列与外部 run 的 hERG **同源同值**，不能当作第二模型/独立验证。

## 8. 参考来源

见 `06_deliverables/SOP-05_report.md` §12；方法学依据列表亦见 `docs/METHODS.md`。
