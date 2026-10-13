# 方法与来源（SOP_5）

## 1. 口径来源

所有阈值与算法口径集中在 `config/frozen_rule.json`（`rules_sha256` 为其规范哈希，
脚本不得内嵌阈值）。本文件只解释**为什么这样定**，不重复数值。

## 2. 步骤与依据

| 步骤 | 做什么 | 依据/理由 |
|---|---|---|
| 1 姿态核验 | SDF 属性与上游排名表逐值比对、`CNNscore×CNNaffinity == CNN_VS`、坐标往返偏差、身份分级（identical/stereo_only/protonation_diff）、InChIKey 去重 | 姿态是后续全部结构证据的物理基础；身份分级显式区分「良性差异」（3D 指派手性、质子化/互变）与真实图不一致，避免把良性差异当成错误 |
| 1c 接触指纹与锚点 | 16 残基重原子最近距离 < 4.0 Å；锚点 N253/E169/F168/W246；盒子几何 | 口袋定义与 run01 `contact_fp.py` 完全一致（逐距离 0 偏差复现），锚点是 A2A 拮抗剂已知药效团接触点 |
| 2 类药性 | RDKit 重算 PAINS/Brenk/Ro5/Veber/QED/CNS MPO/化学型/骨架 | 上游 `pains` 列全空 → 必须重算；不预先圆整描述符（否则边界分子会漂移）；CNS MPO 按原文 Table 1 的分段线性函数实现，并以原文 Table 3 三例回归 |
| 3 门禁 | 姿态 pass ∧ anchor=4（主）；maxTc(479 活性) <0.90；hERG <0.5（**包内冻结 ADMET 全量表，作用域=全池；外部表作逐值对账件**） | anchor 命中率在带标签集上有区分度（活性 vs 诱饵）；maxTc 控制与已知活性的相似度上限；hERG 是 MD 前的安全前置 |
| 4 选序 | CNN_VS 降序 + 成对 Tc<0.65 + 骨架 ≤1 + 化学型 ≤ round(8N/20) + 风险档 ≤ round(3N/20) + 2 对照 | 同时约束「局部相似」「骨架多样性」「化学系列多样性」「参数化成本」，避免名单被单一母核/骨架挤占 |
| 3b′ hERG 第二模型交叉验证（08b，分析步；**步序上位于 10_select_lists 之后**，因其需读主名单算名单级影响） | ChEMBL `CHEMBL240` 活性（IC50/Ki，`standard_relation='='`，pChEMBL 非空）→ InChIKey 去重取中位数 → 标注（pChEMBL≥6 阻断 / ≤5 非阻断 / 之间丢弃）→ 仅 ECFP4 特征训 RF + 逻辑回归 + Tanimoto-kNN → 对四臂并集复算，报阈值一致性与**与阈值无关的排序一致性**、对池/名单的影响 | 外部单模型硬门必须量化其模型依赖度；用「匹配过门率」的对齐比较把阈值/标定因素剥离，避免把标定差异误读为模型分歧；该步只产证据，不改任何门与名单 |
| 4′ ADMET 注释层（10b，分析步） | 读包内冻结的 ADMET-AI f3 全量表（49,963 × 124）→ 池级 43 个标注列（毒性/ADME/CYP/描述符/已上市药百分位）+ 池内重标定切点（避免照搬外部 54,932 库的切点）+ RDKit 独立复算交叉核验（MW/cLogP/TPSA/HBD/HBA/QED/Ro5/PAINS/Brenk）+ MD 前置复核清单 | 把「分子画像」变成可复核证据而非印象：同一化学空间的预测值随表固化，切点按本包池的概率分位重标定；交叉核验让注释表与包内口径互相确认（独立实现给出同一判定） |
| 4″ hERG 阈值锚定（10c，分析步） | 三锚定（① 库内相对仅排序 ② 已上市药参照定软线 ③ 取舍驱动定硬线候选）+ 0.01 步长保留曲线（库/库级头部/池/主池/四臂并集/交付名单 9 个头部）+ 名单档位标注 + 门作用域缺口量化 | 阈值不能「凭感觉取整」：把阈值拆成排序层、解释层、决策层三问，用保留曲线显式给出「改阈值要牺牲多少对接头部」的代价，并在冲突时按合理性窗口消解；该步只产证据，不改门与名单 |
| 5 稳健性 | A σ 扰动、B 阈值矩阵、C 单轴、D 指纹替换、E 并列带宽、F Butina、G 确定性重跑 | 回答「名单结论是否依赖具体阈值/分数/指纹/并列次序」这一可复现性核心问题 |

| 6 MD 批次选择（15_md_batch_select，分析步） | 从主口径 50 名单挑 20 个进一级 MD：第二 hERG 模型（集成）不过门者**保留并逐分子标注模型依赖**（`retain_and_annotate`）、通用骨架重复者剔除 → 按 N=20 缩放配额（成对 Tc<0.65、骨架 ≤1、化学型 ≤8、风险档 ≤3）→ 在可行域内用 MILP（`scipy.optimize.milp` / HiGHS）按**字典序**求可证明最优：① max 分层覆盖下限 t（CNN_VS / cLogP / CNS MPO 三轴各 4 分位，资格池自分布）→ ② 在 t=t\* 内 max CNN_VS 跨度 → ③ max Σ CNN_VS（质量兜底）→ ④ 规范并列打破；v1.6.0 贪心填充降为对照实现留痕（实测它正是「max Σ CNN_VS」口径的全局最优解）；固定 2 位阴性对照 | 名单层门禁已完成，剩下的是机时分配，而**目标口径决定这批 MD 能回答什么问题**：分位覆盖与跨度是为让「对接分数 ⇒ 结合稳定性」可检验（纯分数截断会把跨度压到 0.5 以内、且最低分位 0 席），Σ CNN_VS 只作同信息量下的兜底。**「最优」限定在该目标函数口径内**（其 Σ CNN_VS 反而比质量口径低 1.16%、名册 6 换 6），低分位分子失败概率更高正是要**测量**的量；安全性问题由正交方法解决，不占用 MD 预算 |

## 3. 参考来源（实际检索/读取）

- Liu W, et al. *Science* 2012;337(6091):232–236. PMID 22798613; DOI 10.1126/science.1219218
  —— 4EIY 结构（https://www.rcsb.org/structure/4EIY ），本包受体来源。
- McNutt AT, et al. GNINA 1.0. *J Cheminform* 2021;13:43. PMID 34108002;
  DOI 10.1186/s13321-021-00522-2 —— 姿态打分来源（https://github.com/gnina/gnina ）。
- Trott O, Olson AJ. AutoDock Vina. *J Comput Chem* 2010;31(2):455–461. PMID 19499576
  —— 姿态内部 Vina 经验分来源。
- Wallach I, Heifets A. *J Cheminform* 2018;10:45. PMID 31430292
  —— DUD-E 诱饵偏置，故富集指标不作为本包验收判据。
- Wager TT, Hou X, Verhoest PR, Villalobos A. *ACS Chem Neurosci* 2010;1(6):435–449.
  DOI 10.1021/cn100008c —— CNS MPO Table 1 的分段函数与 Table 2/3 回归用例。
  本包 CNS MPO 的两个输入均为**上游第三方预测**（非本包重算），随上游表列固化：
  `basic_pka`＝MolGpKa 最强碱性中心 pKa（provider version `4dc835219baad18c2e705274ac54d49597703966`，
  **点值、无区间**；`pka_status='no_basic_center'` 者上游记 0，全表 192 个，按 pKa≤8 档取满分 T0=1.0）；
  `clogd74`＝cpLogD v2.0 / ChEMBL33（provider version `sha256:281338b9…24a72`，`clogd_confidence=0.8`）
  点估计，并附 **80% 预测区间** `clogd74_lower_80` / `clogd74_upper_80`。
  决策卡另给出区间界定的 CNS MPO：`cns_mpo_lo80`（悲观下界：cLogD 取区间上界）、
  `cns_mpo_hi80`（乐观上界：取区间下界）、`cns_mpo_pi80_span`，均为标注列、不作门。
- Swanson K, et al. ADMET-AI. *Bioinformatics* 2024;40(7):btae416.
  https://academic.oup.com/bioinformatics/article/40/7/btae416/7698030
  —— hERG 概率与 105 个 ADMET 端点的来源。本包以**外部预测表接口**接入该模型，未内置实现
  （`03_gates/herg/HERG_INTERFACE.md`）：① 落门表 `03_gates/herg/herg_predictions.csv`
  （四臂并集 3,236 个分子，按 hERG < 0.5 落门，由 HPC 应用 `admet-ai` 生成）；
  ② 全量 f3 库表 `03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv`（49,963 × 124，
  含 41 端点 + 8 描述符 + 105 个已上市药（DrugBank approved）百分位），供分析步 10b/10c
  做注释、阈值锚定与门作用域缺口量化。两表同源：交集 3,236 个分子上逐值一致（max|Δ| = 1e-6），
  因此**不构成独立验证**（独立性的量度见 08b 第二模型交叉验证）。
- RDKit 2025.09.5. https://www.rdkit.org —— 描述符、PAINS/Brenk catalog、QED、
  Bemis-Murcko 骨架、Morgan/MACCS 指纹、Butina 聚类。
- ChEMBL REST API（EBI），靶点 `CHEMBL240`（KCNH2 / hERG, Homo sapiens）：
  https://www.ebi.ac.uk/chembl/api/data/activity.json （检索日期 2026-09-23；
  第二模型交叉验证的训练/标注数据来源）

## 4. 复现性工程要点（踩坑记录）

| 问题 | 现象 | 处理 |
|---|---|---|
| SDF 块合法性 | 按 counts 切片写「数据部分」会丢掉 `M  END`，文件虽含 `$$$$` 但块非法（RDKit：`Cannot convert '> <' to int`） | 统一走 `molblock_data_lines()`，以 `M  END` 收尾；`tests/verify_pack.py` 增加「`$$$$` 数 == `M  END` 数」与几何抽样比对 |
| Butina 内存 | RDKit `Butina.ClusterData` 在 ≈5,000 分子池上把距离展成 n×n 矩阵被 OOM 杀掉 | 自实现 `cluster_butina_fast`：只建邻居表，其余步骤逐条对齐官方源码；在 n=200/800/1500 上划分与官方**完全一致** |
| SVG 字节抖动 | matplotlib 把生成日期与随机 hashsalt 写进 SVG，同一数据两次运行字节不同 | 固定 `svg.hashsalt` 且 `metadata={"Date": None}` |
| 比对自指 | 确定性比对报告记录上一轮差异数，自身必然变化 | 将 `determinism_compare.json`（及其交付副本）排除在哈希范围外 |
| 步序末位 | G 臂需覆盖交付层产物 | `12_determinism` 置于整链末步（`config/run_order.json`） |
| 自指闭环 | 报告/索引/验收/MANIFEST 把 G 臂结论写进自身，而结论又由含它们的哈希范围决定 → True/False 振荡 | 4 个「状态载体」+ 2 个比对报告移出比对范围；范围改为对生成封闭（125 个解析型产物） |
| 载体滞后一轮 | `MANIFEST.sha256` 在报告之前生成，记录的是上一轮文件哈希 | 判据取「连续两轮载体一致」；报告不嵌入 `mode` 等逐轮可变字段 |
| 注释列追加不得动成员 | 交付物需「追加列」而非旁挂表，但 pandas 重写会改变数值格式（`0.4864445628151945` → `0.486445`）与边界分子取值 | 追加时按 `dtype=str` 读入原文件、只写出注释列，原列逐字节保留（断言：`cid` 顺序与全部原列逐值不变）；注释值统一 6 位有效数字 |
| 同源表可能被误当独立验证 | f3 全量表与落门表的 hERG 同源（交集逐值一致），若并列陈述会造成「两模型互证」的错觉 | 在规则、报告、README 与 `REPRO_INFO.known_gaps` 四处显式声明「同源同值，不构成第二模型」，并给出交集 max\|Δ\| 作为证据 |
| 外部表路径值含注释 | `frozen_rule.json` 的 `external_table` 曾写成 `"…/herg_predictions.csv（列：cid,hERG）"`，`ROOT / 该值` 永不成立 → 即使表已就位仍判 `pending` | 路径字段只放纯路径，列名提示移入独立字段 `external_table_columns`；`tests/verify_pack.py` 增补 hERG 接口断言（列名、取值 ∈[0,1]、哈希一致） |
| C 臂 hERG 轴自证 | 单轴敏感性若直接按「门开启」重算，结果与主名单恒等（重叠 1.0），等于没测 | 改为**关掉 hERG 门重选一遍**（池 1,960）再与主名单（池 606）比：重叠 0.2292，方为可实测结论 |
| 缺失值显示 | 应用接口不回传权重哈希 → 空列被 pandas 读成 `NaN`，报告里渲染为 `nan…` | 接口读取时把空/`NaN` 归一为 `None`，报告统一写作「未提供（应用接口不回传权重哈希）」 |
| 重复列静默丢列 | `descriptors.csv` 已含 `basic_pka`/`clogd74`/`pka_status`，`build_master` 又合并一次上游同名表 → 产出 `basic_pka_x/_y`；下游按「列存在才写入」取数，交付决策卡里 **静默缺少 CNS MPO 的两个输入列** | 合并前先删除描述符表中的上游同名列（上游表为唯一来源）；`verify_pack` 增加「主表无 `_x`/`_y` 重复列」与决策卡列存在性断言 |
| 区间中点当点估计 | `cns_mpo_bounded` 初版用 80% 区间中点作 cLogD 点估计，虽上下界正确但会改变点估计值（与上游 `clogd74` 口径不一致） | 点估计一律取上游 `clogd74` 点值，区间只用于 `cns_mpo_lo80`/`hi80`；回归用例（4.6179 / 4.8483 / 5.1838）逐位复现 |
| ChEMBL 关系符 URL 编码 | 查询 `standard_relation='='` 时写成 `%27%3D%27`（带引号）返回 **0 条**，被误读成「该靶点无 IC50 数据」 | 关系符必须编码为 `%3D`；首次抓取即比对 `api_total_count`（IC50 11,709 / Ki 1,812）与落盘行数（13,521）是否吻合 |
| 标注阈值解包写反 | `lo, hi = blocker_ge, nonblocker_le` 使 `label = (p ≥ 5 → 1, p ≤ 6 → 0)` **恒为 1**（阻断 7,513 / 非阻断 3,415 / 模糊带 0），模型"性能"看似正常但语义全错 | 显式 `hi = blocker_ge`、`lo = nonblocker_le`，并加 `assert hi > lo`；同时校验「模糊带规模 > 0」与两类样本量级 |
| 混淆矩阵语义错标 | `confusion_matrix(a, b, labels=[0,1]).ravel()` 的顺序是 `TN, FP, FN, TP`，把「模型1过/模型2不过」记成 `FP` → 结论方向相反 | 固定写 `m1_pass_m2_fail = fn`、`m1_fail_m2_pass = fp`，并在验收里断言「四格和 == 池规模」以防再次错位 |
| 标注集规模被清洗暗削 | 13,521 条原始记录 → InChIKey 去重 10,928 个分子 → 可用标注仅 **6,057**（模糊带 (5,6) 丢弃 4,871），若只看原始行数会高估训练集 2 倍以上 | 报告链路上每一步的计数（raw → dedup → labeled → blocker/nonblocker）并写进结论 JSON；`verify_pack` 断言 `blocker + nonblocker == labeled` |
| 跨度线性化退化（假跨度） | 跨度若只用单侧下界式（`u ≥ z_i·x_i`），最大化 `u − l` 会把 `u` 顶到上界、`l` 压到下界，得到与所选名册无关的跨度（首轮报 0.7337，名册实算仅 0.6903） | 改用**独热指派**：`y_i ≤ x_i`、`Σy=1`、`u = Σ z_i y_i`，再加紧约束 `u ≥ z_j − M(1−x_j)`（`l` 同理取下界方向）；脚本内置**自洽断言**（各阶段目标值必须等于名册实算值），自检 [21] 再从批次表独立复算 t / 跨度 / Σ |
| MILP 并列解的确定性 | 同为最优的多个名册会被求解器内部 tie-break 任选，破坏「确定性重跑」判据 | 加第 4 阶段**规范并列打破**（`max Σ w_i·x_i`，`w_i` 按 CNN_VS 降序索引）；同一步内重解一次并断言名册逐位相同（写入证书 `resolved_roster_consistent`），自检 [21] 另做**独立重解**断言 `t` 上界 = `t*` |
| 目标值随上游静默漂移 | 覆盖下限 / 跨度 / Σ CNN_VS / 与对照的互换数会随名单或描述符变化，而脚本只报数字不报错 | 把实测值登记进规则 `step6_md_batch.expected_certificate`；步骤脚本与 `tests/verify_pack.py` 第 [21] 组**双向断言**（t\*、跨度、Σ、互换数、分位边界容差 1e-6），任一处漂移即显式失败并要求人工复核 |
| 整目录删除后无法重跑 | 步骤脚本用 `Path.write_text` 写文件时不建父目录：只删文件无碍，但把产物**目录**整体删掉后会在中段抛 `FileNotFoundError`（实测 `11_robustness` → `05_robustness/reports/`） | `run_all.py` 起链前按 `config/run_order.json` 的 `outputs` 预建全部产物父目录（只建目录、不生成占位文件，因而不改变任何产物哈希），并把补齐清单写进 `logs/run_all.log`；`tests/verify_pack.py` 第 [20] 组用**功能测试**守住该行为 |
| 冻结输入层被连带删除 | `03_gates/admet/raw/`（ADMET 全量表）与 `03_gates/herg/cross_model/raw/`（ChEMBL 活性表 + 抓取元信息）位于产物目录内部，但属**输入层**：随产物一起删掉后，失败表现为深处的裸 `FileNotFoundError` | `08b` 增加前置检查（`SystemExit`），`10_select_lists` 经 `vendor/master.py` 的缺输入检查（`FileNotFoundError`）一并补上「属输入层、请从输入副本恢复（不重抓/不重算）」的指引；`verify_pack` [20] 组断言 raw 层齐备 |

## 5. 与上游工程的关系

本包把上一轮多路线工作中分散的口径（姿态闸门、锚点门、类药性门、选序约束）**重新实现并
自包含化**：不复制、不引用任何包外产物；所有数值（包括参考活性集统计、阴性对照的选取）
都在包内从原始输入与归档重算，因此可以与既有数字逐项对账（见最终报告 §8）。
