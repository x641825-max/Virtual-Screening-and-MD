# SOP_ML 输入换版接口契约（SWAP_SPEC）

> 版本：1.0（2026-09-24）｜适用对象：`SOP_ML` 骨架（派生自 `SOP_5`，规则 v1.7.0 一字未改）
> 本文件是「换数据」的操作口径与验收判据。链上任何一步的失败都不应绕过本文件的登记要求。

## 1. 一句话口径

本包只换**一个**输入文件：全库打分/排名表。位置与文件名固定，35 列 schema 不变，其余 4 个输入（5,000 姿态、三个归档）与全部阈值/门沿用 v1.7.0。换表后必须同步处置 §4 列出的 6 项耦合，否则链会以「显式失败」或「静默丢行」两种方式之一出问题。

| 项目 | 值 |
|---|---|
| 目标路径 | `00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv` |
| 文件名 | **不得更改**（11 个脚本按字面量读取该文件名，改名等于改代码） |
| 当前登记（v1.7.0） | md5 `444f07fbe92cd15572ea459e6e955c66`｜sha256 `dae23885e2748ca8c024c8847e01f74c560527cdc80f527884f8c2d35efb89c4` |
| 当前规模（实测） | 49,963 行 × 35 列｜18.75 MB |
| 登记位置 | `config/frozen_rule.json → inputs.<文件名>`（md5 / sha256 / role） |
| 校验工具 | `tests/check_swap_csv.py`（结构 + 取值 + 耦合覆盖） |
| 登记工具 | `00_inputs/scripts/00b_register_swap_input.py` |

## 2. 谁在读这张表（列级消费者）

| 脚本 | 用到的列 | 用途与硬校验 |
|---|---|---|
| `01_pose_qc/scripts/02_pose_identity_qc.py` | cid, smiles, CNN_VS, CNNscore, CNNaffinity, CNNaffinity_variance, minimizedAffinity, rank_score | **硬校验**：SDF 属性 vs 上游表逐值比对；`\|CNNscore × CNNaffinity − CNN_VS\| ≤ 5e-4` |
| `02_druglikeness/scripts/05_descriptors.py` | cid, CNN_VS, mw, logp, hbd, hba, tpsa, rotb, basic_pka, clogd74, pka_status, clogd74_lower_80, clogd74_upper_80, pka_basic_center_count, clogd_status, clogd_confidence, pka/clogd provider 与 version | RDKit 重算描述符并与上游逐值交叉校验；CNS MPO 的两个输入取 `basic_pka` / `clogd74`（含 80% 区间） |
| `03_gates/scripts/10b0_admet_regenerate.py` | cid, smiles | 生成 admet-ai 应用输入（手工步，不入链） |
| `03_gates/scripts/10b_admet_annotation.py` | cid, CNN_VS | 池级 ADMET 注释与值一致性交叉核验 |
| `04_selection/scripts/10_select_lists.py` | cid, smiles, CNN_VS, qed, molstar | 选序主键：CNN_VS 降序 + 成对 Tanimoto < 0.65 约束贪心 |
| `06_deliverables/scripts/14_report.py` | CNN_VS 区间、`\|CNNscore × CNNaffinity − CNN_VS\|` 越界数、basic_pka、clogd74 及区间、provider | 报告与决策卡 |
| `07_md_batch/scripts/15_md_batch_select.py` | cid, CNN_VS | MILP 目标（分层覆盖 / CNN_VS 跨度 / Σ CNN_VS / 加权序） |
| `vendor/master.py`、`vendor/chem_core.py` | 13 个上游列（见 §3 标 U） | 主表与决策卡的列来源 |
| `tests/verify_pack.py` | 多个 | 自检：pKa/cLogD 列齐备、区间自洽、决策卡字段 |
| `08_visuals/scripts/make_visuals.py` | 多个 | 图册读数 |

## 3. 35 列契约（实测口径，2026-09-24 复核）

`类型` 列：num = 数值（下游做算术），str = 字符串。`空` 列给出 SOP_5 实测空值数，`U` 标记进入 `vendor/master.py` 的 UPSTREAM_COLS（13 列随主表进入交付层）。

| # | 列名 | 类型 | 实测范围 / 取值 | 空 | 换表要求 |
|---|---|---|---|---|---|
| 1 | rank | num | 1 … N，连续无重 | 0 | 必须是 1..N 的排列 |
| 2 | cid | str | 49,963 个唯一（MCULE-/ZINC-…） | 0 | 非空、唯一；下游所有 join 的主键 |
| 3 | smiles | str | 唯一，RDKit 可解析 | 0 | 非空、可解析；重复会在 10b0 告警 |
| 4 | CNN_VS | num | 0.6619 … 8.2632 | 0 | 与第 5、6 列满足恒等式（容差 5e-4） |
| 5 | CNNscore | num | 0.1060 … 0.9844 | 0 | 同上 |
| 6 | CNNaffinity | num | 3.8559 … 9.0029 | 0 | 同上 |
| 7 | CNNaffinity_variance | num | 0 … 2.0756 | 0 | 进入姿态属性比对 |
| 8 | minimizedAffinity | num | −13.9670 … 125.5629 | 0 | 进入姿态属性比对（Vina 口径） |
| 9 | part | str | `part_01` … `part_32` | 0 | 命名建议保持 `part_NN` |
| 10 | mw | num | 180.23 … 649.80 | 0 | 与 RDKit 重算值交叉校验 |
| 11 | logp | num | −3.0 … 6.0 | 0 | 同上（下游用重算值，上游列只做对账） |
| 12 | hbd | num | 0 … 5 | 0 | 同上 |
| 13 | hba | num | 2 … 12 | 0 | 同上 |
| 14 | tpsa | num | 26.0 … 209.9 | 0 | 同上 |
| 15 | rotb | num | 0 … 12 | 0 | 同上 |
| 16 | pains | num | **整列全空** | 49,963 | 允许空：pains 由 RDKit 重算 |
| 17 | brenk | str | 54 种告警名（`imine_1`…） | 41,764 | 允许部分空；下游用重算的整型列 |
| 18 | p_RF | num | 0.5806 … 1.0 | 0 | 与 rank_score 同值（SOP_5 实测） |
| 19 | rank_score | num | 0.5806 … 1.0 | 0 | 选序/报告只读 CNN_VS 与它 |
| 20 | basic_pka | num | 0 … 17.6044 | 0 | **CNS MPO 输入**；无碱性中心时为 0 且 pka_status=no_basic_center |
| 21 | clogd74 | num | −6.032 … 7.366 | 0 | **CNS MPO 输入**；必须落在 32/33 两列之间 |
| 22 | pka_provider | str | `MolGpKa` | 0 | 非空；换 provider 需在报告口径中登记 |
| 23 | pka_provider_version | str | commit hash 形态 | 0 | 非空 |
| 24 | pka_status | str | `ok`(49,771) / `no_basic_center`(192) | 0 | **取值域受限**：其它取值会被下游当作 imputed_0.5 静默处理 |
| 25 | pka_basic_center_count | num | 0 … 11 | 0 | 与第 26、27 列长度一致 |
| 26 | pka_basic_center_atom_indices | str | `"[3,17,19,21]"` | 0 | 非空 |
| 27 | pka_basic_center_values | str | `"[5.6605,2.9297,…]"` | 0 | 非空 |
| 28 | pka_error | num | **整列全空** | 49,963 | 允许空（诊断列） |
| 29 | clogd_provider | str | `cpLogD_v2.0_ChEMBL33` | 0 | 非空 |
| 30 | clogd_provider_version | str | `sha256:281338b9…` | 0 | 非空 |
| 31 | clogd_confidence | num | 0.8（常量） | 0 | 建议 ∈ (0,1] |
| 32 | clogd74_lower_80 | num | −7.029 … 7.197 | 0 | 与第 21 列构成区间（≤ 点估计） |
| 33 | clogd74_upper_80 | num | −5.589 … 7.569 | 0 | ≥ 点估计 |
| 34 | clogd_status | str | `ok` | 0 | 非空 |
| 35 | clogd_error | num | **整列全空** | 49,963 | 允许空（诊断列） |

> 上表「实测范围」是 SOP_5 输入层的事实记录，**不是**对新表的取值要求——新模型的分数量程不同属正常。列为 `U` 的 13 列会随主表进入交付层，语义变更需在报告里写明。

## 4. 换表后必须处置的 6 项耦合（按顺序）

### 4.1 姿态 SDF 的内嵌属性（最容易踩）

`00_inputs/raw/inputs/top5000_best_pose.sdf` 的每个分子自带属性：`cid, CNN_VS, CNNscore, CNNaffinity, CNNaffinity_variance, minimizedAffinity, rank_score, n_fragments`。第 1 步会把它们与排名表**逐值比对**。实测：在当前 SOP_5 数据上两侧逐值一致（校验项 E13 通过）。因此若新表来自不同打分定义，第 1 步会直接判「属性不一致」。三条处置路径，必须选一条并登记：

| 处置 | 做法 | 代价 |
|---|---|---|
| A（正解，超出本包范围） | 用新模型重跑对接，重生成 5,000 姿态及其属性 | 需上游 gnina/CNN 打分作业 |
| B（换选择依据，保留旧姿态） | 保留旧 SDF 属性，把新表当「排序依据」层，在 `README` 与报告登记「姿态属性来自旧打分轮」 | 第 1 步属性门会失败，需在规则里显式登记例外（属口径变更，rules_sha256 会变） |
| C（同名换源） | 新表沿用同名同义列但数值来自新模型，并同步重写 SDF 属性 | 需自己保证两侧同源，否则等于把不一致藏进输入 |

> **实测提醒（2026-09-24）**：如果新表只是「旧表 + 新增分子」（例如工作区里的 `CNN_VS_ranking_merged.csv`），
> 那么 5,000 姿态的 cid 在表内、且内嵌属性与表值逐值同源 → **本节无需任何动作**。
> 只有当新表给**同一批 cid** 换了新分数、或姿态被重生成时，本节才需要处置。实测读数见 §8。

### 4.2 「top-5,000」的选择来源变了

5,000 姿态是按旧 CNN_VS 取头部得到的。新表重排后，头部 5,000 的分子集合会变——此时「姿态集」不再对应新表的 top-5,000。若这层语义重要，同 4.1 的处置 A；若只是「固定一批姿态 + 新排序」，须在文档写明。

### 4.3 hERG 落门表必须 ⊆ 新表 cid

`03_gates/herg/herg_predictions.csv`（实测 3,236 行 × 4 列）的 cid 在 SOP_5 数据上全部落于排名表内。第 3 步按 cid 左连接，**落在表外的行会被静默丢弃**，门决策随之改变。校验项 E14 会 FAIL 拦住这种情况。

### 4.4 ADMET 冻结表 cid 集合必须一致

`03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv`（49,963 行 × 124 列，81 MB）在 SOP_5 数据上与排名表 cid 集合**完全一致**（实测双向差集为 0）。它是外部 admet-ai 作业的产物，不在包内重算。新表 cid 集合一旦变化，必须重生成（`10b0_admet_regenerate.py` 手工执行），否则 10b/10c 的覆盖率下降。校验项 E15 会 WARN。

### 4.5 规则登记与 `rules_sha256`

`config/frozen_rule.json → inputs` 登记了该表的 md5/sha256；`01_verify_inputs` 对不上就会 **exit 4**（显式失败，不静默降级）。因此换表必须走登记：`00b` 工具会更新 `inputs` 登记、把 `rules_sha256` 置为 `PENDING`，下一次 `01_verify_inputs` 重新冻结。说明：

- **阈值与门一字未改**（rule_version 仍为 1.7.0）；变的只是「输入登记」这一处，因此 `rules_sha256` 会变——这是「输入换版」的预期痕迹，不是规则被改。
- 影响面：`PACK_MANIFEST.sha256` 需重生成（`tests/make_index.py`），`FREEZE.json` 由 `01_verify_inputs` 重建。

### 4.6 确定性基线必须重取

`05_robustness/determinism/` 的 G 臂基线绑定旧数据。换表后须走既定节奏：`run_all.py --freeze`（写基线，**不构成通过**）→ 再连续跑 **2–3 轮**普通比对直到 7 个载体逐字节一致。机制：第 17 步 `13_assemble` 早于第 19 步 `12_determinism`，交付层副本永远落后一轮。

> 不受影响：5,098 富集库、诱饵基准、3EML 重对接（三个归档）与 2 个阴性对照（NEG1/NEG2 由 5,098 库按规则重算，与排名表无关）。

## 5. 接入流程（4 步）

```bash
cd SOP_ML

# 1) 放入新表（文件名必须一致；建议同时保留原表以便对比）
#    00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv

# 2) 契约校验（先只查表；有 SDF/hERG/ADMET 时一并带上做耦合覆盖）
python3 tests/check_swap_csv.py --file 00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv \
    --json 00_inputs/manifests/swap_check.json
#    期望：结论「通过」；FAIL 必须清零后再继续

# 3) 登记（先 dry-run 看将改哪里，--confirm 才落盘）
python3 00_inputs/scripts/00b_register_swap_input.py \
    --file 00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv --confirm

# 4) 重组索引 + 重取基线 + 全链
python3 tests/make_index.py            # 重生成 INDEX_PACK.md 与 PACK_MANIFEST.sha256
python3 run_all.py                     # 首轮即完成 rules_sha256 重新冻结 + 全链
python3 run_all.py --freeze            # 写确定性基线（不构成通过）
python3 run_all.py                     # 连续 2–3 轮普通比对 → 7 载体逐字节一致
```

## 6. 校验判据清单（`tests/check_swap_csv.py`）

| 编号 | 判据 | 级别 | 说明 |
|---|---|---|---|
| E1 | 行数 ≥ 5000 | FAIL | 5,000 姿态的主人必须能在表中找到 |
| E2 / E2b | 35 列齐备 / 无多余列 | FAIL / WARN | 列名逐一匹配（顺序不限） |
| E3 | cid 非空且唯一 | FAIL | 全部 join 的主键 |
| E4 / E4b | smiles 全部 RDKit 可解析 / 无重复 | FAIL / WARN | 重复会让 10b0 告警 |
| E5 | rank 为 1..N 的排列 | FAIL | |
| E6 | 20 个必需数值列非空可解析 | FAIL | |
| E7 | `\|CNNscore×CNNaffinity − CNN_VS\| ≤ 5e-4` | FAIL | 与第 1 步同容差；可用 `--allow-identity-break` 降级（须登记） |
| E8 | `clogd74_lower_80 ≤ clogd74 ≤ clogd74_upper_80` | FAIL | |
| E9 | `pka_status ∈ {ok, no_basic_center}` | FAIL | 其它取值会被下游静默按 imputed_0.5 处理 |
| E10 / E10b | `clogd_confidence ∈ (0,1]` / provider 与 status 非空 | WARN / FAIL | |
| E11 | part 命名 `part_NN` | WARN | |
| E12 | SDF 的 cid ⊆ 表 cid | FAIL | 需 `--sdf` |
| E13 | SDF 内嵌属性 vs 新表同源 | WARN | 需 `--sdf`；告警即触发 §4.1 的三选一 |
| E14 | hERG 落门表 cid ⊆ 表 cid | FAIL | 需 `--herg`；否则第 3 步静默丢行 |
| E15 | ADMET 冻结表 cid 集合一致 | WARN | 需 `--admet`；不一致即需重生成 |

退出码：`0` 通过｜`7` 有 FAIL（或 `--strict` 下有 WARN）｜`8` 文件缺失/不可读。
正反验证：`<构建期工作区>/SOP_ML_build/test_checker.py`（6 个受控缺陷全部命中、真实表 16 项通过 0 FAIL，2026-09-24 实测）。

## 7. 边界与不保证

- 本契约只管**输入层契约**：表能进得来、能过第 1–3 步的列级校验、不会静默丢行。它**不保证**新表的科学有效性（打分是否可靠、是否适用同一靶点与化学空间）——那是换表方的责任。
- 阈值与门沿用的前提是「同一靶点、同一库、同一批姿态」。若换表同时换了靶点或库，规则阈值、hERG 作用域、ADMET 冻结表都需重新论证，**不属于**本契约的适用范围。
- 骨架当前未包含大体积输入与链上产物，见 `RESTORE_LIST.md`；缺失状态下跑链会在相应步骤显式失败（设计如此，不静默降级）。


## 8. 实测记录（2026-09-24，针对工作区内的两份候选新表）

工作区里出现过两份候选表（`<构建期工作区>/SOP_ML_build/`，非本包内容，未改动）：

| 判据 | `CNN_VS_ranking_4979.csv` | `CNN_VS_ranking_merged.csv` |
|---|---|---|
| 规模 | 4,979 × 20 | 54,942 × 22（含 `source` / `source_rank`） |
| E1 行数 ≥ 5000 | **FAIL**（4,979） | PASS（54,942） |
| E2 35 列齐备 | **FAIL**（缺 15 列 pKa/cLogD 块） | **FAIL**（同样缺 15 列） |
| E2b 多余列 | WARN（`n_poses`） | WARN（`source` / `source_rank` / `n_poses`） |
| E3 cid 唯一 | PASS（4,979） | PASS（54,942） |
| E4/E4b smiles | PASS / PASS | PASS / WARN（10 处重复） |
| E5 rank 1..N | PASS | PASS |
| E7 恒等式 ≤5e-4 | 4.74e-4 ✓ | 4.95e-4 ✓（贴容差） |
| E12 姿态 cid ⊆ 表 | **FAIL**（0/5,000 命中） | PASS（5,000/5,000 命中） |
| E13 姿态属性同源 | 不适用 | **逐值一致**：6 个属性 max\|Δ\|=5.0e-5（= CSV 四位小数舍入），越界 0 行 |
| E14 hERG 落门表 ⊆ 表 | 未测 | PASS（3,236/3,236 命中） |
| 与 SOP_5 表 cid 关系 | 交集 0（全为新的 `SCALE-*`） | 完全包含原 49,963，另加 4,979 个新 cid |

由此得到两条可操作的结论：

1. **§4.1（姿态属性）在本例中无需动作**——`merged` 表对那 5,000 个 cid 携带的是与姿态 SDF 完全一致的分数
   （第 1 步的 `prop_max_abs_diff_vs_upstream` 会保持 5.0e-5 量级，与 SOP_5 定稿轮相同）。
2. **真正的阻塞项是 pKa/cLogD 缺列（E2）**，以及——若改用 4,979 那份——姿态覆盖（E12）与行数（E1）。
   缺列的两个可选补法：从 SOP_5 表按 cid 回填（命中 49,963 行，4,979 个新分子留空并登记），
   或由上游重算 MolGpKa / cpLogD 两列后一并给出。

### 8.1 缺列对下游的直接影响（按脚本）

| 脚本 | 依赖的缺失列 | 后果 |
|---|---|---|
| `02_druglikeness/scripts/05_descriptors.py` | basic_pka, clogd74, pka_status, clogd74_lower_80/upper_80, pka/clogd provider+version | CNS MPO 输入与 pKa/cLogD 标注列无法生成；代码会取到 `NaN` |
| `02_druglikeness/scripts/06_funnel.py` | clogd74 区间 | 漏斗统计「cLogD 区间宽度中位」缺输入 |
| `vendor/master.py` / `vendor/chem_core.py` | 13 个上游列 | 主表与决策卡缺这 13 列 |
| `06_deliverables/scripts/14_report.py` | basic_pka / clogd74 / 区间 / provider | 报告的 pKa/cLogD 段落与决策卡字段缺失 |
| `tests/verify_pack.py` 第 [11] 组 | 同上 | 自检「CNS MPO 输入与 80% 区间标注列」直接 FAIL |

即：**补列不是可选项**，除非把 CNS MPO 与 pKa/cLogD 标注整体降级（那属于口径变更，需升 `rule_version` 并重设验收）。
