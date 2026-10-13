# 数据来源、版本、许可与获取方式

按赛事《代码提交要求》附件 5 §三（可复现性与溯源）与 §五（特殊情形说明）编写。
**本包不分发上游原始数据**：下表逐项登记其来源、版本、获取时间与 sha256/md5，供核验与恢复。

---

## 1. 上游原始输入（未随包分发，恢复后可跑完整链）

### 1.1 目录化合物链（`src/SOP_5`）

| 文件（恢复位置） | 内容 | md5 | sha256 |
|---|---|---|---|
| `00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv` | 全库 CNN_VS 排名 + 上游描述符 + MolGpKa pKa + cpLogD7.4（49,963 行） | `444f07fbe92cd15572ea459e6e955c66` | `dae23885e2748ca8c024c8847e01f74c560527cdc80f527884f8c2d35efb89c4` |
| `00_inputs/raw/inputs/top5000_best_pose.sdf` | 5,000 个 gnina 最优姿态（每分子 1 条三维记录） | `0ed679e1d61d7f562b27bfcdd1c7e956` | `1040e55a80d96139eb9538871050849015551dfbbdcb831b8c417991406353de` |
| `00_inputs/raw/archives/run01_redock_enrichment.zip` | 受体 4EIY、5,098 带标签富集库（479 活性）、重对接与富集基线 | `82a64a41b06f49dc4327b6bd7fcb135e` | `dc439580ea71d1320a36969a9618e4685fd91201bbcebda447b33a291b2d9fd9` |
| `00_inputs/raw/archives/screening_kit.zip` | ChEMBL 独立已知 A2A 配体集、模型卡、受体副本 | `ff5902239a51536d6bc34e42599ceac3` | `21e859f7a03c5693ca9b9595bf769fcfe20297d05ba43ea16b486f2a64182e6e` |
| `00_inputs/raw/archives/redocking.zip` | 3EML / SOP-02 重对接与诱饵基准（协议溯源） | `cac648447c54f688562ddafed90dcf65` | `562d40fecb692a1aefcd2018808eb3b9f41815b199566cf5a97dae09212a1ea8` |
| `03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv` | ADMET-AI f3 库冻结全量表（49,963 × 124，81.1 MB） | — | `5d693ea7d9d06e29375f7b73fc496911ce514444ef2e243e3b2feb4703a67476` |

上述 md5/sha256 同时登记在 `src/SOP_5/00_inputs/INPUTS_MD5.txt`、`INPUTS_SHA256.txt`
与 `00_inputs/manifests/`；恢复后可用包内脚本对账：

```bash
cd src/SOP_5 && python3 00_inputs/scripts/01_verify_inputs.py
```

### 1.2 ML 生成分子链（`src/SOP_ML_run05`）

该链与 1.1 **共用同一套脚本、阈值与门，只换一个输入表**（数据对象换成本批生成分子的对接结果）：

| 文件（恢复位置） | 内容 | md5 | sha256 |
|---|---|---|---|
| `00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv` | **换表对象**：4,979 个 ML 生成分子（SCALE-）的排名 + 上游描述符 + pKa + cLogD74 | `d9c600d9cbaf5019b359cbc7ce128a92` | `f7ed93c2f4e7ae0fc43fad8d1487d38807070b9d0a6fc03576a59cffc201d9e2` |
| `00_inputs/raw/inputs/top5000_best_pose.sdf` | 同批 4,979 个分子的 gnina 最优姿态 | `0ac1353fa852fd49ce2f34d1e4e148a0` | `f5d8f2075626f060ce048b0dbe44b4e48dff270a4bd4d95aed98598d316842eb` |
| `00_inputs/raw/archives/*.zip` | 与 1.1 三个归档同文件同名同哈希（受体、基准、ChEMBL 独立集） | 同 1.1 | 同 1.1 |

> 注意：`src/SOP_ML_run05/00_inputs/INPUTS_MD5.txt` 内登记的仍是**换表前**的 SOP_5 值
> （骨架生成时的快照）。换表后的新哈希以本表为准，登记流程见该包 `SWAP_SPEC.md §5`。

---

## 2. 第三方模型与外部数据（预测层，本包不重算）

| 项 | 来源 | 版本 / 获取时间 | 在本项目中的角色 |
|---|---|---|---|
| **gnina** | McNutt A, et al. 2021, PMID [34108002](https://pubmed.ncbi.nlm.nih.gov/34108002/)；[github.com/gnina/gnina](https://github.com/gnina/gnina) | 1.3.2；`semi_flexible` + `cnn_scoring=rescore`，`exhaustiveness=8`，`num_modes=9` | 分子对接与 CNN 打分，产出 `CNN_VS = CNNscore × CNNaffinity`（**关键预测指标**） |
| **ADMET-AI**（chemprop） | Swanson K, et al. *Bioinformatics* 2024, PMID [38420927](https://pubmed.ncbi.nlm.nih.gov/38420927/)；TDC ADMET 基准 | 冻结预测表（HPC 应用 `admet-ai`），2026-09 | hERG 分类门（<0.5）与 ADMET 注释列 |
| **ChEMBL REST API** | EMBL-EBI，`target_chembl_id=CHEMBL240`（KCNH2 / hERG，Homo sapiens） | 抓取时刻 2026-09-23T06:47:45Z；原始记录 13,521 条，可用标签 6,057（blocker 2,411 / non-blocker 3,646） | 第二 hERG 模型的训练/标注集（RF 5 折 AUC 0.9489） |
| **MolGpKa** | 上游第三方 pKa 预测 | provider version `4dc835219baad18c2e705274ac54d49597703966` | CNS MPO 的 `basic_pka` 输入 |
| **cpLogD v2.0** | 上游第三方 logD 预测 | 随上游排名表冻结 | CNS MPO 的 `clogd74` 输入（附 80% 区间） |
| **A2A 便携初筛模型（ECFP4-RF）** | 本项目上游自建，训练数据见 `models/MODEL_CARD.md` | 训练日期 2026-09-14 | 在大规模库上做排序富集（**不在本包的选批段重跑**） |
| **DUD-E / DeepCoy 诱饵** | 公开基准，用于富集检验与阴性对照 | — | 2 个共用阴性对照（NEG1/NEG2）与富集评价参考 |

第三方训练/基准数据的许可与引用在 `models/MODEL_CARD.md` 逐项列出。

---

## 3. 许可与再分发边界

- ZINC / MCULE 目录化合物与其对接打分、ADMET 预测表体量大且属第三方分发，**随包不再分发**；
  评审方如需完整复核，可按 §1 的哈希向上游索取对应文件。
- ChEMBL 数据按 EMBL-EBI 的 ChEMBL 许可使用（CC BY-SA 3.0 系）；本包只分发由它派生的
  **模型交叉验证结果表**（`src/*/03_gates/herg/cross_model/`），不分发原始活性记录以外的库表。
- `redocking.zip` / `screening_kit.zip` / `run01_redock_enrichment.zip` 含本项目上游的受体副本与
  基准构建产物，本包不分发，但保留完整哈希与解包白名单（`src/SOP_5/00_inputs/scripts/00_unpack_inputs.py`）。

---

## 4. 随包可运行的最小示例

`data/example/` 提供小规模示例数据（可直接跑通 Notebook 与冒烟测试），
`src/` 下随包保留了「MD 候选选择 → 联立合并」所需的全部中间表（`04_selection/`、`03_gates/`、
`02_druglikeness/tables/fingerprints_ecfp4.npz`、`09_synth_access/`），因此该段**无需任何上游大文件
即可完整复跑**（见包根 `README.md` §6 与 `logs/RUN_RECORD.md`）。

---

## 5. 随包分发的自有数据（M2 初筛模型，2026-10-03 补充）

| 文件 | 内容 | 规模 | 字段 |
|---|---|---|---|
| `data/dataset_train.csv` | M2 训练集（RDKit 规范化 SMILES + 标签 + 来源标记） | **1,307 行**（活性 479 / 阴性 828） | `cid`, `smiles`, `label`, `origin` |
| `data/dataset_external.csv` | M2 外部独立验证集（训练中未见过的化学型 + held-out 诱饵） | **2,167 行**（活性 167 / 阴性 2,000） | `cid`, `smiles`, `label`, `origin` |

配套文件：训练入口 `train.py`（依赖 `src/features.py`、`src/metrics.py`）、训练日志 `logs/train.log`、
导出权重 `models/a2a_rf_portable.npz` / `models/a2a_rf_trained_portable.npz` / `models/a2a_rf_v2_portable.npz`、
模型卡 `models/MODEL_CARD.md` 与 `models/MODEL_CARD_A2AR-RF-v2.md`。

来源与边界：

- 训练集分子来自公开基准（DUD-E 风格诱饵与文献活性/非活性记录，出处见 `models/MODEL_CARD.md` §2），
  本包只分发 SMILES 与标签，不含第三方库表；
- 外部验证集为训练中未出现过的化学型与 held-out 诱饵，用于跨骨架泛化评价；
- 训练数据的划分口径（骨架切分、去重与防泄漏）见 `models/MODEL_CARD.md` §2 与 §4。
