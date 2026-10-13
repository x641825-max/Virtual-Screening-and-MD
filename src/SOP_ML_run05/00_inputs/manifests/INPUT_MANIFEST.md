# 输入固化清单（SOP_ML）

- 规则集：`SOP-05-A2A-antagonist-shortlist` v1.7.0
- rules_sha256：`8e71187653eb21269695d4e80bb1506002db8edb729c5a14e1e9acf5fec63881`
- 全部输入哈希一致：**是**

| 文件 | 角色 | md5 | 校验 | 大小(MB) |
|---|---|---|---|---|
| `00_inputs/raw/inputs/top5000_best_pose.sdf` | 4,979 个 gnina 最优姿态（run05 · ML 生成分子 · CNN 口径工作池） | `0ac1353fa852fd49ce2f34d1e4e148a0` | 一致 | 13.63 |
| `00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv` | 全库 CNN_VS 排名 + 上游描述符 + MolGpKa 碱性 pKa + cpLogD7.4 | `d9c600d9cbaf5019b359cbc7ce128a92` | 一致 | 2.01 |
| `00_inputs/raw/archives/run01_redock_enrichment.zip` | run01 归档：受体 4EIY、5,098 带标签富集库、富集基线、重对接验证 | `82a64a41b06f49dc4327b6bd7fcb135e` | 一致 | 26.15 |
| `00_inputs/raw/archives/screening_kit.zip` | 筛选套件：ChEMBL 独立已知配体集、模型卡、受体副本 | `ff5902239a51536d6bc34e42599ceac3` | 一致 | 66.87 |
| `00_inputs/raw/archives/redocking.zip` | 3EML/SOP-02 重对接与诱饵基准（协议溯源与对照物来源） | `cac648447c54f688562ddafed90dcf65` | 一致 | 39.95 |
| `00_inputs/extracted/run01/01_redocking/receptor/4EIY_protein.pdb` | 派生输入（解包产物） | `aef775ee83735af8a878d42ca1e5fa85` | 一致 | 0.25 |

## 关键派生输入

- 5,098 带标签富集库（479 活性 / 4,619 诱饵）：`00_inputs/extracted/run01/02_enrichment_5098/library/quick_5100_library.csv`（5100 行）
- ChEMBL 独立已知 A2A 配体：`00_inputs/extracted/screening_kit/07_independent_validation/independent_a2a_chembl.csv`（9130 行）

> 本包内所有引用均位于包根之下；脚本禁止引用包外路径（`tests/check_isolation.py` 静态扫描）。
