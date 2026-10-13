# 路径映射与隔离约定（SOP_5）

## 1. 包根定位

每个脚本以三行自举定位包根（向上查找 `REPRO_INFO.json`）：

```python
_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p
```

`vendor/pack_paths.py` 的 `pack_root()` / `pth()` / `guard_inside()` 提供同一逻辑，
并对越界路径直接抛错。

## 2. 硬约定

1. 脚本只能通过 `ROOT / "<包内相对路径>"` 取文件；**禁止**绝对路径与 `../` 上溯。
2. 临时文件（若需要）放在 `logs/` 或对应步骤目录，不写 `/tmp`（沙箱每次调用会清空）。
3. `tests/check_isolation.py` 静态扫描全部 `.py`，命中禁例即退出码 6。
4. 端口/外部服务：本包不访问网络；hERG 门通过包内文件接口（外部离线产出）完成。
5. **可删边界**：各步产物子目录（`tables/`、`reports/`、`figures/`、`outputs/`）可整体删除后重跑
   （`run_all.py` 起链前按 `run_order.json` 的 `outputs` 预建父目录）；但 `raw/` **输入层**
   （`03_gates/admet/raw/`、`03_gates/herg/cross_model/raw/`）、`03_gates/herg/herg_predictions.csv`
   与 `03_gates/herg/hpc_herg/`、以及 `config/ vendor/ tests/ docs/ 00_inputs/ scripts/` 不得删。

## 3. 输入 → 产物映射

| 输入 | 步骤 | 关键产物 |
|---|---|---|
| `00_inputs/raw/inputs/top5000_best_pose.sdf` | 02–04 | `01_pose_qc/outputs/SOP-05_poses_normalized.sdf`、`01_pose_qc/tables/pose_gate.csv` |
| `00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv` | 02、03、05 | `01_pose_qc/tables/pose_pool.csv`、`02_druglikeness/tables/descriptors.csv`（pKa/cLogD 列）、`03_gates/tables/candidate_universe.csv` |
| `.../run01/01_redocking/receptor/4EIY_protein.pdb` | 04 | 16 残基接触指纹、锚点命中、盒子几何 |
| `.../run01/02_enrichment_5098/library/quick_5100_library.csv` | 04、07 | 5,098 集 IFP（锚点门证据链）、`03_gates/tables/maxTc_actives.csv`（479 参考活性） |
| `.../run01/02_enrichment_5098/scores/gnina_ranking.csv` + `poses_raw/poses.sdf` | 04、10 | 阴性对照 NEG1/NEG2（包内重算） |
| `.../screening_kit/07_independent_validation/independent_a2a_chembl.csv` | 09 | `03_gates/tables/known_ligand_similarity.csv`（已知物核验） |
| `03_gates/herg/herg_predictions.csv`（**外部模型产出的落门表**，已登记哈希） | 08、09、10 | `herg_gate` 列、门叠加链的 `6_hERG<0.5` 步、名单 `herg_prob` 与区间 |
| `03_gates/herg/hpc_herg/*`（外部运行的原始产出，**链上不读取**） | — | 惰性审计证据，见该目录 `README.md` |
| `03_gates/herg/cross_model/raw/*`（ChEMBL hERG 活性表 + 抓取元信息，**输入层**） | 08b | 自建第二模型的训练/评估数据；哈希由 `fetch_meta.json` 登记，步骤**不重抓** |
| `03_gates/herg/cross_model/{tables,reports,figures}/*`（产物） | 08b | `herg_cross_model.json`、一致性表、偏离分子清单 |
| `03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv`（外部模型全量表，**输入层**） | 10、10b、10c | 池级 ADMET 注释、hERG 阈值锚定与保留曲线 |
| `04_selection/tables/list_main_QED0.50_anchor4.csv`（主口径 50 名单）+ `03_gates/herg/cross_model/tables/herg_cross_model_predictions.csv`（第二模型 hERG） | 15 | `07_md_batch/tables/md_batch20.csv`（MD 一级批次 20）、`md_batch_excluded.csv`、`reports/md_batch.json`、`MD_BATCH.md` |
| `.../redocking/00_3EML_structure/3EML_crystal_zma.pdb`、SOP-02 基准 | 溯源 | 报告中的协议来源说明（不参与选序） |

## 4. 步骤目录约定

每个步骤目录含 `scripts/`（可执行）、`tables/`（结构化表）、`figures/`（图）、
`reports/`（JSON 摘要）；（步骤 1 另有 `outputs/` 存 SDF）。**报告只读取这些产物**，
不重新计算，保证「报告数字 ↔ 产物」同源。若某步依赖第三方原始表，则在该步目录下另设
`raw/` 存**输入**（如 `03_gates/admet/raw/`、`03_gates/herg/cross_model/raw/`）——
`raw/` 与 `scripts/` 同属不可删层，其余子目录均为可重建产物。

## 5. CNS MPO 标注列（步骤 05 产出）

| 列 | 含义 |
|---|---|
| `basic_pka` | MolGpKa 预测的最强碱性中心 pKa（点值，无区间）；无碱性中心者上游记 0 |
| `clogd74` | cpLogD v2.0 / ChEMBL33 的 cLogD7.4 点估计 |
| `clogd74_lower_80` / `clogd74_upper_80` | 同模型的 80% 预测区间（`clogd_confidence=0.8`） |
| `cns_mpo` | 点估计（cLogD 取 `clogd74`），Wager 2010 Table 1 |
| `cns_mpo_lo80` / `cns_mpo_hi80` | 悲观界（cLogD 取区间上界）/ 乐观上界（取区间下界） |
| `cns_mpo_pi80_span` | `hi80 − lo80`，即该分子的区间宽度 |

均为**标注列**（不作门），口径冻结在 `config/frozen_rule.json → step2_druglikeness.cns_mpo`。
