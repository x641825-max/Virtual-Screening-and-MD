# RESTORE_LIST · 未随骨架的文件清单与恢复方式

`SOP_ML` 是**精简复制**的骨架：为了让包体停在几 MB 级，复制时排除了大体积输入与全部链上产物。
被排除的文件共 **413 个 / 429.1 MB**，逐条哈希见 `RESTORE_LIST.sha256`（sha256 取自派生时的源包，可用于恢复后校验）。

> 说明：本清单里的 sha256 是**源包 `SOP_5` 当时的值**。换数据版只替换排名表一个文件，
> 若你后来换了表，该表的哈希应以 `00b_register_swap_input.py` 登记的新值为准。

## 1. 概览（按组）

| 组 | 文件数 | 体积 |
|---|---|---|
| 白名单解包产物中的大件（>400 KB） | 9 | 134.2 MB |
| 原始归档（3 个 zip） | 3 | 126.8 MB |
| ADMET-AI 冻结全量表（77 MB，换表后需重生成） | 1 | 77.4 MB |
| 原始输入表（含换数据对象与姿态 SDF） | 2 | 31.8 MB |
| 步骤 1 产物 | 13 | 18.2 MB |
| 步骤 3 产物 | 30 | 12.4 MB |
| 可视化图册（手工步，非判据） | 21 | 9.7 MB |
| hERG 落门表 HPC 审计件 | 1 | 4.7 MB |
| 步骤 2 产物 | 12 | 4.7 MB |
| 运行日志（含历史快照） | 205 | 4.3 MB |
| 交付层产物 | 61 | 3.2 MB |
| 步骤 4 产物 | 20 | 0.7 MB |
| 步骤 5 产物（含 G 臂确定性基线） | 21 | 0.6 MB |
| 步骤 16 产物 | 6 | 0.3 MB |
| 输入清单与冻结登记（首次跑链自动重建） | 5 | 0.0 MB |
| 其它 | 3 | 0.0 MB |

## 2. 体积最大的 14 项

| 文件 | 体积 | 源包 sha256（前 16 位） |
|---|---|---|
| `00_inputs/extracted/run01/02_enrichment_5098/poses_raw/poses.sdf` | 123.6 MB | `a2df15d0270ed50e…` |
| `03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv` | 77.4 MB | `5d693ea7d9d06e29…` |
| `00_inputs/raw/archives/screening_kit.zip` | 63.8 MB | `21e859f7a03c5693…` |
| `00_inputs/raw/archives/redocking.zip` | 38.1 MB | `562d40fecb692a1a…` |
| `00_inputs/raw/archives/run01_redock_enrichment.zip` | 24.9 MB | `dc439580ea71d132…` |
| `00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv` | 18.8 MB | `dae23885e2748ca8…` |
| `01_pose_qc/outputs/SOP-05_poses_normalized.sdf` | 13.8 MB | `64ad6fd5d1640b96…` |
| `00_inputs/raw/inputs/top5000_best_pose.sdf` | 13.1 MB | `1040e55a80d96139…` |
| `03_gates/herg/hpc_herg/predictions.csv` | 4.7 MB | `56e989f96e62d83a…` |
| `03_gates/admet/tables/admet_annotation_pool.csv` | 4.4 MB | `1a69ec9c8c36e46f…` |
| `03_gates/tables/candidate_universe.csv` | 4.2 MB | `4a85fea05f8ff4ed…` |
| `02_druglikeness/tables/descriptors.csv` | 3.0 MB | `13a0bfacd7b781f5…` |
| `00_inputs/extracted/redocking/SOP-02/02_cleaned_benchmark/decoys_cleaned.csv` | 2.8 MB | `3026c8bef3eb00b4…` |
| `00_inputs/extracted/run01/02_enrichment_5098/scores/gnina_all_poses.csv` | 2.7 MB | `6ae9ce5a4280f926…` |

## 3. 恢复方式

### 3.1 从源包复制（推荐，保持逐字节一致）

```bash
# 在源包所在工作区执行；SOP_ML 与 SOP_5 平级
for f in $(awk '{print $2}' SOP_ML/RESTORE_LIST.sha256); do
  src="SOP_5/$f"; dst="SOP_ML/$f"
  [ -e "$src" ] || continue
  mkdir -p "$(dirname "$dst")"
  cp -p "$src" "$dst"
done
```

### 3.2 校验

```bash
cd SOP_ML && sha256sum -c RESTORE_LIST.sha256
```

全部 `OK` 即与派生时的源包一致（`00_inputs/manifests/` 与 `logs/` 会因本轮之前的状态不同而可能不等，属正常；
这两类也可直接跳过，让首次运行重建）。

### 3.3 可再生成的部分（不必复制）

| 组 | 谁重建 |
|---|---|
| 各步骤 `tables/`、`reports/`、`figures/`、`outputs/` | `python3 run_all.py` |
| `00_inputs/manifests/`（INPUT_MANIFEST / input_verify / FREEZE / EXTRACTED_MANIFEST） | `01_verify_inputs` 与 `00_unpack_inputs` |
| `05_robustness/determinism/`（G 臂基线） | `python3 run_all.py --freeze`（写基线，**不构成通过**）后连续跑 2–3 轮 |
| `08_visuals/`（图册） | `python3 08_visuals/scripts/make_visuals.py`（手工步，不在链上） |
| `06_deliverables/` 交付件 | 步骤 13/14 与 `run_all.py` |
| `03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv` | 外部 admet-ai 作业（`10b0_admet_regenerate.py` 生成输入、手工执行）；**换表后必须重生成** |
| `logs/` | 每轮运行自动滚动重建 |

## 4. 缺文件时的链上行为

- `00_inputs/scripts/01_verify_inputs.py`：把缺失输入登记为 `md5_got = MISSING`、`match = false`，
  以 **exit 4** 显式失败并打印恢复指引（不再以 traceback 终止）。
- `00_inputs/scripts/00_unpack_inputs.py`：缺归档时无法解包，相关步骤随后报错。
- `03_gates/scripts/10b_admet_annotation.py` / `10c_herg_threshold.py`：缺 ADMET 冻结表时前置检查失败。
- 设计原则：**显式失败，不静默降级**（沿用 SOP_5 的口径）。
