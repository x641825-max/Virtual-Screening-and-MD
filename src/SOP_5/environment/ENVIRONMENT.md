# 环境快照与依赖（SOP_5）

## 1. 运行时

| 项 | 值 |
|---|---|
| Python | 3.12.2 |
| 平台 | Linux x86_64（glibc 2.36） |
| 依赖 | RDKit 2025.9.5、pandas 2.3.3、numpy 2.4.6、scikit-learn 1.9.1、scipy 1.15.3、matplotlib 3.11.1、Pillow 12.3.0 |
| 中文字体 | Noto Sans CJK SC（绘图用；缺失时自动回退并提示） |

完整快照见 `environment/env_snapshot.json`，冻结版本见 `environment/requirements.txt`。

## 2. 明确不包含的能力（以及本包如何处置）

| 能力 | 状态 | 处置 |
|---|---|---|
| hERG 预测（ADMET-AI / chemprop / torch） | **不在本环境** | 门做成外部接口：生成补算清单，表缺失时 `herg_gate=pending` 并写入名单、报告、验收；不伪造预测值 |
| `gmx_MMPBSA` | 缺失 | 下游二级判据降级为「锚点接触保留率 + H 键占据率 + 能量分解」，并在 `06_deliverables/SOP-05_downstream_MD_plan.md` 显式标注 |
| gnina / AutoDock Vina 可执行文件 | 不在本环境 | 本包不重跑对接；姿态与打分由输入提供，本包只核验一致性与几何 |
| zip/unzip、ps/time 等系统工具 | 缺失 | 用 Python `zipfile`、`/proc` 与标准库替代 |

## 3. 确定性设置

- 主链无随机数；稳健性臂使用 `numpy.random.default_rng(PROJECT_SEED)`（PCG64）。
- `run_all.py` 强制 `PYTHONHASHSEED=0`、`OMP_NUM_THREADS=1`（避免线程调度引入浮点归约差异）。
- 排序一律带次级键（cid 字典序），杜绝并列次序抖动。
- G 臂对全量产物做 sha256 比对（`05_robustness/determinism/`）。

## 4. 资源与耗时

| 步骤 | 主要成本 | 说明 |
|---|---|---|
| 00 解包 | 约 300 MB 写盘 | 大体积成员由 `extract_large` 开关控制 |
| 02–04 姿态/指纹 | 5,000 姿态 × 16 残基距离矩阵；5,098 姿态 IFP | 单核数分钟 |
| 05 描述符 | 4,990 × RDKit（含 PAINS/Brenk catalog 匹配） | 单核 1–3 分钟 |
| 11 稳健性 | A/E 各 3 档 × 100 次重抽样选序、F 全池 Butina | 单核 3–8 分钟 |
| 合计 | — | 约 10–20 分钟（视机器） |
