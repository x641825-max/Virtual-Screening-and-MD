# SOP_ML 骨架跑通与确定性基线记录（2026-09-24）

本文件记录 `SOP_ML` 骨架「恢复输入层 → 跑通全链 → 取确定性不动点」这一轮的实测读数。
所有读数均取自本包 `logs/` 与 `05_robustness/determinism/`，可独立复算。

## 1. 输入层恢复（15 件，逐条 sha256 校验）

来源：源包 `SOP_5`（同工作区，平级）。恢复范围＝跑链必需的大体积件，
**审计件与链上产物未恢复**（按设计由本轮自动重建）。

| 组 | 文件数 | 说明 |
|---|---|---|
| 原始输入表 | 2 | 全库排名表（换表对象）、5,000 姿态 SDF |
| 原始归档 | 3 | screening_kit / redocking / run01_redock_enrichment |
| 白名单解包产物大件 | 9 | 含 run01 的 `poses_raw/poses.sdf`（第 1 步需要） |
| ADMET 冻结全量表 | 1 | 49,963 × 124 |

校验方式：按 `RESTORE_LIST.sha256` 过滤出上述 15 条后 `sha256sum -c` → **15/15 OK**。
另做包级比对：`00_inputs` 下与源包共有的 63 个文件逐字节比对，**仅 1 处不同**
（`00_inputs/scripts/01_verify_inputs.py`，即骨架的三处已登记适配之一）。

## 2. 四轮全链（每轮 19 步 + 末尾自检，均 exit 0）

| 轮 | 起跑 UTC | 确定性比对读数 | 说明 |
|---|---|---|---|
| r1 | 2026-09-24T03:54:31Z | 写基线 186 文件（`mode=baseline_written`） | **不构成通过** |
| r2 | 2026-09-24T03:58:20Z | 186 文件，差异 0（`mode=compared`） | 首次判据 |
| r3 | 2026-09-24T04:04:26Z | 186 文件，差异 0 | |
| r4 | 2026-09-24T04:09:31Z | 186 文件，差异 0 | 载体不动点确认轮 |

哈希范围＝186 个解析型产物文件（`06_deliverables` 58、`03_gates` 45、`04_selection` 21、
`05_robustness` 20、`01_pose_qc` 16、`02_druglikeness` 14、`07_md_batch` 7、`00_inputs` 5）；
排除 6 个自指文件（2 个比对报告 + 4 个交付层状态载体）与全部 `logs/`。

## 3. 7 载体收敛

载体清单与快照：`logs/history/carrier_round_r2.md5`、`carrier_round_r3.md5`、`carrier_round_r4.md5`。

- r2 → r3：7 项中 2 项变化（`06_deliverables/MANIFEST.sha256`、`PACK_MANIFEST.sha256`）
- r3 → r4：**7/7 逐字节一致** → 载体不动点

机制与源包一致：`13_assemble` 早于链末的 `12_determinism`，交付层副本永远落后一轮；
索引/清单又落后交付层一轮，故稳定需要「写基线 + 连续两轮以上」。

## 4. 与源包 `SOP_5` 定稿轮的差异（逐条）

对共有产物做逐字节比对，差异仅以下三类，均属预期：

1. 确定性读数口径：本包哈希范围 186 文件，源包 187 文件——差 1 个为
   `03_gates/herg/hpc_herg/predictions.csv`（HPC 审计件，骨架未随包，见 `RESTORE_LIST.md`）。
   该差异只影响报告/验收件中引用这一读数的那一行。
2. `00_inputs/manifests/FREEZE.json`：仅 5 个脚本哈希不同——
   `01_verify_inputs.py`、`tests/make_index.py`、`tests/verify_pack.py`（已登记适配）
   与新增的 `00b_register_swap_input.py`、`tests/check_swap_csv.py`。
3. 包名相关文本（`INDEX_PACK.md`、`PACK_MANIFEST.sha256` 抬头与文件清单）：
   归一化包名后与源包逐行一致。

其余全部产物（名单、决策卡、门栈、稳健性、MD 批次与最优性证书、图、姿态 SDF 等）
与源包**逐字节相同** → 输入层恢复正确、规则与代码行为一致。

## 5. 末尾自检（r4 后复核）

| 自检 | 结果 |
|---|---|
| `tests/verify_pack.py` | 全部通过（含 MILP 最优性证书、独立重解、交付层副本同源等） |
| `tests/check_isolation.py` | 34 个包内 Python 文件，越界引用 0 处 |
| `tests/make_index.py` | 生成 `INDEX_PACK.md` 与 `PACK_MANIFEST.sha256`（225 个文件） |

## 6. 过程中的一次中断（如实登记）

r2 首次启动被中断，残留锁 `logs/.run_all.lock` 内的 pid 被沙箱回收后由新进程占用，
`acquire_lock()` 的存活探测因此误判「已有链在运行」而 ABORT（exit 3）。
处置：确认无活动链后删除残留锁，重跑 r2 成功。锁机制本身按设计工作，
本例暴露的是**沙箱 pid 回收导致的存活探测假阳性**，不影响产物正确性
（r2 之后的 r3、r4 均给出 186/186 一致）。

## 7. 结论

骨架可用：输入层可恢复且校验通过，全链可跑通，确定性基线已建立并连续两轮以上一致，
7 载体达到不动点。本轮**不构成**对任何分子结论的重建——它跑的是源包那批数据，
结论与 `SOP_5` 定稿轮相同；ML 数据换版后必须重取基线（见 `SWAP_SPEC.md` §4.6）。
