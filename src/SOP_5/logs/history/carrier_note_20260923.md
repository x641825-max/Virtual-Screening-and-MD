# 收尾修复后的载体不动点证据（2026-09-23，19 步版）

脚本层定稿（`run_all.py` 预建产物目录、`08b` raw 输入层前置检查、`vendor/master.py` 缺输入指引、
`tests/verify_pack.py` 第 [20] 组、`08_visuals` 图册异常表条件化）后：

- 轮次 B：`run_all.py --refreeze`（写基线，不构成通过）
- 轮次 C：`run_all.py` → G 臂 `mode=compared` / `n_diff=0`（187/187）
- 轮次 D：`run_all.py` → `n_diff=0`；与 C 比，5/7 载体变化（报告、验收、交付层比对副本、
  `MANIFEST.sha256`、`PACK_MANIFEST.sha256`）——即 README §6 所述「载体滞后一轮」
- 轮次 E：`run_all.py` → `n_diff=0`；与 D 比，仅 `MANIFEST.sha256`、`PACK_MANIFEST.sha256` 变化
- 轮次 F：`run_all.py` → `n_diff=0`；**与 E 比 7/7 载体逐字节一致（不动点）**

证据文件：`carrier_roundE.md5`、`carrier_roundF.md5`（7 个载体/自检文件）。
注：`08_visuals/` 属链外目录（不进 G 臂），本日最后一次重建图册后其哈希进入 `PACK_MANIFEST.sha256`，
故根清单随后一次重建保持一致，不影响 187 个链上产物的不动点。

---

## 附：v1.7.0（MD 批次可证明最优）的收敛轮次（2026-09-23 续）

规则层升 v1.7.0（`step6_md_batch.objective`，旧贪心降为对照）+ 步骤层 15 号脚本接入 MILP +
自检新增第 [21] 组 + 文档/图册同步后：

- 轮次 G：`run_all.py --refreeze`（写基线，不构成通过）
- 轮次 H：`run_all.py` → `mode=compared` / `n_diff=0`（187/187）；与 G 比 5/7 载体变化（滞后一轮）
- 轮次 I：`run_all.py` → `n_diff=0`；与 H 比仅 `MANIFEST.sha256`、`PACK_MANIFEST.sha256` 变化
- 轮次 J：`run_all.py` → `n_diff=0`；**与 I 比 7/7 载体逐字节一致（不动点）**
- 图册重建（链外 `08_visuals/make_visuals.py`）+ `tests/make_index.py` 重生成根清单 → 轮次 K0
- 轮次 K：`run_all.py` → `n_diff=0`；**与 K0 比 7/7 载体逐字节一致**

证据文件：`carrier_roundG/H/I/J/K0/K.md5`。终态读数：`verify_pack` 193 项 / 21 组全过、
`check_isolation` 0 越界、硬门 9 条（6 通过 / 2 待 G 臂时点 / 1 未通过＝σ=0.10 主判据）。
MD 批次证书：t*=4、CNN_VS 跨度 0.7337（池内全量程）、Σ CNN_VS 136.5031、与 v1.6.0 贪心对照 6 换 6。

---

## 附二：显示层修正后的最终收敛轮次（2026-09-23 收尾）

最后又修了两处**显示层**问题（① 证书阶段表把求解器返回的 `fun` 直接展示成 `--4` 之类的双负号；
② 「保留 + 标注」机时代价多乘了 10×（应为 0.2 卡日而非 2.0）；③ 报告 §10 把轴数写成 `len(字符串)`＝20），
故再做一轮完整重冻与收敛：

- 轮次 Z0：`run_all.py --refreeze`（写基线，不构成通过）
- 轮次 Z1 / Z2：普通跑 → `n_diff=0`；与上一轮比各 5/7 载体变化（滞后两轮：报告/验收/交付层比对副本/两份 MANIFEST）
- 轮次 Z3：普通跑 → `n_diff=0`；与 Z2 比仅 `MANIFEST.sha256`、`PACK_MANIFEST.sha256` 变化
- 轮次 Z4：普通跑 → `n_diff=0`；**与 Z3 比 7/7 逐字节一致（不动点）**
- 图册重建 + `tests/make_index.py` 重生成根清单 → W0；轮次 W：**与 W0 比 7/7 一致**

证据文件：`carrier_round{Z0,Z1,Z2,Z3,Z4,W0,W}.md5`（选轮证据中 `Z0→Z1` 5 处、`Z1→Z2` 5 处、
`Z2→Z3` 2 处、`Z3→Z4` 0 处，与 README §6 描述的「载体滞后」规律一致）。

---

## 附三：链外写入事件与根清单重建后的收敛轮次（2026-09-23 深夜 → 09-24 00:15 UTC）

**事件（2026-09-23 19:08 UTC）**：链外目录 `08_visuals/` 的脚本、9 张图、`gallery.html`、`README.md`、
`visuals_meta.json` 被重新生成（19:08:51–19:09:08），另有 `vendor/__pycache__`、根 `__pycache__` 残留。
根清单 `PACK_MANIFEST.sha256`（18:52:52 写入）对这 12 个条目随即过期：
`sha256sum -c` 报 12 项 FAILED、231 项一致。链上产物层 7 个载体当时仍与 `carrier_roundX.md5` 逐字节一致。

**处置与复验**：
- `tests/make_index.py` 重建 `INDEX_PACK.md` + `PACK_MANIFEST.sha256` → `sha256sum -c` **243/243 一致**
  （此前 12 项 FAILED 全部消除）。
- 包内 `__pycache__` 残留清除（0 处）；`sha256sum -c`、`verify_pack`、`check_isolation` 现场复跑均通过。

**重建后的收敛轮次**（每轮 19 步 + 自检，约 304–306 s）：

| 轮次 | 起跑（UTC） | 参数 | G 臂读数 | 与上一轮的 7 个载体比对 |
|---|---|---|---|---|
| Y0 | 2026-09-23T23:53:37Z | `--refreeze`（写基线；**不构成通过**，含根清单重建） | `mode=baseline_written`，`n_files=187` | —（基线轮） |
| Y | 2026-09-23T23:58:46Z | 默认 | `mode=compared`、`n_diff=0`（187/187） | Y0≡Y **2/7**（报告、验收、交付层比对副本、MANIFEST、PACK_MANIFEST 变） |
| Z | 2026-09-24T00:04:06Z | 默认 | `mode=compared`、`n_diff=0` | Y≡Z **2/7**（同上 5 项继续随滞后机制更新） |
| AA | 2026-09-24T00:09:34Z | 默认 | `mode=compared`、`n_diff=0` | Z≡AA **5/7**（仅 `MANIFEST.sha256`、`PACK_MANIFEST.sha256` 变） |
| BB | 2026-09-24T00:14:45Z | 默认 | `mode=compared`、`n_diff=0`（187/187） | **AA≡BB 7/7 逐字节一致（不动点）** |

证据文件：`carrier_roundY0.md5`、`carrier_roundY.md5`、`carrier_roundZ.md5`、`carrier_roundAA.md5`、
`carrier_roundBB.md5`（各 7 个载体/自检文件）。

**机制说明（与 README §6「载体滞后一轮」一致）**：`13_assemble`（第 17 步）早于 `12_determinism`（第 19 步），
交付层的确定性副本与验收表登记的「上一轮读数」永远落后一轮，故一个冻结周期需要连续两轮以上普通比对才收敛；
本次因先做根清单重建（相当于一次冻结事件），滞后链条比常规多走了一轮（2/7 → 2/7 → 5/7 → 7/7）。

**终态读数**：`mode=compared` / `n_diff=0`（187/187）、`verify_pack` 193 项 / 21 组全过、0 FAIL、
`check_isolation` 0 越界、根清单 243/243 一致、包内无 `__pycache__` 残留；
硬门仍为 6 通过 / 2 待时点 / 1 未通过（σ=0.10 主判据）。
