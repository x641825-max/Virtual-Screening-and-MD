#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""执行序第 16 步（步骤 6）：MD 批次选择 —— 从主口径名单（候选 + 2 对照）中挑一批进 MD
（名单与批次规模均见 frozen_rule.step4_selection / step6_md_batch，脚本不内嵌规模）。

为什么需要这一步
  名单层已经完成「结构身份 → 类药性 → 相似度/活性 → hERG → 多样性」的全部门禁，
  候选彼此骨架互不相同（Bemis-Murcko 骨架 1/骨架），因此「先跑哪几个」不再是
  过滤问题，而是**机时分配**问题：在有限的 MD 预算下，怎样的一批分子能同时回答
  「谁能过 MD」与「对接分数是否预测 MD 稳定性」。

选序（frozen_rule.step6_md_batch；全部阈值来自规则集，脚本不内嵌阈值）
  ① 资格（v1.6.0 起）：第二 hERG 模型（集成口径）判不过门的候选**保留入池并逐分子标注**
     「模型依赖」（不剔除）；再做 Bemis-Murcko 通用骨架去重（保留 CNN_VS 高者），
     避免同系物重复占位。
  ② 可行域：|S| = n_candidates、成对 ECFP4 Tanimoto < 0.65、骨架 ≤1、化学型 ≤ 8、
     AM1-BCC 风险档 ≤ 3、≥5 席来自 CNN_VS 下半区（相对名单候选中位）。
  ③ 目标（v1.7.0 起，frozen_rule.step6_md_batch.objective）：在②的可行域内用 MILP
     （scipy.optimize.milp，HiGHS 后端）按**字典序**求可证明最优解 ——
       阶段 1  max 分层覆盖下限 t（CNN_VS / cLogP / CNS MPO 三轴 × 4 分位，资格池自分布）
       阶段 2  在 t=t* 内 max CNN_VS 跨度（max − min）
       阶段 3  仍并列时 max Σ CNN_VS（对接质量兜底）
       阶段 4  仍并列时 max Σ w_i·x_i（w_i = n − i，按 CNN_VS 降序）——规范化并列打破，
               使名册与求解器内部 tie-break 无关
     每阶段必须 status=0（optimal）且 mip_gap=0；且各阶段目标值必须等于**名册实算值**
     （防「假跨度」建模错误：跨度若只用单侧下界式 u ≥ z_i·x_i，优化器会把 u 顶到上界、
     把 l 压到下界，得到与所选分子无关的假跨度。本脚本用独热指派把 u/l 钉在所选分子上）。
  ④ 对照：v1.6.0 的贪心填充（化学型锚点 → 下半区补齐 → 高分填充）作为**对照实现**
     (heuristic_reference) 在同一份报告里留痕；其「是否为 max Σ CNN_VS 口径的全局最优」
     由本脚本独立求解该口径最优值后**按实测判定**（SOP_5 批次实测是，换数据后重判）。
  ⑤ 固定 2 位阴性对照（NEG1/NEG2），不占候选名额。

本步是分析步：只做取舍与留痕，不修改名单、不改任何门。若 MD 结果导致换人，
必须改规则并重跑整链（禁止手工编辑名单文件）。

输出：07_md_batch/tables/md_batch20.csv、tables/md_batch_excluded.csv、
      reports/md_batch.json（含最优性证书 optimality / 对照 heuristic_reference）、
      figures/fig_md_batch.png(+.svg)、MD_BATCH.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from rdkit import RDLogger
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from scipy.optimize import Bounds, LinearConstraint, milp

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import plots
from vendor.io_utils import (StepLog, read_json, r4, r6, sha256_file, write_csv, write_json,
                             write_text)

RDLogger.DisableLog("rdApp.*")
OUT = ROOT / "07_md_batch"
LIST_MAIN = ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv"
XMODEL = ROOT / "03_gates/herg/cross_model/tables/herg_cross_model_predictions.csv"
FPS = ROOT / "02_druglikeness/tables/fingerprints_ecfp4.npz"

# 对照实现（v1.6.0 贪心填充）的三类标签，仅用于 heuristic_reference 留痕
HEUR_ANCHOR = "化学型覆盖锚点（该类最高分）"
HEUR_LOW = "分数分层（CNN_VS 下半区）"
HEUR_FILL = "高分填充"
# 本版（v1.7.0）三类标签：由**反事实重解**判定，见 classify_roles()
REASON_CTRL = "阴性对照（固定 2 位）"
REASON_COVER = "覆盖下限约束席位（移除此席 t 会下降）"
REASON_SPAN = "跨度端点席位（移除此席跨度收窄）"
REASON_FILL = "质量兜底填充（移除后 t 与跨度均不变）"

OPT_OPTS = {"presolve": True, "disp": False}
MIP_GAP_TOL = 1e-9


def load_fps() -> dict:
    z = np.load(FPS, allow_pickle=True)
    out = {}
    for cid, bits in zip(z["cid"], z["bits"]):
        bv = ExplicitBitVect(2048)
        bv.SetBitsFromList(np.flatnonzero(bits).tolist())
        out[str(cid)] = bv
    return out


def round_half(x: float) -> int:
    """Python round 的银行家舍入不适用于名额口径；此处用四舍五入（与名单层一致）。"""
    return int(np.floor(x + 0.5))


def bins_for(vals: np.ndarray, ref: np.ndarray, k: int):
    """按参考分布 ref 的 100/k 分位切 k 层；返回（每值层号, 层边界）。"""
    edges = np.array([float(np.percentile(ref, 100 * i / k)) for i in range(1, k)])
    return np.searchsorted(edges, vals, side="left"), edges


def greedy_reference(pool: pd.DataFrame, pool_ix: pd.DataFrame, med: float, n_slots: int,
                     min_low: int, chem_cap: int, risk_cap: int) -> list[dict]:
    """v1.6.0 贪心填充（对照实现，逐字沿用其选序）：① 化学型锚点 → ② 下半区补齐 → ③ 高分填充。"""
    sel: list[dict] = []
    used_chem: dict[str, int] = {}
    used_risk = 0

    def can_take(row) -> bool:
        return (used_chem.get(row.chemotype, 0) < chem_cap
                and used_risk + int(row.risk_am1bcc) <= risk_cap)

    def take(row, reason: str) -> None:
        nonlocal used_risk
        sel.append({"cid": str(row.cid), "reason": reason})
        used_chem[row.chemotype] = used_chem.get(row.chemotype, 0) + 1
        used_risk += int(row.risk_am1bcc)

    for _chem, grp in pool.groupby("chemotype", sort=True):
        if len(sel) >= n_slots:
            break
        take(grp.iloc[0], HEUR_ANCHOR)
    n_low = sum(1 for s in sel if float(pool_ix.loc[s["cid"], "CNN_VS"]) < med)
    chosen = {s["cid"] for s in sel}
    for _, row in pool.iterrows():
        if n_low >= min_low or len(sel) >= n_slots:
            break
        if row.cid in chosen or row.CNN_VS >= med or not can_take(row):
            continue
        take(row, HEUR_LOW)
        chosen.add(str(row.cid))
        n_low += 1
    for _, row in pool.sort_values(["CNN_VS", "cid"], ascending=[False, True],
                                   kind="mergesort").iterrows():
        if len(sel) >= n_slots:
            break
        if row.cid in chosen or not can_take(row):
            continue
        take(row, HEUR_FILL)
        chosen.add(str(row.cid))
    return sel


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    mb = rule["step6_md_batch"]
    obj = mb["objective"]
    with StepLog(ROOT, "15_md_batch_select") as log:
        for d in ("tables", "reports", "figures"):
            (OUT / d).mkdir(parents=True, exist_ok=True)

        listed = pd.read_csv(LIST_MAIN, low_memory=False)
        cand = listed[listed.candidate_type == "candidate"].copy()
        ctrl = listed[listed.candidate_type == "negative_control"].copy()
        n_expected_cand = int(mb["list"]["n_candidates"])
        n_expected_ctrl = int(mb["list"]["n_negative_controls"])
        assert len(cand) == n_expected_cand, f"候选数 {len(cand)} != 规则 {n_expected_cand}"
        assert len(ctrl) == n_expected_ctrl, f"对照数 {len(ctrl)} != 规则 {n_expected_ctrl}"
        log.info(f"名单：{len(cand)} 候选 + {len(ctrl)} 对照（主口径 "
                 f"{mb['list']['main_arm']}）")

        # ---- 资格①：第二 hERG 模型（集成口径）不过门者
        # v1.6.0 起：**不再剔除**，而是保留入池并逐分子标注「模型依赖」。
        # 理由：该判定来自一个与主模型 κ≈0 的第二模型（近乎随机一致），且两个模型一致判高风险的
        # 候选实测很少；用一个未校准信号决定机时分配，纪律上最弱。机时成本 +0.2–0.3 卡日即可买回。
        xm = pd.read_csv(XMODEL, low_memory=False)
        xm = xm[xm.cid.isin(cand.cid)].copy()
        assert len(xm) == len(cand), f"第二模型预测缺失 {len(cand) - len(xm)} 条"
        ens_col = mb["eligibility"]["hERG_cross_model_column"]
        xm["hERG_cross_ens"] = xm[ens_col]
        _thr = float(rule["step3_gates"]["herg"]["threshold"])
        fail = xm[xm.hERG_cross_ens >= _thr]
        n_exp_fail = int(mb["eligibility"]["expected_n_hERG_cross_model_fail"])
        assert len(fail) == n_exp_fail, (f"第二模型不过门候选 {len(fail)} != 规则登记 "
                                        f"{n_exp_fail}；规则或上游已变，须显式复核")
        cand = cand.merge(xm[["cid", "hERG_cross_ens"]], on="cid", how="left")
        cand["herg_cross_model_flag"] = np.where(
            cand.hERG_cross_ens >= _thr, "cross_model_fail", "cross_model_pass")
        _retain = bool(mb["eligibility"].get("exclude_hERG_cross_model_fail", True))
        excl = []
        if _retain is False:
            cand["md_eligibility_note"] = np.where(
                cand.herg_cross_model_flag == "cross_model_fail",
                f"保留入池 + 标注：第二 hERG 模型（集成）不过门（hERG_cross_ens ≥ {_thr}），"
                f"主模型口径下过门；本层为机时分配，不产生安全性结论", "")
            log.info(f"资格①（保留+标注，不再剔除）：第二模型不过门候选 {len(fail)} 个 "
                     f"（{', '.join(str(c) for c in fail.cid)}）继续参与抽选；"
                     f"机时代价 +{len(fail) * float(mb['md_cost']['level1_10ns_cardday_per_molecule']):.2f} "
                     f"卡日量级")
        else:
            for r in fail.itertuples():
                excl.append({"cid": str(r.cid), "reason": "第二 hERG 模型（集成）不过门",
                             "detail": f"hERG_cross_ens={r.hERG_cross_ens:.6f} ≥ {_thr}"})
            cand = cand[~cand.cid.isin(fail.cid)].copy()
            log.info(f"资格剔除（hERG 第二模型不过门）：{len(fail)} 个 "
                     f"（{', '.join(str(c) for c in fail.cid)}）")

        pool = cand.copy()
        # ---- 资格②：通用骨架去重（保留 CNN_VS 高者）
        pool = pool.sort_values("CNN_VS", ascending=False)
        dup = pool[pool.duplicated("scaffold_generic", keep="first")]
        for r in dup.itertuples():
            keep = pool[pool.scaffold_generic == r.scaffold_generic].iloc[0]
            excl.append({"cid": str(r.cid), "reason": "通用骨架重复（同系物占位）",
                         "detail": f"scaffold_generic 已有更高分者 {keep.cid}"
                                   f"（{keep.CNN_VS:.4f} > {r.CNN_VS:.4f}）"})
        pool = pool[~pool.duplicated("scaffold_generic", keep="first")].reset_index(drop=True)
        log.info(f"资格集：{len(pool)} 个候选（通用骨架去重再剔除 {len(dup)} 个）")

        # ---- 配额与分层参数
        n_batch = int(mb["batch_size"])
        n_slots = int(mb["n_candidates"])
        caps = mb["caps"]
        chem_cap = round_half(8 * n_batch / 20)
        risk_cap = round_half(3 * n_batch / 20)
        assert chem_cap == int(caps["chem_cap_value"]), \
            f"化学型上限 {chem_cap} != 规则 {caps['chem_cap_value']}"
        assert risk_cap == int(caps["risk_cap_value"]), \
            f"风险档上限 {risk_cap} != 规则 {caps['risk_cap_value']}"
        med = float(cand.CNN_VS.median())
        min_low = int(mb["stratification"]["min_below_median_candidates"])
        tc_cap = float(caps["pairwise_tc_max"])
        # ---- 资格池化学型 vs 帽 vs 名额的**咬合度**（决定可行域冗余与反事实可解性）
        # 名额可填上限 = Σ min(帽, 该化学型池内计数)；差 0 ⇒ 每个化学型都必须取满，
        # 于是「计数 ≤ 帽」的化学型其**全部**成员被结构性锁定（禁掉任一个可行域即为空）。
        chem_counts = {str(k): int(v) for k, v in pool.chemotype.value_counts().items()}
        fill_ub = int(sum(min(chem_cap, v) for v in chem_counts.values()))
        slot_slack = int(fill_ub - n_slots)
        n_forced = int(sum(v for v in chem_counts.values() if v <= chem_cap)) \
            if slot_slack == 0 else 0
        log.info(f"化学型咬合：池内 {chem_counts}、帽 {chem_cap} ⇒ 可填上限 {fill_ub} "
                 f"vs 名额 {n_slots}（冗余 {slot_slack} 席）"
                 + (f"；**零冗余** ⇒ {n_forced} 席被结构性锁定" if slot_slack == 0 else ""))
        log.info(f"批次 {n_batch} = {n_slots} 候选 + {n_expected_ctrl} 对照；"
                 f"化学型 ≤{chem_cap}、风险档 ≤{risk_cap}、骨架 ≤1、成对 Tc <{tc_cap}；"
                 f"CNN_VS 中位（相对 {len(cand)} 候选）={med:.4f}，下半区下限 {min_low} 个")

        # ---- 目标函数参数（全部来自规则集的 objective 块）
        axes = [(a["column"], a["label"]) for a in obj["coverage_axes"]]
        K = int(obj["n_bins_per_axis"])
        stage_specs = obj["lexicographic"]
        assert len(stage_specs) == 4, f"规则登记的字典序阶段数 {len(stage_specs)} != 4"
        assert [s["stage"] for s in stage_specs] == [1, 2, 3, 4], "阶段编号必须为 1–4"
        assert str(obj["bin_reference"]).startswith("eligible_pool"), \
            f"分位参照必须是资格池自分布，规则给的是 {obj['bin_reference']}"
        exp_cert = mb["expected_certificate"]
        assert int(exp_cert["bins_per_axis"]) == K \
            and int(exp_cert["n_bins"]) == K * len(axes), \
            f"登记分位格数 {exp_cert['n_bins']} != {K}×{len(axes)}"
        log.info(f"目标 v1.7.0（{obj['kind']}）：字典序 " + " → ".join(
            f"阶段{s['stage']} {s['name']}" for s in stage_specs) +
            f"；分位 {len(axes)} 轴 × {K} 层（参照 {obj['bin_reference']}）")

        # ---- 分层：按资格池自身分布切分位
        ids = [str(c) for c in pool.cid]
        idx_of = {c: i for i, c in enumerate(ids)}
        pool_ix = pool.set_index("cid")
        zvals = pool.CNN_VS.to_numpy(float)
        bin_mem: dict[tuple[str, int], list[str]] = {}
        bin_edges: dict[str, list[float]] = {}
        for col, name in axes:
            idx, edges = bins_for(pool[col].to_numpy(float), pool[col].to_numpy(float), K)
            bin_edges[name] = [float(e) for e in edges]
            for b in range(K):
                bin_mem[(name, b)] = [ids[i] for i in range(len(ids)) if idx[i] == b]
        n_bins = sum(1 for m in bin_mem.values() if m)
        empty_bins = [f"{k[0]}|Q{k[1] + 1}" for k, m in bin_mem.items() if not m]
        log.info(f"分位边界（资格池自分布）：" + "；".join(
            f"{k} {[round(e, 4) for e in v]}" for k, v in bin_edges.items()) +
            f"；有效分位格 {n_bins}（空格 {len(empty_bins)}）")

        # ---- 指纹与成对 Tc（约束与留痕共用）
        fps = load_fps()
        Tc = np.zeros((len(ids), len(ids)))
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                v = DataStructs.TanimotoSimilarity(fps[ids[i]], fps[ids[j]])
                Tc[i, j] = Tc[j, i] = v
        tc_conflict = [(i, j) for i in range(len(ids)) for j in range(i + 1, len(ids))
                       if Tc[i, j] >= tc_cap]

        # ---- MILP：变量布局 x(0..n-1) | t | u(max) | l(min) | y(0..n-1) | w(0..n-1)
        n = len(ids)
        n_var = 3 * n + 3
        I_T, I_U, I_L = n, n + 1, n + 2
        I_Y, I_W = n + 3, 2 * n + 3
        Mbig = float(zvals.max() - zvals.min()) + 1.0

        base: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        A = np.zeros((1, n_var)); A[0, :n] = 1.0
        base.append((A, np.array([float(n_slots)]), np.array([float(n_slots)])))
        for c in sorted(pool.chemotype.unique()):
            A = np.zeros((1, n_var))
            A[0, [idx_of[i] for i in ids if pool_ix.loc[i, "chemotype"] == c]] = 1.0
            base.append((A, np.array([-np.inf]), np.array([float(chem_cap)])))
        A = np.zeros((1, n_var)); A[0, :n] = pool.risk_am1bcc.to_numpy(float)
        base.append((A, np.array([-np.inf]), np.array([float(risk_cap)])))
        A = np.zeros((1, n_var)); A[0, :n] = (zvals < med).astype(float)
        base.append((A, np.array([float(min_low)]), np.array([np.inf])))
        for i, j in tc_conflict:
            A = np.zeros((1, n_var)); A[0, i] = A[0, j] = 1.0
            base.append((A, np.array([-np.inf]), np.array([1.0])))
        # 跨度（max_S z − min_S z）的严格线性化：独热指派 + 紧约束
        A = np.zeros((1, n_var)); A[0, I_Y:I_Y + n] = 1.0
        base.append((A, np.array([1.0]), np.array([1.0])))
        A = np.zeros((1, n_var)); A[0, I_W:I_W + n] = 1.0
        base.append((A, np.array([1.0]), np.array([1.0])))
        A = np.zeros((1, n_var)); A[0, I_U] = 1.0; A[0, I_Y:I_Y + n] = -zvals
        base.append((A, np.array([0.0]), np.array([0.0])))            # u = Σ z_i y_i
        A = np.zeros((1, n_var)); A[0, I_L] = 1.0; A[0, I_W:I_W + n] = -zvals
        base.append((A, np.array([0.0]), np.array([0.0])))            # l = Σ z_i w_i
        for i in range(n):
            A = np.zeros((1, n_var)); A[0, I_Y + i] = 1.0; A[0, i] = -1.0
            base.append((A, np.array([-np.inf]), np.array([0.0])))    # y_i ≤ x_i
            A = np.zeros((1, n_var)); A[0, I_W + i] = 1.0; A[0, i] = -1.0
            base.append((A, np.array([-np.inf]), np.array([0.0])))    # w_i ≤ x_i
            A = np.zeros((1, n_var)); A[0, I_U] = 1.0; A[0, i] = -Mbig
            base.append((A, np.array([zvals[i] - Mbig]), np.array([np.inf])))   # u ≥ max_S z
            A = np.zeros((1, n_var)); A[0, I_L] = 1.0; A[0, i] = Mbig
            base.append((A, np.array([-np.inf]), np.array([zvals[i] + Mbig])))  # l ≤ min_S z
        # 「每层席位 ≥ t」：Σ_{i∈层} x_i − t ≥ 0
        bin_rows: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        for b, members in bin_mem.items():
            if not members:
                continue
            A = np.zeros((1, n_var))
            A[0, [idx_of[c] for c in members]] = 1.0
            A[0, I_T] = -1.0
            bin_rows.append((A, np.array([0.0]), np.array([np.inf])))

        zmin, zmax = float(zvals.min()), float(zvals.max())
        lb = [0.0] * n + [0.0, zmin, zmin] + [0.0] * (2 * n)
        ub = [1.0] * n + [float(n_bins), zmax, zmax] + [1.0] * (2 * n)
        bounds = Bounds(np.array(lb), np.array(ub))
        integrality = np.array([1] * n + [0, 0, 0] + [1] * (2 * n))

        def c_t() -> np.ndarray:
            c = np.zeros(n_var); c[I_T] = -1.0; return c

        def c_span() -> np.ndarray:
            c = np.zeros(n_var); c[I_U] = -1.0; c[I_L] = 1.0; return c

        def c_sum() -> np.ndarray:
            c = np.zeros(n_var); c[:n] = -zvals; return c

        def c_canon() -> np.ndarray:
            # 并列打破：按 CNN_VS 降序索引给权重（高分者权重大），保证名册唯一且与求解器无关
            c = np.zeros(n_var); c[:n] = -np.arange(n, 0, -1).astype(float); return c

        def solve(rows, c, *, allow_infeasible: bool = False):
            A_ = np.vstack([a for a, _b, _c in rows])
            b_ = np.concatenate([b for _a, b, _c in rows])
            c_ = np.concatenate([cc for _a, _b, cc in rows])
            r = milp(c=c, integrality=integrality, bounds=bounds,
                     constraints=LinearConstraint(A_, b_, c_), options=OPT_OPTS)
            if int(r.status) == 2 and allow_infeasible:
                # 反事实重解专用：禁掉一席后**可行域为空**不是求解失败，而是该席
                # 「结构性不可替代」的最强证据（本批化学型帽与名额恰好零冗余咬合）。
                # 字典序四阶段仍走默认 require_optimal 路径，status≠0 一律硬失败。
                return None
            if r.status != 0:
                raise SystemExit(f"[15] MILP 未取得最优解：status={r.status} {r.message}"
                                 "（规则要求每阶段 status=0 且 mip_gap=0）")
            gap = float(getattr(r, "mip_gap", np.nan))
            if not (gap <= MIP_GAP_TOL):
                raise SystemExit(f"[15] MILP 最优性间隙 {gap} > {MIP_GAP_TOL}；"
                                 "规则要求 gap=0，请检查求解器选项或问题规模")
            return r

        def lex_solve(extra_rows=()) -> tuple[list[str], list[dict]]:
            """按规则登记的字典序四阶段求解；返回（名册, 各阶段日志）。"""
            fixed: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
            stages = [("分层覆盖下限 t", c_t()), ("CNN_VS 跨度", c_span()),
                      ("Σ CNN_VS", c_sum()), ("规范并列打破", c_canon())]
            res, journal = None, []
            for k, ((name, c), spec) in enumerate(zip(stages, stage_specs), 1):
                assert name == spec["name"], f"阶段 {k} 名称不一致：{name} vs {spec['name']}"
                r = solve(list(base) + list(bin_rows) + list(extra_rows) + fixed, c)
                journal.append({"stage": k, "name": name, "objective": float(r.fun),
                                "status": int(r.status),
                                "mip_gap": float(getattr(r, "mip_gap", np.nan)),
                                "nodes": int(getattr(r, "mip_node_count", -1))})
                res = r
                if k == 1:
                    fixed.append((np.array([[1.0 if i == I_T else 0.0 for i in range(n_var)]]),
                                  np.array([float(round(r.x[I_T]))]),
                                  np.array([float(round(r.x[I_T]))])))
                elif k == 2:
                    row = np.zeros((1, n_var)); row[0, I_U] = 1.0; row[0, I_L] = -1.0
                    v = float(r.fun) * -1.0
                    fixed.append((row, np.array([v]), np.array([v])))
                elif k == 3:
                    row = np.zeros((1, n_var)); row[0, :n] = zvals
                    v = float(r.fun) * -1.0
                    fixed.append((row, np.array([v]), np.array([v])))
            roster = [ids[i] for i in range(n) if round(res.x[i]) == 1]
            return roster, journal

        def evaluate(members) -> dict:
            counts = {b: sum(1 for c in members if c in bin_mem[b]) for b in bin_mem}
            live = [v for b, v in counts.items() if bin_mem[b]]
            zv = np.array([float(pool_ix.loc[c, "CNN_VS"]) for c in members])
            return {"t": int(min(live)) if live else 0, "span": float(zv.max() - zv.min()),
                    "sum": float(zv.sum()), "counts": counts,
                    "n_occupied_bins": int(sum(1 for b in bin_mem
                                               if counts[b] > 0 and bin_mem[b]))}

        sel_ids, journal = lex_solve()
        ev_o = evaluate(sel_ids)

        # ---- 自洽断言（护栏 1）：各阶段目标值必须等于名册实算值
        assert ev_o["t"] == round(-journal[0]["objective"]), \
            f"阶段1 目标 {journal[0]['objective']} 与名册实算 t={ev_o['t']} 不一致"
        assert abs(ev_o["span"] - (-journal[1]["objective"])) < 1e-6, \
            (f"阶段2 跨度 {-journal[1]['objective']} 与名册实算 {ev_o['span']} 不一致"
             "——跨度线性化可能退化（假跨度）")
        assert abs(ev_o["sum"] - (-journal[2]["objective"])) < 1e-6, \
            f"阶段3 Σ CNN_VS 与名册实算不一致"

        # ---- 重解一致（护栏 2）：同输入再解一次，名册与目标向量必须逐位相同
        sel_ids_2, journal_2 = lex_solve()
        resolved_consistent = (sel_ids_2 == sel_ids
                               and all(abs(a["objective"] - b["objective"]) < 1e-9
                                       for a, b in zip(journal, journal_2)))
        assert resolved_consistent, "同一输入两次求解名册不一致（求解器不具确定性）"

        # ---- 登记值核对（上游漂移检测）
        assert ev_o["t"] == int(exp_cert["t_star"]), \
            (f"覆盖下限 t*={ev_o['t']} != 规则登记 {exp_cert['t_star']}；"
             "名单/描述符/规则已漂移，须显式复核并同步 expected_certificate")
        assert abs(round(ev_o["span"], 4) - float(exp_cert["span"])) < 1e-9, \
            f"跨度 {r4(ev_o['span'])} != 规则登记 {exp_cert['span']}"
        assert abs(round(ev_o["sum"], 4) - float(exp_cert["sum_CNN_VS"])) < 1e-6, \
            f"Σ CNN_VS {r4(ev_o['sum'])} != 规则登记 {exp_cert['sum_CNN_VS']}"
        log.info(f"最优解：t*={ev_o['t']}（占据 {ev_o['n_occupied_bins']}/{n_bins} 个分位格）、"
                 f"CNN_VS 跨度 {r4(ev_o['span'])}、Σ CNN_VS {r4(ev_o['sum'])}；"
                 f"四阶段 status 全 0、mip_gap=0、重解一致")

        # ---- 对照实现（v1.6.0 贪心）与对照的目标值
        heur_sel = greedy_reference(pool, pool_ix, med, n_slots, min_low, chem_cap, risk_cap)
        assert len(heur_sel) == n_slots, f"对照贪心只选出 {len(heur_sel)} 个（配额未填满）"
        heur_ids = [s["cid"] for s in heur_sel]
        ev_h = evaluate(heur_ids)
        # ---- 对照论断的独立核验：可行域内**纯** max Σ CNN_VS（不固定 t/span；
        # bin_rows 中的 t 自由取 0 即不构成约束）。贪心是否「max Σ 口径全局最优」
        # 一律由这次求解的实测值判定，不在文案里写死。
        r_sum_pure = solve(list(base) + list(bin_rows), c_sum())
        sum_pure_max = float(-r_sum_pure.fun)
        greedy_is_max_sum = bool(abs(ev_h["sum"] - sum_pure_max) <= 1e-6)
        log.info(f"纯 max Σ CNN_VS 口径 = {r4(sum_pure_max)}；对照贪心 Σ = "
                 f"{r4(ev_h['sum'])} → 贪心{'' if greedy_is_max_sum else '不'}是该口径全局最优")
        # 字典序支配：对照解本身可行，最优解只能不劣于它；逐阶段条件断言（不可跨阶段比 Σ）
        assert ev_o["t"] >= ev_h["t"], "对照贪心在分层覆盖上胜过 MILP 最优解（建模错误）"
        if ev_o["t"] == ev_h["t"]:
            assert ev_o["span"] >= ev_h["span"] - 1e-9, \
                "同覆盖下限下对照贪心跨度更大（建模错误）"
            if abs(ev_o["span"] - ev_h["span"]) < 1e-9:
                assert ev_o["sum"] >= ev_h["sum"] - 1e-9, \
                    "同为覆盖/跨度最优时对照贪心 Σ CNN_VS 更大（建模错误）"
        n_swaps = len(set(sel_ids) ^ set(heur_ids)) // 2
        assert n_swaps == int(exp_cert["n_swaps_vs_heuristic"]), \
            (f"与对照的名册互换数 {n_swaps} != 规则登记 {exp_cert['n_swaps_vs_heuristic']}；"
             "上游或规则已漂移，须显式复核")
        log.info(f"对照（v1.6.0 贪心）：t={ev_h['t']}、跨度 {r4(ev_h['span'])}、"
                 f"Σ CNN_VS {r4(ev_h['sum'])}；与最优解重叠 "
                 f"{n_slots - n_swaps} 席、互换 {n_swaps} 席")

        # ---- 席位角色（① 跨度端点 ② 覆盖下限席位 ③ 质量兜底）：留一法（roster 内移除该席）
        roles: dict[str, str] = {}
        for cid in sel_ids:
            rest = [c for c in sel_ids if c != cid]
            zv_rest = np.array([float(pool_ix.loc[c, "CNN_VS"]) for c in rest])
            span_rest = float(zv_rest.max() - zv_rest.min())
            live_rest = [sum(1 for c in rest if c in bin_mem[b])
                         for b in bin_mem if bin_mem[b]]
            if span_rest < ev_o["span"] - 1e-9:
                roles[cid] = REASON_SPAN
            elif min(live_rest) < ev_o["t"]:
                roles[cid] = REASON_COVER
            else:
                roles[cid] = REASON_FILL
        n_cover = sum(1 for c in sel_ids if roles[c] == REASON_COVER)
        n_span = sum(1 for c in sel_ids if roles[c] == REASON_SPAN)
        n_fill = sum(1 for c in sel_ids if roles[c] == REASON_FILL)
        log.info(f"席位角色（留一法）：跨度端点 {n_span}、覆盖下限 {n_cover}、质量兜底 {n_fill}")

        # ---- 反事实重解（更强的口径）：把该席从池中禁掉后，覆盖/跨度是否仍可达标
        # 留一法说明「这一席在已定名册里承担什么」；反事实重解说明「这一席是否不可替代」。
        # 可行域为空（status=2）计为覆盖不可替代的**最强形式**并单独计数披露：它来自
        # 化学型帽与名额的零冗余咬合，而非求解器或建模问题。
        cf_cover = cf_span = cf_infeasible = 0
        for cid in sel_ids:
            i = idx_of[cid]
            ban = np.zeros((1, n_var)); ban[0, i] = 1.0
            ban_row = [(ban, np.array([-np.inf]), np.array([0.0]))]
            r1 = solve(list(base) + list(bin_rows) + ban_row, c_t(), allow_infeasible=True)
            if r1 is None:
                # 连「凑满 n_candidates 席」都不可能 ⇒ 覆盖与跨度自然都不可达
                cf_cover += 1
                cf_infeasible += 1
                continue
            if float(round(r1.x[I_T])) < ev_o["t"]:
                cf_cover += 1
                continue
            t_eq = np.array([[1.0 if k == I_T else 0.0 for k in range(n_var)]])
            r2 = solve(list(base) + list(bin_rows) + ban_row
                       + [(t_eq, np.array([float(ev_o["t"])]), np.array([float(ev_o["t"])]))],
                       c_span(), allow_infeasible=True)
            if r2 is None:
                cf_span += 1
                cf_infeasible += 1
                continue
            if float(-r2.fun) < ev_o["span"] - 1e-9:
                cf_span += 1
        cf_replaceable = len(sel_ids) - cf_cover - cf_span
        log.info(f"反事实重解（禁掉该席后仍可达标 ⇒ 可替代）：覆盖不可替代 {cf_cover}、"
                 f"跨度不可替代 {cf_span}、可替代 {cf_replaceable}；其中可行域为空 "
                 f"{cf_infeasible} 席（结构性锁定，非求解失败）")

        def compliance(members) -> dict:
            """四条硬约束的实测值 + 信息量代理（证书/自检共用，独立于求解器）。"""
            idxs = [idx_of[c] for c in members]
            zv = np.array([zvals[i] for i in idxs])
            sub = Tc[np.ix_(idxs, idxs)]
            mx = float(np.max(sub[~np.eye(len(idxs), dtype=bool)])) if len(idxs) > 1 else 0.0
            chem = pd.Series([pool_ix.loc[c, "chemotype"] for c in members]).value_counts()
            chem = {k: int(v) for k, v in chem.items()}
            risk = int(sum(int(pool_ix.loc[c, "risk_am1bcc"]) for c in members))
            below = int(sum(1 for v in zv if v < med))
            return {"n_selected": len(members), "chemotype_counts": chem, "risk_count": risk,
                    "below_median_48_count": below, "max_pairwise_tc": r6(mx),
                    "min_bin_count": evaluate(members)["t"],
                    "n_occupied_bins": evaluate(members)["n_occupied_bins"],
                    "slope_info_sum_sq_dev": r6(float(np.sum((zv - zv.mean()) ** 2))),
                    "sigma_CNN_VS": r6(float(zv.std(ddof=1))),
                    "pass": {"n_selected": len(members) == n_slots,
                             "chem_cap": max(chem.values()) <= chem_cap,
                             "risk_cap": risk <= risk_cap,
                             "below_median_min": below >= min_low,
                             "pairwise_tc": mx < tc_cap}}

        # ---- 批次表（候选按 CNN_VS 降序，同分按 cid；对照固定 2 位）
        keep_cols = ["cid", "pick_order", "smiles", "CNN_VS", "qed", "MW", "cLogP", "TPSA",
                     "HBD", "HBA", "RotB", "cns_mpo", "cns_mpo_lo80", "cns_mpo_hi80",
                     "chemotype", "scaffold", "scaffold_generic", "anchor_hits_4A",
                     "n_contacts_4A", "pose_gate", "maxTc_actives", "nearest_active_id",
                     "herg_prob", "herg_gate", "risk_am1bcc", "basic_pka", "clogd74",
                     "clogd74_lower_80", "clogd74_upper_80"]
        batch = cand[cand.cid.isin(set(sel_ids))].copy()
        batch = batch.sort_values(["CNN_VS", "cid"], ascending=[False, True], kind="mergesort")
        batch = batch[[c for c in keep_cols if c in batch.columns]].copy()
        batch = batch.merge(xm[["cid", "hERG_cross_ens", "m2_pass_ens"]], on="cid", how="left")
        _ann = cand[["cid", "herg_cross_model_flag", "md_eligibility_note"]] \
            if "md_eligibility_note" in cand.columns else cand[["cid", "herg_cross_model_flag"]]
        batch = batch.merge(_ann, on="cid", how="left")
        if "md_eligibility_note" not in batch.columns:
            batch["md_eligibility_note"] = ""
        batch["candidate_type"] = "candidate"
        batch["control_role"] = ""
        batch["selection_reason"] = batch.cid.map(roles)
        # 覆盖分位（三轴各 1 列，便于人读与报告转述；对应证书里的分位边界）
        for col, name in axes:
            _vals = batch[col].to_numpy(float)
            _edges = np.array(bin_edges[name])
            batch[f"bin_{name}"] = ["Q" + str(int(np.searchsorted(_edges, v, side="left")) + 1)
                                    for v in _vals]
        batch["coverage_bins"] = [";".join(f"{nm}={batch.loc[i, 'bin_' + nm]}"
                                           for _c, nm in axes)
                                  for i in batch.index]
        ctrl_out = ctrl.copy()
        if "smiles" not in ctrl_out.columns:
            ctrl_out["smiles"] = ""
        for c in batch.columns:
            if c not in ctrl_out.columns:
                ctrl_out[c] = np.nan
        ctrl_out["candidate_type"] = "negative_control"
        ctrl_out["selection_reason"] = REASON_CTRL
        batch = pd.concat([batch, ctrl_out[batch.columns]], ignore_index=True)
        batch["batch_order"] = range(1, len(batch) + 1)
        assert len(batch) == n_batch, f"批次规模 {len(batch)} != 规则 {n_batch}"

        # ---- 配额与多样性自检（硬断言，防止规则与实现漂移）
        b_cand = batch[batch.candidate_type == "candidate"]
        chem_cnt = b_cand.chemotype.value_counts().to_dict()
        assert max(chem_cnt.values()) <= chem_cap, f"化学型超限：{chem_cnt}"
        risk_n = int(b_cand.risk_am1bcc.sum())
        assert risk_n <= risk_cap, f"风险档超限：{risk_n}"
        assert b_cand.scaffold_generic.nunique() == len(b_cand), "批次内通用骨架重复"
        assert b_cand.scaffold.nunique() == len(b_cand), "批次内 Bemis-Murcko 骨架重复"
        n_low_final = int((b_cand.CNN_VS < med).sum())
        assert n_low_final >= min_low, f"下半区候选 {n_low_final} < 规则 {min_low}"
        cids = b_cand.cid.tolist()
        tc_max, tc_pair = 0.0, None
        for i in range(len(cids)):
            for j in range(i + 1, len(cids)):
                t = DataStructs.TanimotoSimilarity(fps[cids[i]], fps[cids[j]])
                if t > tc_max:
                    tc_max, tc_pair = t, (cids[i], cids[j])
        assert tc_max < tc_cap, f"批次内成对 Tc {tc_max} ≥ {tc_cap}"
        nan_cns = int(b_cand.cns_mpo.isna().sum())
        assert nan_cns == 0, f"CNS MPO 缺失 {nan_cns} 行（决策卡列丢失）"
        assert "hERG_cross_ens" in batch.columns and \
            batch[batch.candidate_type == "candidate"].hERG_cross_ens.notna().all(), \
            "第二模型 hERG 列缺失"
        assert set(e["cid"] for e in excl).isdisjoint(set(b_cand.cid)), "被剔除分子混入批次"
        # 证书 vs 批次表交叉核对（两条独立路径）
        ev_csv = evaluate(cids)
        assert ev_csv["t"] == ev_o["t"] and abs(ev_csv["span"] - ev_o["span"]) < 1e-9, \
            "批次表与求解结果不一致"

        write_csv(OUT / "tables/md_batch20.csv", batch)
        write_csv(OUT / "tables/md_batch_excluded.csv", pd.DataFrame(excl))

        spans = {"CNN_VS_min": r4(b_cand.CNN_VS.min()), "CNN_VS_max": r4(b_cand.CNN_VS.max()),
                 "CNN_VS_span": r4(b_cand.CNN_VS.max() - b_cand.CNN_VS.min()),
                 "cLogP_min": r4(b_cand.cLogP.min()), "cLogP_max": r4(b_cand.cLogP.max()),
                 "cns_mpo_min": r4(b_cand.cns_mpo.min()), "cns_mpo_max": r4(b_cand.cns_mpo.max()),
                 "maxTc_min": r4(b_cand.maxTc_actives.min()),
                 "maxTc_max": r4(b_cand.maxTc_actives.max())}
        level1 = float(mb["md_cost"]["level1_10ns_cardday_per_molecule"])
        level2 = float(mb["md_cost"]["level2_150ns_cardday_per_molecule"])

        swapped_out = [c for c in heur_ids if c not in set(sel_ids)]
        swapped_in = [c for c in sel_ids if c not in set(heur_ids)]
        roster_diff = [
            {"direction": "被替换（对照贪心）", "cid": c,
             "CNN_VS": r4(float(pool_ix.loc[c, "CNN_VS"])),
             "bin_CNN_VS": "Q" + str(int(np.searchsorted(
                 np.array(bin_edges["CNN_VS"]), float(pool_ix.loc[c, "CNN_VS"]),
                 side="left")) + 1),
             "cLogP": r4(float(pool_ix.loc[c, "cLogP"])),
             "cns_mpo": r4(float(pool_ix.loc[c, "cns_mpo"])),
             "chemotype": pool_ix.loc[c, "chemotype"]} for c in swapped_out] + [
            {"direction": "新增（v1.7.0 最优解）", "cid": c,
             "CNN_VS": r4(float(pool_ix.loc[c, "CNN_VS"])),
             "bin_CNN_VS": "Q" + str(int(np.searchsorted(
                 np.array(bin_edges["CNN_VS"]), float(pool_ix.loc[c, "CNN_VS"]),
                 side="left")) + 1),
             "cLogP": r4(float(pool_ix.loc[c, "cLogP"])),
             "cns_mpo": r4(float(pool_ix.loc[c, "cns_mpo"])),
             "chemotype": pool_ix.loc[c, "chemotype"]} for c in swapped_in]

        comp_o = compliance(sel_ids)
        comp_h = compliance(heur_ids)
        rep = {
            "rule_version": rule["rule_version"], "batch_size": n_batch,
            "n_candidates": len(b_cand), "n_controls": n_expected_ctrl,
            "source_list": str(LIST_MAIN.relative_to(ROOT)),
            "source_list_sha256": sha256_file(LIST_MAIN),
            "hERG_cross_model_table_sha256": sha256_file(XMODEL),
            "pool": {"n_listed_candidates": len(cand),
                     "n_excluded_hERG_cross_model_fail": (0 if _retain is False
                                                          else int(len(fail))),
                     "n_retained_hERG_cross_model_fail": (int(len(fail))
                                                          if _retain is False else 0),
                     "hERG_cross_model_policy": ("retain_and_annotate" if _retain is False
                                                 else "exclude"),
                     "n_excluded_generic_scaffold_dup": int(len(dup)),
                     "n_eligible": len(pool),
                     # 化学型帽 vs 名额的咬合（决定可行域冗余与反事实可解性）
                     "chemotype_counts": chem_counts,
                     "fill_capacity_upper_bound": fill_ub,
                     "slot_slack": slot_slack,
                     "n_structurally_forced_seats": n_forced},
            "caps": {"pairwise_tc_max": tc_cap,
                     "chem_cap": chem_cap, "risk_cap": risk_cap, "scaffold_cap": 1,
                     "chemotype_counts": {k: int(v) for k, v in chem_cnt.items()},
                     "risk_count": risk_n},
            "stratification": {"reference_median_CNN_VS": r4(med),
                               "min_below_median": min_low,
                               "n_below_median_selected": n_low_final},
            "diversity": {"n_unique_scaffold": int(b_cand.scaffold.nunique()),
                          "n_unique_generic_scaffold": int(b_cand.scaffold_generic.nunique()),
                          "max_pairwise_tc": r4(tc_max), "max_pairwise_tc_pair": tc_pair},
            "spans": spans,
            "objective": {
                "kind": obj["kind"],
                "solver": "scipy.optimize.milp（HiGHS 后端）",
                "scipy_version": scipy.__version__,
                "coverage_axes": [{"column": c, "label": nm} for c, nm in axes],
                "n_bins_per_axis": K, "n_bins": n_bins, "empty_bins": empty_bins,
                "bin_reference": obj["bin_reference"],
                "bin_edges": {k: [r6(e) for e in v] for k, v in bin_edges.items()},
                "objective_convention": ("各阶段 objective 为求解器返回的最小化形式（fun）；"
                                         "实际最大化目标 = −objective"),
                "lexicographic": [{"stage": j["stage"], "name": j["name"],
                                   "objective": r6(j["objective"]), "status": j["status"],
                                   "mip_gap": j["mip_gap"], "nodes": j["nodes"]}
                                  for j in journal],
                "reference_impl": ("heuristic_reference（v1.6.0 贪心填充，fill_order 见规则）"
                                   "——它不驱动选批，仅作对照")},
            "optimality": {
                "certificate": True,
                "t_star": ev_o["t"], "span": r6(ev_o["span"]),
                "sum_CNN_VS": r6(ev_o["sum"]),
                "n_occupied_bins": ev_o["n_occupied_bins"], "n_bins": n_bins,
                "bin_counts": {f"{b[0]}|Q{b[1] + 1}": {"pool_n": len(bin_mem[b]),
                                                      "optimal": ev_o["counts"][b],
                                                      "heuristic": ev_h["counts"][b]}
                               for b in bin_mem},
                "roster": sel_ids,
                "resolved_roster_consistent": bool(resolved_consistent),
                "selection_role_counts": {REASON_SPAN: n_span, REASON_COVER: n_cover,
                                          REASON_FILL: n_fill},
                "counterfactual_resolve": {
                    "n_irreplaceable_for_coverage": cf_cover,
                    "n_irreplaceable_for_span": cf_span,
                    "n_replaceable": cf_replaceable,
                    "n_infeasible_after_ban": cf_infeasible,
                    "note": ("留一法说明该席在名册内承担什么；本块说明「禁掉该席后是否仍可达标」。"
                             "可行域为空（n_infeasible_after_ban）计为覆盖不可替代的最强形式，"
                             "来源是化学型帽与名额的零冗余咬合（可填上限=名额），不是求解失败。"
                             f"实测：覆盖不可替代 {cf_cover}、跨度不可替代 {cf_span}、"
                             f"可替代 {cf_replaceable}；结论一律以上方计数为准，"
                             "不在本文案里预设「有富余/端点约束」的口径。")},
                "compliance": comp_o},
            "heuristic_reference": {
                "role": "v1.6.0 贪心填充（对照实现，不驱动选批）",
                "fill_order_used": mb["fill_order"],
                "t": ev_h["t"], "span": r6(ev_h["span"]), "sum_CNN_VS": r6(ev_h["sum"]),
                "n_occupied_bins": ev_h["n_occupied_bins"],
                "roster": heur_ids,
                "overlap_with_optimal": n_slots - n_swaps, "n_swaps": n_swaps,
                "roster_diff": roster_diff,
                "compliance": comp_h,
                "pure_max_sum_CNN_VS": r6(sum_pure_max),
                "greedy_is_max_sum_optimum": greedy_is_max_sum,
                "note": (f"可行域内纯 max Σ CNN_VS = {r4(sum_pure_max)}（独立求解，不固定 t/span）；"
                         f"对照贪心 Σ = {r4(ev_h['sum'])} → 贪心"
                         + ("恰是该口径的全局最优解——本版修的是目标函数，不是算法"
                            if greedy_is_max_sum else
                            f"**不是**该口径的全局最优（差 {r4(sum_pure_max - ev_h['sum'])}）——"
                            "其价值在于给出可复现的贪心基线，而非最优性"))},
            "md_cost_estimate": {
                "note": "按 06_deliverables/SOP-05_downstream_MD_plan.md §3 的单分子卡时口径线性外推，仅供预算参考",
                "level1_cardday": r4(n_batch * level1),
                "level2_cardday_if_all_pass": r4(n_batch * level2),
                "level2_cardday_at_65pct_pass": r4(round_half(0.65 * n_batch) * level2),
                "level1_assumption": f"{level1} 卡日/分子/10 ns",
                "level2_assumption": f"{level2} 卡日/分子/150 ns（通过率按 60–70% 估计）"},
            "selection_reason_counts": {REASON_COVER: n_cover, REASON_SPAN: n_span,
                                         REASON_FILL: n_fill},
            "note": ("分析步：只做机时分配与留痕；不改名单、不改门。「最优」的限定口径是 "
                     "step6_md_batch.objective（信息量/设计优先，v1.7.0）：Σ CNN_VS 反而略低于"
                     "质量口径，价值在于让一级 MD 能检验分数-稳定性关系并覆盖性质空间；"
                     "低分位分子失败概率更高，这是要测量的量，不得事后解读为「选坏了」。"
                     "MD 后若换人须改规则并重跑整链"),
        }
        write_json(OUT / "reports/md_batch.json", rep)
        log.info(f"批次配额：化学型 {chem_cnt}、风险档 {risk_n}、下半区 {n_low_final}；"
                 f"批内成对 Tc 上限 {r4(tc_max)}（{tc_pair[0] if tc_pair else '-'} × "
                 f"{tc_pair[1] if tc_pair else '-'}）")
        log.info(f"CNN_VS 跨度 {spans['CNN_VS_min']}–{spans['CNN_VS_max']}"
                 f"（宽 {spans['CNN_VS_span']}）")

        # ---- 图：批次在名单中的位置 + 批次构成 + 分位覆盖（最优 vs 对照）
        plots.setup_cjk()
        import matplotlib.pyplot as plt
        fig, axes_ = plt.subplots(1, 3, figsize=(17.0, 4.4))
        ax = axes_[0]
        rank = cand.sort_values("CNN_VS", ascending=False).reset_index(drop=True)
        inb = rank.cid.isin(set(sel_ids))
        excl_ids = [e["cid"] for e in excl]
        is_ex = rank.cid.isin(excl_ids)
        ax.scatter(rank.index[is_ex] + 1, rank.CNN_VS[is_ex], s=42, c="#b0b7c3",
                   label=f"未入批（剔除 {int(is_ex.sum())}）", zorder=3)
        ax.scatter(rank.index[~inb & ~is_ex] + 1, rank.CNN_VS[~inb & ~is_ex], s=22,
                   c="#8fa8c8", alpha=0.6, label=f"未入批（保留 {int((~inb & ~is_ex).sum())}）",
                   zorder=2)
        ax.scatter(rank.index[inb] + 1, rank.CNN_VS[inb], s=52, c="#1f6fb2",
                   label=f"MD 批次 {int(inb.sum())}", zorder=4)
        ax.axhline(med, ls="--", lw=1, c="#c0504d")
        ax.text(1, med, f" 中位 {med:.3f}", color="#c0504d", fontsize=8, va="bottom")
        ax.set_xlabel(f"名单内 CNN_VS 名次（{len(cand)} 候选）")
        ax.set_ylabel("CNN_VS")
        ax.set_title("MD 批次在名单中的分布")
        ax.legend(fontsize=8, loc="lower left")
        ax = axes_[1]
        labels = list(chem_cnt.keys()) + [f"对照 {n_expected_ctrl}"]
        vals = list(chem_cnt.values()) + [n_expected_ctrl]
        ax.bar(range(len(vals)), vals,
               color=["#1f6fb2", "#5b8fc9", "#9dbbdd", "#b0b7c3"][:len(vals)])
        for i, v in enumerate(vals):
            ax.text(i, v + 0.1, str(v), ha="center", fontsize=9)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels([l.replace(" + ", "+\n") for l in labels], fontsize=8)
        ax.set_ylabel("分子数")
        ax.set_title(f"批次构成（化学型上限 {chem_cap}）")
        ax = axes_[2]
        bl = list(bin_mem)
        xs = np.arange(len(bl))
        ax.bar(xs - 0.26, [len(bin_mem[b]) for b in bl], 0.26, color="#cfd6de",
               label="资格池")
        ax.bar(xs, [ev_h["counts"][b] for b in bl], 0.26, color="#d08c1e",
               label="对照（v1.6.0 贪心）")
        ax.bar(xs + 0.26, [ev_o["counts"][b] for b in bl], 0.26, color="#1f6fb2",
               label="最优解（v1.7.0）")
        ax.axhline(ev_o["t"], ls="--", lw=1, c="#4f8a63")
        ax.text(0.1, ev_o["t"] + 0.12, f"覆盖下限 t*={ev_o['t']}", color="#3d6b4d",
                fontsize=8)
        ax.set_xticks(xs)
        ax.set_xticklabels([f"Q{k[1] + 1}" for k in bl], fontsize=8)
        for a_i, (_c, nm) in enumerate(axes):
            ax.text(a_i * K + (K - 1) / 2, -0.155 * max(len(bin_mem[b]) for b in bl),
                    nm, ha="center", va="top", fontsize=9, transform=ax.transData,
                    color="#2347C4" if a_i == 0 else "#4a4a4a")
            if a_i:
                ax.axvline(a_i * K - 0.5, ls=":", lw=0.8, c="#c3cad3")
        ax.set_ylabel("席位数")
        ax.set_title(f"分位覆盖：{len(axes)} 轴 × {K} 层（{n_bins} 格）", pad=10)
        ax.set_ylim(0, max(len(bin_mem[b]) for b in bl) * 1.32)
        ax.legend(fontsize=8, loc="upper right", frameon=False, ncol=3,
                  columnspacing=0.8, handlelength=1.1)
        fig.tight_layout()
        plots.save(fig, OUT / "figures/fig_md_batch.png")
        plt.close(fig)

        # ---- 人读版说明
        t = []
        A = t.append
        A("# MD 批次（20 个）选择说明")
        A("")
        A(f"- 来源名单：`{rep['source_list']}`（{len(cand)} 候选 + {n_expected_ctrl} 对照，主口径 "
          f"{mb['list']['main_arm']}），sha256 `{rep['source_list_sha256'][:16]}…`")
        A(f"- 批次规模：**{n_batch}** = {len(b_cand)} 候选 + {n_expected_ctrl} 阴性对照"
          f"（NEG1/NEG2，固定占 2 位）")
        A(f"- 规则版本：{rule['rule_version']}（阈值与目标全部来自 `config/frozen_rule.json → "
          f"step6_md_batch`，脚本不内嵌阈值）")
        A("")
        A("## 资格与配额")
        A("")
        if _retain is False:
            A(f"1. **不再剔除**第二 hERG 模型（集成口径）不过门的候选：**{len(fail)}** 个 "
              + "、".join(f"`{c}`" for c in fail.cid)
              + f"（hERG_cross_ens ≥ {_thr}）**保留入池并逐分子标注「模型依赖」**。"
                "理由：该判定来自与主模型 κ≈0（近乎随机一致）的第二模型，"
                "用一个未校准信号决定机时分配在纪律上最弱；机时代价仅 +"
              + f"{len(fail) * level1:.2f} 卡日量级（{level1} 卡日/分子/10 ns）。"
                "**本层为机时分配，不产生安全性结论**；两个模型一致判高风险者才应优先复核。")
        else:
            A(f"1. 剔除第二 hERG 模型（集成口径）不过门的候选 **{len(fail)}** 个："
              + "、".join(f"`{c}`" for c in fail.cid)
              + "（MD 回答结合稳定性，不回答心脏安全性；安全性未决者不占 MD 预算，"
                "应先用正交方法解决）")
        A(f"2. 通用骨架去重剔除 **{len(dup)}** 个同系物占位者（保留 CNN_VS 高者）："
          + "、".join(f"`{c}`" for c in dup.cid))
        A(f"3. 配额：成对 ECFP4 Tc < {tc_cap}、骨架 ≤1、"
          f"化学型 ≤ round(8×{n_batch}/20)={chem_cap}、AM1-BCC 风险档 ≤ "
          f"round(3×{n_batch}/20)={risk_cap}")
        A(f"4. 分层：≥{min_low} 个候选来自 CNN_VS 下半区（相对 {len(cand)} 候选的中位 "
          f"{r4(med)}）；实测 **{n_low_final}** 个")
        A("")
        A("## 目标与最优性证书（v1.7.0）")
        A("")
        A(f"在配额可行域内用 MILP（`scipy.optimize.milp`，HiGHS 后端；scipy "
          f"{scipy.__version__}）按**字典序**求解，可证明最优：")
        A("")
        A("| 阶段 | 目标 | 解值（最大化目标） | status | mip_gap | 分支节点 |")
        A("|---|---|---|---|---|---|")
        for j in journal:
            A(f"| {j['stage']} | {j['name']} | {-j['objective']:.6g} | {j['status']} | "
              f"{j['mip_gap']} | {j['nodes']} |")
        A("")
        A("> 求解器内部以最小化形式 `fun = −目标` 求解（证书 `objective.lexicographic[].objective` "
          "记的即 `fun`），上表已换算为最大化目标的解值。")
        A("")
        A(f"- 阶段 1 的分层覆盖：{len(axes)} 条性质轴（"
          + "、".join(f"{nm}" for _c, nm in axes) + f"）各按**资格池自分布**切 {K} 分位，"
          f"共 {n_bins} 格；最优解使**每格至少 {ev_o['t']} 席**（{ev_o['n_occupied_bins']}/"
          f"{n_bins} 格被占据）")
        A(f"- 阶段 2 的 CNN_VS 跨度：**{r4(ev_o['span'])}**（分数轴量程，"
          f"资格池全量程 {r4(zmax - zmin)}）；阶段 3 的 Σ CNN_VS：**{r4(ev_o['sum'])}**")
        A(f"- 重解一致：{resolved_consistent}（同一输入两次求解名册逐位相同）；"
          f"各阶段目标值均等于名册实算值（跨度用独热指派线性化，防「假跨度」）")
        A("")
        A("**分位覆盖（资格池 / 对照贪心 / 最优解）**")
        A("")
        A("| 分位 | 池内分子数 | 对照贪心席位 | 最优解席位 |")
        A("|---|---|---|---|")
        for b in bin_mem:
            A(f"| {b[0]} Q{b[1] + 1} | {len(bin_mem[b])} | {ev_h['counts'][b]} | "
              f"{ev_o['counts'][b]} |")
        A("")
        A("## 与 v1.6.0 贪心填充对照（同一步内的对照实现）")
        A("")
        A(f"- 对照贪心（① 化学型覆盖锚点 → ② 下半区补齐至 ≥{min_low} 席 → "
          f"③ CNN_VS 降序填满，受 caps 限制）实测：t={ev_h['t']}、"
          f"跨度 {r4(ev_h['span'])}、Σ CNN_VS {r4(ev_h['sum'])}、"
          f"占据 {ev_h['n_occupied_bins']}/{n_bins} 格")
        A(f"- 与最优解重叠 **{n_slots - n_swaps}** 席、互换 **{n_swaps}** 席"
          f"（代价 Σ CNN_VS {r4(ev_o['sum'] - ev_h['sum'])}，换来跨度 "
          f"+{r4(ev_o['span'] - ev_h['span'])}、下半区席位 "
          f"{comp_h['below_median_48_count']}→{comp_o['below_median_48_count']}）")
        if greedy_is_max_sum:
            A("- **独立求解确认：该贪心恰是纯「max Σ CNN_VS」口径的全局最优解**——"
              f"（纯口径 Σ = {r4(sum_pure_max)} = 贪心 Σ）本版修的是目标函数，不是算法。")
        else:
            A(f"- **独立求解显示该贪心不是纯「max Σ CNN_VS」口径的全局最优**：纯口径 "
              f"Σ = {r4(sum_pure_max)}、贪心 Σ = {r4(ev_h['sum'])}（差 "
              f"{r4(sum_pure_max - ev_h['sum'])}）；本对照只作可复现基线，不作最优性依据。")
        A("")
        A("| 方向 | cid | CNN_VS | 分位 | cLogP | CNS MPO | 化学型 |")
        A("|---|---|---|---|---|---|---|")
        for d in roster_diff:
            A(f"| {d['direction']} | `{d['cid']}` | {d['CNN_VS']} | {d['bin_CNN_VS']} | "
              f"{d['cLogP']} | {d['cns_mpo']} | {d['chemotype']} |")
        A("")
        A("## 批次构成（按 CNN_VS 降序，对照在末 2 位）")
        A("")
        A("| # | cid | 类型 | 名单序 | CNN_VS | QED | CNS MPO | maxTc | hERG | 第二模型 hERG | 模型依赖标注 | 风险档 | 覆盖分位 | 入选理由 |")
        A("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in batch.itertuples():
            typ = "候选" if r.candidate_type == "candidate" else r.control_role
            A(f"| {r.batch_order} | `{r.cid}` | {typ} | "
              f"{'—' if r.candidate_type != 'candidate' else int(r.pick_order)} | "
              f"{'—' if r.candidate_type != 'candidate' else f'{r.CNN_VS:.3f}'} | "
              f"{'—' if r.candidate_type != 'candidate' else f'{r.qed:.3f}'} | "
              f"{'—' if r.candidate_type != 'candidate' else f'{r.cns_mpo:.2f}'} | "
              f"{'—' if r.candidate_type != 'candidate' else f'{r.maxTc_actives:.3f}'} | "
              f"{'—' if r.candidate_type != 'candidate' else f'{r.herg_prob:.3f}'} | "
              f"{'—' if r.candidate_type != 'candidate' else f'{r.hERG_cross_ens:.3f}'} | "
              f"{'—' if r.candidate_type != 'candidate' else (str(r.md_eligibility_note) or '')} | "
              f"{'—' if r.candidate_type != 'candidate' else int(r.risk_am1bcc)} | "
              f"{'—' if r.candidate_type != 'candidate' else r.coverage_bins} | "
              f"{r.selection_reason} |")
        A("")
        A(f"席位角色（留一法：在已定名册内移除该席后 t 是否下降 / 跨度是否收窄）："
          f"跨度端点 **{n_span}**、覆盖下限 **{n_cover}**、质量兜底 **{n_fill}**。"
          f"另做**反事实重解**（把该席从资格池禁掉后重新求解）：覆盖不可替代 {cf_cover}、"
          f"跨度不可替代 {cf_span}、可替代 {cf_replaceable}（其中可行域为空 "
          f"**{cf_infeasible}** 席）。可行域为空来自化学型帽与名额的零冗余咬合"
          f"（池内化学型 {chem_counts}、帽 {chem_cap} ⇒ 可填上限 {fill_ub} = 名额 {n_slots}），"
          "属**结构性锁定**，不是求解器失败；具体哪几席不可替代以批次表 selection_reason 与"
          " 上述计数为准。")
        A("")
        A(f"批次自检：化学型分布 {chem_cnt}；风险档命中 {risk_n}（上限 {risk_cap}）；"
          f"唯一骨架 {int(b_cand.scaffold.nunique())}；批内成对 Tc 上限 "
          f"**{r4(tc_max)}**（名单整体成对 Tc 上限见 selection_summary）；"
          f"CNN_VS 跨度 {spans['CNN_VS_min']}–{spans['CNN_VS_max']}"
          f"（宽 {spans['CNN_VS_span']}）；cLogP {spans['cLogP_min']}–{spans['cLogP_max']}、"
          f"CNS MPO {spans['cns_mpo_min']}–{spans['cns_mpo_max']}、"
          f"maxTc {spans['maxTc_min']}–{spans['maxTc_max']}。")
        A("")
        A("## 执行建议（判据须先预注册）")
        A("")
        A("- **先跑 2 个对照**（各 10 ns）确认一级判据有区分度（对照应漂出）。对照若不漂，"
          "整批判读无效，及时止损。")
        A("- 一级判据（10 ns/分子）：锚点接触保留 ≥3/4 且配体重原子 RMSD < 0.3 nm；"
          "二级（100–200 ns）：锚点 4/4 保留率 ≥80%。`gmx_MMPBSA` 缺失时二级退化为"
          "「锚点保留率 + 氢键占据率 + 能量分解」，不得表述为自由能排序。")
        A(f"- 机时（按下游计划单分子卡时线性外推，仅供预算参考）：一级 "
          f"≈{rep['md_cost_estimate']['level1_cardday']} 卡日"
          f"（{rep['md_cost_estimate']['level1_assumption']}）；二级按 60–70% 通过率估算 "
          f"≈{rep['md_cost_estimate']['level2_cardday_at_65pct_pass']} 卡日"
          f"（全通过则 ≈{rep['md_cost_estimate']['level2_cardday_if_all_pass']} 卡日）"
          f"（{rep['md_cost_estimate']['level2_assumption']}）。")
        A("- 边界分子（hERG 0.44–0.49、CNS MPO<4 者）优先**加重复**而不是加新分子。")
        A("- 分位覆盖意图：本批刻意纳入最低分位（CNN_VS Q1）分子——它们 MD 失败概率更高，"
          "这正是要**测量**的量（分数-稳定性关系需要量程与分位两端的样本），"
          "不得事后解读为「选坏了」。")
        A("")
        A("> 本步为分析步：不修改名单与任何门。「最优」的限定口径是 "
          "`step6_md_batch.objective`（信息量/设计优先，v1.7.0），**不是**「这批分子比其余候选"
          "更好」——最优解的 Σ CNN_VS 反而略低。MD 后若需换人，须改规则并重跑整链"
          "（`python3 run_all.py --freeze` → 再跑至不动点），禁止手工编辑名单文件。")
        A("")
        write_text(OUT / "MD_BATCH.md", "\n".join(t))
        log.info(f"MD 批次 → 07_md_batch/tables/md_batch20.csv（{n_batch} 个）、"
                 f"MD_BATCH.md、reports/md_batch.json（含最优性证书）、"
                 f"figures/fig_md_batch.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
