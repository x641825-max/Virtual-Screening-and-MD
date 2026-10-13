#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 17 · 合并两分支 MD 批次 → 选 20 席最终 MD 候选（派生分析层）

为什么需要这一步
  SOP_5（目录化合物 ZINC/MCULE）与 SOP_ML_run05（ML 生成分子 SCALE-）是同一套规则
  v1.7.0 下的两条平行链，各自已用 MILP 字典序优化出 MD 一级批次（18 席 / 17 席）。
  两份名册按 InChIKey **结构零重叠**（唯二相同的两个是共用阴性对照），因此合并是
  「两个互斥分子来源取并集」，不是两套方法的重复计数。

本层做的事
  把两条链的 MD 一级批次候选并成 35 个池，在合并池上**重解**一次名额 20 的批次，
  目标与约束沿用 frozen_rule.step6_md_batch（不改规则、不改门、不改上游名单）。

与主链的三处口径差异（本层登记，不写回冻结规则）
  1. 名额：n_candidates 由 18 改为 20（用户指定），名册 22 行 = 20 候选 + 2 共用对照。
  2. 配额：化学型帽取冻结登记值 8、风险档取 3。规则自带的按席位缩放式
     round(8N/20) 在 N=22 下会给出 9，本层取**更严**的一侧，不新增登记值。
  3. 分层参照：「CNN_VS 下半区」的中位参照取**合并候选池（35）中位**，
     与规则「用本层候选名单自身中位」的行为一致（原链为 48 候选的 7.414）。

必须随结果转述的边界
  · 跨来源 CNN_VS 不可直接比：SOP_5 的池是目录池按 CNN_VS 截断出的 top-4,990
    （下限 6.554），ML 池是生成库全范围（0.963–8.536）。阶段 3 的 Σ CNN_VS 会系统性
    偏向目录来源 —— 这是选池口径差异而非化学优劣。本层因此额外跑一个
    「来源池内百分位」口径的敏感性对照（sensitivity），只作留痕，不改主名册。
  · hERG 只是单一模型的分类线，只能说「在 ADMET-AI 口径下过门」；第二模型不过门者
    保留并逐分子标注「模型依赖」（v1.6.0 起口径），本层不产生安全性结论。
  · 本层是机时分配层，不产生有效性结论。

输出：tables/final_md_20.csv、tables/final_md_excluded.csv、reports/final_md.json、
      reports/FINAL_MD_20.md（图由 18_final_md_figures.py 另出）
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy  # noqa: F401  （版本写入 JSON，便于复现）
from rdkit import DataStructs, RDLogger
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from scipy.optimize import Bounds, LinearConstraint, milp

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
ROOT = _p
sys.path.insert(0, str(ROOT))

from vendor.io_utils import StepLog, read_json, r4, r6, sha256_file, write_csv, write_json, write_text  # noqa: E402

RDLogger.DisableLog("rdApp.*")

OUT = ROOT / "10_final_md"
ML_ROOT = Path(os.environ.get("SOPML_ROOT", str(ROOT.parent / "SOP_ML_run05")))
# 另一条链（SOP_ML_run05）的包根：默认取与本包**同级**的目录，故原工作区布局
# （SOP_5 与 SOP_ML_run05 同级）与本提交包布局
# （src/SOP_5 与 src/SOP_ML_run05 同级）都直接可用；特殊布局用环境变量 SOPML_ROOT 覆盖。

SRC = {
    "SOP_5（目录 ZINC/MCULE）": {
        "batch": ROOT / "07_md_batch/tables/md_batch20.csv",
        "pool": ROOT / "03_gates/tables/candidate_universe.csv",
        "fps": ROOT / "02_druglikeness/tables/fingerprints_ecfp4.npz",
        "synth": ROOT / "09_synth_access/tables/synth_access_list.csv",
        "short": "目录",
    },
    "SOP_ML_run05（ML 生成 SCALE-）": {
        "batch": ML_ROOT / "07_md_batch/tables/md_batch20.csv",
        "pool": ML_ROOT / "03_gates/tables/candidate_universe.csv",
        "fps": ML_ROOT / "02_druglikeness/tables/fingerprints_ecfp4.npz",
        "synth": ML_ROOT / "09_synth_access/tables/synth_access_list.csv",
        "short": "生成",
    },
}
RULE = ROOT / "config/frozen_rule.json"


def _rel(p) -> str:
    """输入登记写**包根相对路径**，不把构建期宿主绝对路径写进报告（赛事附件 5 §六）。"""
    p = Path(p)
    for base, pre in ((ROOT, ""), (ROOT.parent, "../")):
        try:
            return pre + str(p.relative_to(base))
        except ValueError:
            continue
    return p.name

N_SLOTS = 20
N_CTRL = 2
CHEM_CAP = 8          # = frozen_rule.step6_md_batch.caps.chem_cap_value（取更严的一侧）
RISK_CAP = 3          # = frozen_rule.step6_md_batch.caps.risk_cap_value
TC_CAP = 0.65         # = caps.pairwise_tc_max
MIN_LOW = 5           # = stratification.min_below_median_candidates
N_BINS = 4            # = objective.n_bins_per_axis
AXES = [("CNN_VS", "CNN_VS"), ("cLogP", "cLogP"), ("cns_mpo", "CNS MPO")]
OPT_OPTS = {"presolve": True, "disp": False}
MIP_GAP_TOL = 1e-9

REASON_CTRL = "阴性对照（固定 2 位，不占候选名额）"
REASON_COVER = "覆盖下限约束席位（移除此席 t 会下降）"
REASON_SPAN = "跨度端点席位（移除后跨度收窄）"
REASON_FILL = "质量兜底填充（移除后 t 与跨度均不变）"
HEUR_ANCHOR = "化学型覆盖锚点（该类最高分）"
HEUR_LOW = "分数分层（CNN_VS 下半区）"
HEUR_FILL = "高分填充"

KEEP_COLS = ["cid", "pick_order", "smiles", "CNN_VS", "qed", "MW", "cLogP", "TPSA",
             "HBD", "HBA", "RotB", "cns_mpo", "cns_mpo_lo80", "cns_mpo_hi80",
             "chemotype", "scaffold", "scaffold_generic", "anchor_hits_4A",
             "n_contacts_4A", "pose_gate", "maxTc_actives", "nearest_active_id",
             "herg_prob", "herg_gate", "risk_am1bcc", "basic_pka", "clogd74",
             "clogd74_lower_80", "clogd74_upper_80"]


def round_half(x: float) -> int:
    return int(np.floor(x + 0.5))


def bins_for(vals: np.ndarray, ref: np.ndarray, k: int):
    edges = np.array([float(np.percentile(ref, 100 * i / k)) for i in range(1, k)])
    return np.searchsorted(edges, vals, side="left"), edges


def load_fps(path: Path) -> tuple[dict, dict]:
    z = np.load(path, allow_pickle=True)
    bv, bits = {}, {}
    for cid, b in zip(z["cid"], z["bits"]):
        v = ExplicitBitVect(2048)
        v.SetBitsFromList(np.flatnonzero(b).tolist())
        bv[str(cid)] = v
        bits[str(cid)] = b.astype(np.float32)
    return bv, bits


def tanimoto_matrix(bits: np.ndarray) -> np.ndarray:
    inter = bits @ bits.T
    sz = bits.sum(1)
    t = inter / (sz[:, None] + sz[None, :] - inter)
    np.fill_diagonal(t, 0.0)
    return t


def greedy_reference(pool: pd.DataFrame, med: float, n_slots: int, min_low: int) -> list[dict]:
    """v1.6.0 贪心填充（对照实现，逐字沿用其选序）：① 化学型锚点 → ② 下半区补齐 → ③ 高分填充。"""
    sel: list[dict] = []
    used_chem: dict[str, int] = {}
    used_risk = 0

    def can_take(row) -> bool:
        return (used_chem.get(row.chemotype, 0) < CHEM_CAP
                and used_risk + int(row.risk_am1bcc) <= RISK_CAP)

    def take(row, reason: str) -> None:
        nonlocal used_risk
        sel.append({"cid": str(row.cid), "reason": reason})
        used_chem[row.chemotype] = used_chem.get(row.chemotype, 0) + 1
        used_risk += int(row.risk_am1bcc)

    for _chem, grp in pool.groupby("chemotype", sort=True):
        if len(sel) >= n_slots:
            break
        take(grp.iloc[0], HEUR_ANCHOR)
    n_low = sum(1 for s in sel if float(pool.set_index("cid").loc[s["cid"], "CNN_VS"]) < med)
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
    log = StepLog(OUT, "17_final_md_merge")
    rule = read_json(RULE)
    assert rule["rule_version"] == "1.7.0", f"规则版本 {rule['rule_version']} != 1.7.0"
    mb = rule["step6_md_batch"]
    assert int(mb["caps"]["chem_cap_value"]) == CHEM_CAP
    assert int(mb["caps"]["risk_cap_value"]) == RISK_CAP
    assert abs(float(mb["caps"]["pairwise_tc_max"]) - TC_CAP) < 1e-12
    assert int(mb["stratification"]["min_below_median_candidates"]) == MIN_LOW
    assert int(mb["objective"]["n_bins_per_axis"]) == N_BINS

    # ------------------------------------------------------------ 1. 合并两分支批次
    frames, excl = [], []
    for name, spec in SRC.items():
        d = pd.read_csv(spec["batch"])
        d["source_branch"] = name
        frames.append(d)
        log.info(f"{name}：MD 一级批次 {len(d)} 行"
                 f"（候选 {int((d.candidate_type == 'candidate').sum())} + "
                 f"对照 {int((d.candidate_type != 'candidate').sum())}）")
    m = pd.concat(frames, ignore_index=True)
    cand = m[m.candidate_type == "candidate"].copy()
    ctrl = m[m.candidate_type != "candidate"].copy()
    ctrl = ctrl.drop_duplicates(subset=["cid"]).copy()   # 两个共用对照只留一份
    assert len(ctrl) == N_CTRL, f"共用阴性对照 {len(ctrl)} 个 != {N_CTRL}"
    assert cand.cid.is_unique, "合并候选池有重复 cid"
    log.info(f"合并候选池 {len(cand)} 个（目录 {int((cand.source_branch == list(SRC)[0]).sum())}"
             f" + 生成 {int((cand.source_branch == list(SRC)[1]).sum())}）"
             f"；共用阴性对照 {len(ctrl)} 个")

    # ------------------------------------------------------------ 2. 资格判定
    # ① hERG 第二模型：不过门者**保留入池并标注**（v1.6.0 口径），不剔除
    fail = cand[cand.herg_cross_model_flag.fillna("").eq("cross_model_fail")].copy()
    cand["md_eligibility_note"] = cand["md_eligibility_note"].fillna("")
    log.info(f"第二 hERG 模型不过门：{len(fail)} 个 → 保留入池并逐分子标注「模型依赖」")
    # ② 通用骨架去重（保留 CNN_VS 高者）
    cand = cand.sort_values("CNN_VS", ascending=False, kind="mergesort").reset_index(drop=True)
    dup = cand[cand.duplicated("scaffold_generic", keep="first")]
    for r in dup.itertuples():
        keep = cand[cand.scaffold_generic == r.scaffold_generic].iloc[0]
        excl.append({"cid": str(r.cid), "source_branch": r.source_branch,
                     "CNN_VS": r4(float(r.CNN_VS)), "reason": "通用骨架重复（同系物占位）",
                     "detail": f"scaffold_generic 已有更高分者 {keep.cid}"
                               f"（{keep.CNN_VS:.4f} > {r.CNN_VS:.4f}，"
                               f"{keep.source_branch}）"})
    pool = cand[~cand.duplicated("scaffold_generic", keep="first")].reset_index(drop=True)
    log.info(f"资格池 {len(pool)} 个（通用骨架去重剔除 {len(dup)} 个）")

    # ------------------------------------------------------------ 3. 分层参数
    med = float(cand.CNN_VS.median())
    n_bins = N_BINS * len(AXES)
    log.info(f"名额 {N_SLOTS} 候选 + {N_CTRL} 对照 = {N_SLOTS + N_CTRL} 行；"
             f"化学型 ≤{CHEM_CAP}、风险档 ≤{RISK_CAP}、骨架 ≤1、成对 Tc <{TC_CAP}；"
             f"CNN_VS 中位（相对合并候选池 {len(cand)}）={med:.4f}，下半区下限 {MIN_LOW} 个；"
             f"分位 {len(AXES)} 轴 × {N_BINS} 层 = {n_bins} 格")

    ids = [str(c) for c in pool.cid]
    idx_of = {c: i for i, c in enumerate(ids)}
    pool_ix = pool.set_index("cid")
    zvals = pool.CNN_VS.to_numpy(float)

    bin_mem: dict[tuple[str, int], list[str]] = {}
    bin_edges: dict[str, list[float]] = {}
    for col, name in AXES:
        idx, edges = bins_for(pool[col].to_numpy(float), pool[col].to_numpy(float), N_BINS)
        bin_edges[name] = [float(e) for e in edges]
        for b in range(N_BINS):
            bin_mem[(name, b)] = [ids[i] for i in range(len(ids)) if idx[i] == b]
    live_bins = [k for k, v in bin_mem.items() if v]
    empty_bins = [f"{k[0]}|Q{k[1] + 1}" for k, v in bin_mem.items() if not v]
    log.info("分位边界（资格池自分布）：" + "；".join(
        f"{k} {[round(e, 4) for e in v]}" for k, v in bin_edges.items()) +
        f"；有效分位格 {len(live_bins)}（空格 {len(empty_bins)}）")

    # ------------------------------------------------------------ 4. 指纹与成对 Tc
    fps_all, bits_all = {}, {}
    for spec in SRC.values():
        bv, bits = load_fps(spec["fps"])
        fps_all.update(bv)
        bits_all.update(bits)
    missing = [c for c in ids if c not in fps_all]
    assert not missing, f"指纹缺失 {missing}"
    X = np.vstack([bits_all[c] for c in ids])
    Tc = tanimoto_matrix(X)
    tc_conflict = [(i, j) for i in range(len(ids)) for j in range(i + 1, len(ids))
                   if Tc[i, j] >= TC_CAP]
    tc_max_pool = float(Tc.max()) if len(ids) > 1 else 0.0
    log.info(f"资格池成对 ECFP4 Tc 最大 {tc_max_pool:.4f}；Tc ≥ {TC_CAP} 的冲突对 {len(tc_conflict)}")

    # ------------------------------------------------------------ 5. MILP（字典序四阶段）
    n = len(ids)
    n_var = 3 * n + 3
    I_T, I_U, I_L = n, n + 1, n + 2
    I_Y, I_W = n + 3, 2 * n + 3
    Mbig = float(zvals.max() - zvals.min()) + 1.0

    base: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    A = np.zeros((1, n_var)); A[0, :n] = 1.0
    base.append((A, np.array([float(N_SLOTS)]), np.array([float(N_SLOTS)])))
    for c in sorted(pool.chemotype.unique()):
        A = np.zeros((1, n_var))
        A[0, [idx_of[i] for i in ids if pool_ix.loc[i, "chemotype"] == c]] = 1.0
        base.append((A, np.array([-np.inf]), np.array([float(CHEM_CAP)])))
    A = np.zeros((1, n_var)); A[0, :n] = pool.risk_am1bcc.to_numpy(float)
    base.append((A, np.array([-np.inf]), np.array([float(RISK_CAP)])))
    A = np.zeros((1, n_var)); A[0, :n] = (zvals < med).astype(float)
    base.append((A, np.array([float(MIN_LOW)]), np.array([np.inf])))
    for i, j in tc_conflict:
        A = np.zeros((1, n_var)); A[0, i] = A[0, j] = 1.0
        base.append((A, np.array([-np.inf]), np.array([1.0])))
    # 跨度（max_S z − min_S z）严格线性化：独热指派 + 紧约束（禁止单侧下界式的假跨度）
    A = np.zeros((1, n_var)); A[0, I_Y:I_Y + n] = 1.0
    base.append((A, np.array([1.0]), np.array([1.0])))
    A = np.zeros((1, n_var)); A[0, I_W:I_W + n] = 1.0
    base.append((A, np.array([1.0]), np.array([1.0])))
    A = np.zeros((1, n_var)); A[0, I_U] = 1.0; A[0, I_Y:I_Y + n] = -zvals
    base.append((A, np.array([0.0]), np.array([0.0])))
    A = np.zeros((1, n_var)); A[0, I_L] = 1.0; A[0, I_W:I_W + n] = -zvals
    base.append((A, np.array([0.0]), np.array([0.0])))
    for i in range(n):
        A = np.zeros((1, n_var)); A[0, I_Y + i] = 1.0; A[0, i] = -1.0
        base.append((A, np.array([-np.inf]), np.array([0.0])))
        A = np.zeros((1, n_var)); A[0, I_W + i] = 1.0; A[0, i] = -1.0
        base.append((A, np.array([-np.inf]), np.array([0.0])))
        A = np.zeros((1, n_var)); A[0, I_U] = 1.0; A[0, i] = -Mbig
        base.append((A, np.array([zvals[i] - Mbig]), np.array([np.inf])))
        A = np.zeros((1, n_var)); A[0, I_L] = 1.0; A[0, i] = Mbig
        base.append((A, np.array([-np.inf]), np.array([zvals[i] + Mbig])))
    bin_rows: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for members in bin_mem.values():
        if not members:
            continue
        A = np.zeros((1, n_var))
        A[0, [idx_of[c] for c in members]] = 1.0
        A[0, I_T] = -1.0
        bin_rows.append((A, np.array([0.0]), np.array([np.inf])))

    zmin, zmax = float(zvals.min()), float(zvals.max())
    lb = [0.0] * n + [0.0, zmin, zmin] + [0.0] * (2 * n)
    ub = [1.0] * n + [float(len(live_bins)), zmax, zmax] + [1.0] * (2 * n)
    bounds = Bounds(np.array(lb), np.array(ub))
    integrality = np.array([1] * n + [0, 0, 0] + [1] * (2 * n))

    def c_t():
        c = np.zeros(n_var); c[I_T] = -1.0; return c

    def c_span():
        c = np.zeros(n_var); c[I_U] = -1.0; c[I_L] = 1.0; return c

    def c_sum():
        c = np.zeros(n_var); c[:n] = -zvals; return c

    def c_canon():
        c = np.zeros(n_var); c[:n] = -np.arange(n, 0, -1).astype(float); return c

    def solve(rows, c):
        A_ = np.vstack([a for a, _b, _c in rows])
        b_ = np.concatenate([b for _a, b, _c in rows])
        c_ = np.concatenate([cc for _a, _b, cc in rows])
        r = milp(c=c, integrality=integrality, bounds=bounds,
                 constraints=LinearConstraint(A_, b_, c_), options=OPT_OPTS)
        if r.status != 0:
            raise SystemExit(f"[17] MILP 未取得最优解：status={r.status} {r.message}")
        gap = float(getattr(r, "mip_gap", np.nan))
        if not (gap <= MIP_GAP_TOL):
            raise SystemExit(f"[17] MILP 最优性间隙 {gap} > {MIP_GAP_TOL}")
        return r

    def lex_solve():
        fixed: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        stages = [("分层覆盖下限 t", c_t()), ("CNN_VS 跨度", c_span()),
                  ("Σ CNN_VS", c_sum()), ("规范并列打破", c_canon())]
        res, journal = None, []
        for k, (name, c) in enumerate(stages, 1):
            r = solve(list(base) + list(bin_rows) + fixed, c)
            journal.append({"stage": k, "name": name, "objective": float(r.fun),
                            "status": int(r.status),
                            "mip_gap": float(getattr(r, "mip_gap", np.nan)),
                            "nodes": int(getattr(r, "mip_node_count", -1))})
            res = r
            if k == 1:
                v = float(round(r.x[I_T]))
                fixed.append((np.array([[1.0 if i == I_T else 0.0 for i in range(n_var)]]),
                              np.array([v]), np.array([v])))
            elif k == 2:
                row = np.zeros((1, n_var)); row[0, I_U] = 1.0; row[0, I_L] = -1.0
                v = float(r.fun) * -1.0
                fixed.append((row, np.array([v]), np.array([v])))
            elif k == 3:
                row = np.zeros((1, n_var)); row[0, :n] = zvals
                v = float(r.fun) * -1.0
                fixed.append((row, np.array([v]), np.array([v])))
        return [ids[i] for i in range(n) if round(res.x[i]) == 1], journal

    def evaluate(members) -> dict:
        counts = {b: sum(1 for c in members if c in bin_mem[b]) for b in bin_mem}
        live = [v for b, v in counts.items() if bin_mem[b]]
        zv = np.array([float(pool_ix.loc[c, "CNN_VS"]) for c in members])
        return {"t": int(min(live)) if live else 0, "span": float(zv.max() - zv.min()),
                "sum": float(zv.sum()), "counts": counts,
                "n_occupied_bins": int(sum(1 for b in bin_mem if counts[b] > 0 and bin_mem[b]))}

    sel_ids, journal = lex_solve()
    ev_o = evaluate(sel_ids)
    assert ev_o["t"] == round(-journal[0]["objective"]), "阶段1 目标与名册实算 t 不一致"
    assert abs(ev_o["span"] - (-journal[1]["objective"])) < 1e-6, \
        "阶段2 跨度与名册实算不一致 —— 跨度线性化可能退化（假跨度）"
    assert abs(ev_o["sum"] - (-journal[2]["objective"])) < 1e-6, "阶段3 Σ CNN_VS 与名册实算不一致"
    sel_ids_2, journal_2 = lex_solve()
    assert (sel_ids_2 == sel_ids and
            all(abs(a["objective"] - b["objective"]) < 1e-9
                for a, b in zip(journal, journal_2))), "同一输入两次求解名册不一致"
    log.info(f"最优解：t*={ev_o['t']}（占据 {ev_o['n_occupied_bins']}/{len(live_bins)} 个分位格）、"
             f"CNN_VS 跨度 {r4(ev_o['span'])}、Σ CNN_VS {r4(ev_o['sum'])}；"
             f"四阶段 status 全 0、mip_gap=0、重解一致")

    # ------------------------------------------------------------ 6. 对照实现
    heur_sel = greedy_reference(pool, med, N_SLOTS, MIN_LOW)
    assert len(heur_sel) == N_SLOTS, f"对照贪心只选出 {len(heur_sel)} 个"
    heur_ids = [s["cid"] for s in heur_sel]
    ev_h = evaluate(heur_ids)
    assert ev_o["t"] >= ev_h["t"], "对照贪心在分层覆盖上胜过 MILP 最优解（建模错误）"
    if ev_o["t"] == ev_h["t"]:
        assert ev_o["span"] >= ev_h["span"] - 1e-9, "同覆盖下限下对照贪心跨度更大（建模错误）"
        if abs(ev_o["span"] - ev_h["span"]) < 1e-9:
            assert ev_o["sum"] >= ev_h["sum"] - 1e-9, "同为覆盖/跨度最优时对照贪心 Σ 更大（建模错误）"
    n_swaps = len(set(sel_ids) ^ set(heur_ids)) // 2
    log.info(f"对照（v1.6.0 贪心）：t={ev_h['t']}、跨度 {r4(ev_h['span'])}、"
             f"Σ CNN_VS {r4(ev_h['sum'])}；与最优解重叠 {N_SLOTS - n_swaps} 席、互换 {n_swaps} 席")

    # ------------------------------------------------------------ 7. 席位角色 + 反事实重解
    roles: dict[str, str] = {}
    for cid in sel_ids:
        rest = [c for c in sel_ids if c != cid]
        zv_rest = np.array([float(pool_ix.loc[c, "CNN_VS"]) for c in rest])
        live_rest = [sum(1 for c in rest if c in bin_mem[b]) for b in bin_mem if bin_mem[b]]
        if float(zv_rest.max() - zv_rest.min()) < ev_o["span"] - 1e-9:
            roles[cid] = REASON_SPAN
        elif min(live_rest) < ev_o["t"]:
            roles[cid] = REASON_COVER
        else:
            roles[cid] = REASON_FILL
    cf_cover = cf_span = 0
    for cid in sel_ids:
        i = idx_of[cid]
        ban = np.zeros((1, n_var)); ban[0, i] = 1.0
        ban_row = [(ban, np.array([-np.inf]), np.array([0.0]))]
        r1 = solve(list(base) + list(bin_rows) + ban_row, c_t())
        if float(round(r1.x[I_T])) < ev_o["t"]:
            cf_cover += 1
            continue
        t_eq = np.array([[1.0 if k == I_T else 0.0 for k in range(n_var)]])
        r2 = solve(list(base) + list(bin_rows) + ban_row +
                   [(t_eq, np.array([float(ev_o["t"])]), np.array([float(ev_o["t"])]))], c_span())
        if float(-r2.fun) < ev_o["span"] - 1e-9:
            cf_span += 1
    cf_replaceable = len(sel_ids) - cf_cover - cf_span
    log.info(f"席位角色（留一法）：跨度端点 {sum(1 for c in sel_ids if roles[c] == REASON_SPAN)}、"
             f"覆盖下限 {sum(1 for c in sel_ids if roles[c] == REASON_COVER)}、"
             f"质量兜底 {sum(1 for c in sel_ids if roles[c] == REASON_FILL)}；"
             f"反事实重解：覆盖不可替代 {cf_cover}、跨度不可替代 {cf_span}、可替代 {cf_replaceable}")

    # ------------------------------------------------------------ 8. 敏感性对照：来源池内百分位
    pct = {}
    for name, spec in SRC.items():
        pu = pd.read_csv(spec["pool"])
        vals = np.sort(pu.CNN_VS.to_numpy(float))
        for cid, v in zip(pu.cid.astype(str), pu.CNN_VS.to_numpy(float)):
            pct[cid] = 100.0 * float(np.searchsorted(vals, v, side="right")) / len(vals)
    pool_pct = pool.copy()
    pool_pct["pct_src"] = [pct[c] for c in pool_pct.cid]
    assert pool_pct.pct_src.notna().all(), "来源池百分位缺失"
    pvals = pool_pct.pct_src.to_numpy(float)
    pmed = float(pool_pct.pct_src.median())
    pbin_mem, pbin_edges = {}, {}
    for col, name in [("pct_src", "来源池百分位"), ("cLogP", "cLogP"), ("cns_mpo", "CNS MPO")]:
        idx, edges = bins_for(pool_pct[col].to_numpy(float), pool_pct[col].to_numpy(float), N_BINS)
        pbin_edges[name] = [float(e) for e in edges]
        for b in range(N_BINS):
            pbin_mem[(name, b)] = [ids[i] for i in range(n) if idx[i] == b]
    plive = [k for k, v in pbin_mem.items() if v]

    pM = float(pvals.max() - pvals.min()) + 1.0
    pbase: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    A = np.zeros((1, n_var)); A[0, :n] = 1.0
    pbase.append((A, np.array([float(N_SLOTS)]), np.array([float(N_SLOTS)])))
    for c in sorted(pool.chemotype.unique()):
        A = np.zeros((1, n_var))
        A[0, [idx_of[i] for i in ids if pool_ix.loc[i, "chemotype"] == c]] = 1.0
        pbase.append((A, np.array([-np.inf]), np.array([float(CHEM_CAP)])))
    A = np.zeros((1, n_var)); A[0, :n] = pool.risk_am1bcc.to_numpy(float)
    pbase.append((A, np.array([-np.inf]), np.array([float(RISK_CAP)])))
    A = np.zeros((1, n_var)); A[0, :n] = (pvals < pmed).astype(float)
    pbase.append((A, np.array([float(MIN_LOW)]), np.array([np.inf])))
    for i, j in tc_conflict:
        A = np.zeros((1, n_var)); A[0, i] = A[0, j] = 1.0
        pbase.append((A, np.array([-np.inf]), np.array([1.0])))
    A = np.zeros((1, n_var)); A[0, I_Y:I_Y + n] = 1.0
    pbase.append((A, np.array([1.0]), np.array([1.0])))
    A = np.zeros((1, n_var)); A[0, I_W:I_W + n] = 1.0
    pbase.append((A, np.array([1.0]), np.array([1.0])))
    A = np.zeros((1, n_var)); A[0, I_U] = 1.0; A[0, I_Y:I_Y + n] = -pvals
    pbase.append((A, np.array([0.0]), np.array([0.0])))
    A = np.zeros((1, n_var)); A[0, I_L] = 1.0; A[0, I_W:I_W + n] = -pvals
    pbase.append((A, np.array([0.0]), np.array([0.0])))
    for i in range(n):
        A = np.zeros((1, n_var)); A[0, I_Y + i] = 1.0; A[0, i] = -1.0
        pbase.append((A, np.array([-np.inf]), np.array([0.0])))
        A = np.zeros((1, n_var)); A[0, I_W + i] = 1.0; A[0, i] = -1.0
        pbase.append((A, np.array([-np.inf]), np.array([0.0])))
        A = np.zeros((1, n_var)); A[0, I_U] = 1.0; A[0, i] = -pM
        pbase.append((A, np.array([pvals[i] - pM]), np.array([np.inf])))
        A = np.zeros((1, n_var)); A[0, I_L] = 1.0; A[0, i] = pM
        pbase.append((A, np.array([-np.inf]), np.array([pvals[i] + pM])))
    pbin_rows: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for members in pbin_mem.values():
        if not members:
            continue
        A = np.zeros((1, n_var))
        A[0, [idx_of[c] for c in members]] = 1.0
        A[0, I_T] = -1.0
        pbin_rows.append((A, np.array([0.0]), np.array([np.inf])))

    plb = [0.0] * n + [0.0, float(pvals.min()), float(pvals.min())] + [0.0] * (2 * n)
    pub = [1.0] * n + [float(len(plive)), float(pvals.max()), float(pvals.max())] + [1.0] * (2 * n)
    pbounds = Bounds(np.array(plb), np.array(pub))

    def psolve(rows, c):
        A_ = np.vstack([a for a, _b, _c in rows])
        b_ = np.concatenate([b for _a, b, _c in rows])
        c_ = np.concatenate([cc for _a, _b, cc in rows])
        r = milp(c=c, integrality=integrality, bounds=pbounds,
                 constraints=LinearConstraint(A_, b_, c_), options=OPT_OPTS)
        if r.status != 0:
            raise SystemExit(f"[17] 敏感性口径 MILP 未取得最优解：status={r.status}")
        return r

    def p_c_sum():
        c = np.zeros(n_var); c[:n] = -pvals; return c

    pfixed: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    pres, pjournal = None, []
    for k, (name, c) in enumerate([("分层覆盖下限 t", c_t()), ("来源池百分位跨度", c_span()),
                                   ("Σ 来源池百分位", p_c_sum()),
                                   ("规范并列打破", c_canon())], 1):
        r = psolve(list(pbase) + list(pbin_rows) + pfixed, c)
        pjournal.append({"stage": k, "name": name, "objective": float(r.fun),
                         "status": int(r.status),
                         "mip_gap": float(getattr(r, "mip_gap", np.nan))})
        pres = r
        if k == 1:
            v = float(round(r.x[I_T]))
            pfixed.append((np.array([[1.0 if i == I_T else 0.0 for i in range(n_var)]]),
                           np.array([v]), np.array([v])))
        elif k == 2:
            row = np.zeros((1, n_var)); row[0, I_U] = 1.0; row[0, I_L] = -1.0
            v = float(r.fun) * -1.0
            pfixed.append((row, np.array([v]), np.array([v])))
        elif k == 3:
            row = np.zeros((1, n_var)); row[0, :n] = pvals
            v = float(r.fun) * -1.0
            pfixed.append((row, np.array([v]), np.array([v])))
    sens_ids = [ids[i] for i in range(n) if round(pres.x[i]) == 1]
    overlap = len(set(sens_ids) & set(sel_ids))
    log.info(f"敏感性对照（来源池内百分位口径）：t={round(-pjournal[0]['objective'])}、"
             f"与主名册重叠 {overlap}/{N_SLOTS} 席、互换 {N_SLOTS - overlap} 席")

    def src_counts(members) -> dict:
        return pd.Series([pool_ix.loc[c, "source_branch"] for c in members]).value_counts().to_dict()

    sens = {
        "axis": "CNN_VS → 该分子在**自己来源池**内的 CNN_VS 百分位（目录池 4,990 / 生成池 4,979）",
        "why": "SOP_5 的池是目录池按 CNN_VS 截断出的 top-4,990（下限 6.554），ML 池是生成库"
               "全范围（0.963–8.536）；原始 CNN_VS 跨来源不可直接比，本对照把尺度差异中性化",
        "n_bins": len(plive), "bin_edges": {k: [round(e, 4) for e in v] for k, v in pbin_edges.items()},
        "median_pct": r4(pmed),
        "journal": pjournal,
        "roster": sens_ids,
        "roster_ordered": sorted(sens_ids, key=lambda c: (-pool_pct.set_index("cid").loc[c, "pct_src"], c)),
        "source_counts": src_counts(sens_ids),
        "overlap_with_main": overlap,
        "swapped": {"out_of_sensitivity": [c for c in sel_ids if c not in set(sens_ids)],
                    "in_of_sensitivity": [c for c in sens_ids if c not in set(sel_ids)]},
        "note": "本对照只作留痕，**不改主名册**。重叠度低说明主名册的来源构成确实受原始"
                "CNN_VS 尺度影响；重叠度高说明结论对尺度不敏感。",
    }

    # ------------------------------------------------------------ 9. 可合成性分诊（只标注，不参与选席）
    synth = []
    for spec in SRC.values():
        s = pd.read_csv(spec["synth"])
        synth.append(s[["cid", "sa_score", "sa_band", "complex_flags", "triage_band",
                        "catalog_source"]])
    synth = pd.concat(synth, ignore_index=True).drop_duplicates("cid")
    synth_ix = synth.set_index("cid")

    # ------------------------------------------------------------ 10. 名册表
    batch = cand[cand.cid.isin(set(sel_ids))].copy()
    batch = batch.sort_values(["CNN_VS", "cid"], ascending=[False, True], kind="mergesort")
    batch = batch[[c for c in KEEP_COLS if c in batch.columns]].copy()
    _ann = ["cid", "hERG_cross_ens", "m2_pass_ens", "herg_cross_model_flag",
            "source_branch", "md_eligibility_note"]
    batch = batch.merge(cand[[c for c in _ann if c in cand.columns]], on="cid", how="left")
    batch["candidate_type"] = "candidate"
    batch["control_role"] = ""
    batch["selection_reason"] = batch.cid.map(roles)
    for col, name in AXES:
        _e = np.array(bin_edges[name])
        batch[f"bin_{name}"] = ["Q" + str(int(np.searchsorted(_e, v, side="left")) + 1)
                                for v in batch[col].to_numpy(float)]
    batch["coverage_bins"] = [";".join(f"{nm}={batch.loc[i, 'bin_' + nm]}" for _c, nm in AXES)
                              for i in batch.index]
    batch = batch.merge(synth_ix[["sa_score", "sa_band", "complex_flags", "triage_band",
                                  "catalog_source"]], on="cid", how="left")
    batch["source_short"] = batch.source_branch.map({k: v["short"] for k, v in SRC.items()})

    ctrl_out = ctrl.copy().merge(
        synth_ix[["sa_score", "sa_band", "complex_flags", "triage_band",
                  "catalog_source"]].reset_index(), on="cid", how="left")
    for c in batch.columns:
        if c not in ctrl_out.columns:
            ctrl_out[c] = np.nan
    ctrl_out["candidate_type"] = "negative_control"
    ctrl_out["selection_reason"] = REASON_CTRL
    ctrl_out["source_branch"] = "两分支共用"
    ctrl_out["source_short"] = "对照"
    roster = pd.concat([batch, ctrl_out[batch.columns]], ignore_index=True)
    roster["batch_order"] = range(1, len(roster) + 1)
    roster["md_card_id"] = [f"A2AR-MD-{i:02d}" for i in roster.batch_order]
    assert len(roster) == N_SLOTS + N_CTRL, f"名册 {len(roster)} 行 != {N_SLOTS + N_CTRL}"

    # ------------------------------------------------------------ 11. 自检（硬断言）
    bc = roster[roster.candidate_type == "candidate"]
    chem_cnt = bc.chemotype.value_counts().to_dict()
    assert max(chem_cnt.values()) <= CHEM_CAP, f"化学型超限：{chem_cnt}"
    assert int(bc.risk_am1bcc.sum()) <= RISK_CAP, f"风险档超限：{int(bc.risk_am1bcc.sum())}"
    assert bc.scaffold_generic.nunique() == len(bc), "名册内通用骨架重复"
    assert bc.scaffold.nunique() == len(bc), "名册内 Bemis-Murcko 骨架重复"
    n_low_final = int((bc.CNN_VS < med).sum())
    assert n_low_final >= MIN_LOW, f"下半区候选 {n_low_final} < {MIN_LOW}"
    cids = bc.cid.tolist()
    tc_max, tc_pair = 0.0, None
    for i in range(len(cids)):
        for j in range(i + 1, len(cids)):
            t = DataStructs.TanimotoSimilarity(fps_all[cids[i]], fps_all[cids[j]])
            if t > tc_max:
                tc_max, tc_pair = t, (cids[i], cids[j])
    assert tc_max < TC_CAP, f"名册内成对 Tc {tc_max} ≥ {TC_CAP}"
    assert bc.cns_mpo.notna().all(), "CNS MPO 缺失"
    assert bc.hERG_cross_ens.notna().all(), "第二模型 hERG 列缺失"
    assert bc.sa_score.notna().all(), "可合成性读数缺失（17 个生成分子应已在 run05 派生层算过）"
    assert set(excl_cid := {e["cid"] for e in excl}).isdisjoint(set(cids)), "被剔除分子混入名册"
    ev_csv = evaluate(cids)
    assert ev_csv["t"] == ev_o["t"] and abs(ev_csv["span"] - ev_o["span"]) < 1e-9, \
        "名册表与求解结果不一致"

    write_csv(OUT / "tables/final_md_20.csv", roster)
    write_csv(OUT / "tables/final_md_excluded.csv",
              pd.DataFrame(excl + [{"cid": c, "source_branch": pool_ix.loc[c, "source_branch"],
                                    "CNN_VS": r4(float(pool_ix.loc[c, "CNN_VS"])),
                                    "reason": "资格池内未入 20 席（名额限制）",
                                    "detail": f"t 分层覆盖 / 跨度 / Σ CNN_VS 字典序下不占席"
                                              f"（角色：{roles.get(c, '未入席')}）"}
                                   for c in pool.cid if c not in set(cids)]))

    # ------------------------------------------------------------ 12. 证书
    comp = {"n_selected": len(cids), "chemotype_counts": chem_cnt,
            "risk_count": int(bc.risk_am1bcc.sum()),
            "below_median_count": n_low_final, "median_reference": r4(med),
            "max_pairwise_tc": r6(tc_max), "max_pairwise_tc_pair": tc_pair,
            "min_bin_count": ev_o["t"], "n_occupied_bins": ev_o["n_occupied_bins"],
            "sigma_CNN_VS": r6(float(bc.CNN_VS.std(ddof=1))),
            "pass": {"n_selected": len(cids) == N_SLOTS, "chem_cap": max(chem_cnt.values()) <= CHEM_CAP,
                     "risk_cap": int(bc.risk_am1bcc.sum()) <= RISK_CAP,
                     "below_median_min": n_low_final >= MIN_LOW, "pairwise_tc": tc_max < TC_CAP,
                     "scaffold_unique": bc.scaffold_generic.nunique() == len(bc)}}
    comp_h = {"t": ev_h["t"], "span": r4(ev_h["span"]), "sum": r4(ev_h["sum"]),
              "source_counts": src_counts(heur_ids)}
    role_counts = {REASON_COVER: sum(1 for c in sel_ids if roles[c] == REASON_COVER),
                   REASON_SPAN: sum(1 for c in sel_ids if roles[c] == REASON_SPAN),
                   REASON_FILL: sum(1 for c in sel_ids if roles[c] == REASON_FILL)}

    rep = {
        "step": "17_final_md_merge",
        "layer": "派生分析层（零侵入：不改门、不改名单、不改规则版本；不在 12_determinism 的 INCLUDE_DIRS 内）",
        "rule_version": rule["rule_version"],
        "purpose": "把 SOP_5（目录）与 SOP_ML_run05（ML 生成）两条链的 MD 一级批次候选并成"
                   " 35 个池，在合并池上重解名额 20 的 MD 批次",
        "inputs": {name: {"batch": _rel(spec["batch"]), "batch_sha256": sha256_file(spec["batch"]),
                          "pool": _rel(spec["pool"]), "pool_sha256": sha256_file(spec["pool"]),
                          "synth_access": _rel(spec["synth"]),
                          "synth_access_sha256": sha256_file(spec["synth"]),
                          "n_batch_rows": int(len(pd.read_csv(spec["batch"])))}
                   for name, spec in SRC.items()},
        "frozen_rule_sha256": sha256_file(RULE),
        "merge": {"n_candidates": int(len(cand)), "n_controls": int(len(ctrl)),
                  "n_controls_dedupe": "两个共用阴性对照（C09462894 / C60359994）只留一份",
                  "structural_overlap_between_branches": 0,
                  "overlap_note": "两份名册按 InChIKey 结构零重叠；唯二相同的两个就是共用对照",
                  "source_counts": {k: int(v) for k, v in
                                    cand.source_branch.value_counts().items()}},
        "eligibility": {"hERG_cross_model_policy": "retain_and_annotate（保留入池 + 逐分子标注）",
                        "n_hERG_cross_model_fail": int(len(fail)),
                        "n_hERG_cross_model_fail_ids": [str(c) for c in fail.cid],
                        "dedupe_generic_scaffold": "scaffold_generic ≤1，重复时保留 CNN_VS 高者",
                        "n_deduped": int(len(dup)), "n_pool": int(len(pool)),
                        "pool_source_counts": src_counts(ids)},
        "parameters": {
            "n_slots": N_SLOTS, "n_controls": N_CTRL, "batch_rows": N_SLOTS + N_CTRL,
            "chem_cap": CHEM_CAP, "risk_cap": RISK_CAP, "pairwise_tc_max": TC_CAP,
            "scaffold_cap": 1, "min_below_median_candidates": MIN_LOW,
            "median_reference": f"合并候选池（{len(cand)}）的 CNN_VS 中位数 = {r4(med)}",
            "median_reference_note": "原链参照 48 候选的 7.414；本层按规则「用本层候选名单"
                                     "自身中位」的行为取合并池中位",
            "chem_cap_note": "冻结登记值为 8（对应 20 行批次）。规则自带的按席位缩放式"
                             "round(8N/20) 在 N=22 下给出 9，本层取**更严**的 8，不新增登记值",
            "risk_cap_note": "round(3N/20)：N=20 与 N=22 均得 3，两种取法一致",
            "objective": "frozen_rule.step6_md_batch.objective（字典序四阶段，逐字沿用）",
        },
        "bins": {"axes": [n for _c, n in AXES], "n_bins_per_axis": N_BINS,
                 "n_live_bins": len(live_bins), "empty_bins": empty_bins,
                 "bin_edges": {k: [round(e, 4) for e in v] for k, v in bin_edges.items()},
                 "bin_members": {f"{k[0]}|Q{k[1] + 1}": v for k, v in bin_mem.items()},
                 "bin_reference": "eligible_pool 自分布（25/50/75 分位，线性插值）"},
        "optimality": {"journal": journal, "t_star": ev_o["t"], "span": r4(ev_o["span"]),
                       "sum_CNN_VS": r4(ev_o["sum"]), "n_occupied_bins": ev_o["n_occupied_bins"],
                       "bin_counts": {f"{k[0]}|Q{k[1] + 1}": v
                                      for k, v in ev_o["counts"].items() if bin_mem[k]},
                       "resolved_consistent": True, "all_status_optimal": True,
                       "mip_gap": 0.0, "compliance": comp},
        "roles": {"leave_one_out": role_counts,
                  "counterfactual": {"cover_irreplaceable": cf_cover, "span_irreplaceable": cf_span,
                                     "replaceable": cf_replaceable},
                  "definition": "留一法说明该席在已定名册里承担什么；反事实重解说明该席是否不可替代"},
        "heuristic_reference": {"kind": "v1.6.0 贪心填充（化学型锚点 → 下半区补齐 → 高分填充）",
                                "t": ev_h["t"], "span": r4(ev_h["span"]), "sum_CNN_VS": r4(ev_h["sum"]),
                                "n_swaps": n_swaps, "compliance": comp_h,
                                "swapped_out": [c for c in heur_ids if c not in set(sel_ids)],
                                "swapped_in": [c for c in sel_ids if c not in set(heur_ids)]},
        "sensitivity": sens,
        "roster": {"n": len(roster), "order": [{"md_card_id": r.md_card_id, "cid": r.cid,
                                                "source_short": r.source_short,
                                                "CNN_VS": r4(float(r.CNN_VS)),
                                                "selection_reason": r.selection_reason}
                                               for r in roster.itertuples()],
                   "source_counts": src_counts(cids),
                   "md_cost": {"level1_10ns_cardday_per_molecule":
                               float(mb["md_cost"]["level1_10ns_cardday_per_molecule"]),
                               "level2_150ns_cardday_per_molecule":
                               float(mb["md_cost"]["level2_150ns_cardday_per_molecule"]),
                               "level1_total_cardday": r4(N_SLOTS * float(
                                   mb["md_cost"]["level1_10ns_cardday_per_molecule"]))}},
        "synth_access": {"policy": "只作标注列，不参与筛选或加权",
                         "triage_counts": bc.triage_band.value_counts().to_dict(),
                         "note": "两分支的 09_synth_access 口径与 SOP_5 步骤 16 逐条一致；"
                                 "绿带是「启发式未报警」，不是「确认好合成」"},
        "caveats": [
            "跨来源 CNN_VS 不可直接比：SOP_5 池是目录池按 CNN_VS 截断出的 top-4,990（下限 6.554），"
            "ML 池是生成库全范围（0.963–8.536）——见 sensitivity 对照",
            "名册排序是「合并池内按本层目标函数最优」，不是「这 20 个比落选的 13 个更好」",
            "合并后的 CNN_VS 跨度 3.5015 远大于 SOP_5 单链登记值 0.7337，差量全部来自生成侧池"
            "未被截断（最小 CNN_VS 4.659 的 SCALE-103966 是跨度端点席位）；这不是「更好的批次」，"
            "是合并把两条链的量程差异带进来了",
            "hERG 只是单一模型的分类线，只能说「在 ADMET-AI 口径下过门」；第二模型不过门者"
            "已保留并标注「模型依赖」，本层不产生安全性结论",
            "ML 生成分子无实验数据；目录来源不等于现货可得",
            "2 个共用阴性对照是 DUD-E 诱饵，只作参照，不占候选名额",
            "本层是机时分配层，不产生有效性结论；MD 后若需换人须改规则并重跑整链",
        ],
        "toolchain": {"scipy": scipy.__version__, "numpy": np.__version__, "pandas": pd.__version__,
                      "timestamp_policy": "不写时间戳，保证同一份输入两次运行输出逐字节一致"},
    }
    write_json(OUT / "reports/final_md.json", rep)

    # ------------------------------------------------------------ 13. Markdown 小结
    src5, srcm = list(SRC)[0], list(SRC)[1]
    sc = src_counts(cids)
    lines = []
    A = lines.append
    A("# SOP-05 · 最终 MD 候选 20 席（步骤 17 · 合并两分支批次 · 派生分析层）\n")
    A(f"**结论一句话**：把 SOP_5（目录化合物）与 SOP_ML_run05（ML 生成分子）的 MD 一级批次"
      f"候选并成 {len(cand)} 个池，通用骨架去重后资格池 {len(pool)} 个，在其中重解出 "
      f"**{N_SLOTS} 席**最终 MD 候选（目录 {sc.get(src5, 0)} 席 + 生成 {sc.get(srcm, 0)} 席），"
      f"加 2 个共用阴性对照共 {len(roster)} 行。字典序最优：分层覆盖下限 t\\*={ev_o['t']}"
      f"（占据 {ev_o['n_occupied_bins']}/{len(live_bins)} 个分位格）、CNN_VS 跨度 "
      f"{r4(ev_o['span'])}、Σ CNN_VS {r4(ev_o['sum'])}。\n")
    A("## 1. 合并口径\n")
    A(f"| 项 | 值 |\n| --- | --- |\n"
      f"| 合并候选池 | {len(cand)}（目录 {int((cand.source_branch == src5).sum())} + "
      f"生成 {int((cand.source_branch == srcm).sum())}） |\n"
      f"| 两分支结构重叠 | 0（按 InChIKey；唯二相同的是共用阴性对照） |\n"
      f"| 通用骨架去重剔除 | {len(dup)} 个（均为生成侧） |\n"
      f"| 资格池 | {len(pool)}（目录 {int((pool.source_branch == src5).sum())} + "
      f"生成 {int((pool.source_branch == srcm).sum())}） |\n"
      f"| 第二 hERG 模型不过门（保留并标注） | {len(fail)} 个 |\n"
      f"| 成对 ECFP4 Tc 最大 | {tc_max_pool:.4f}（< {TC_CAP}，本层无 Tc 冲突） |\n")
    A("## 2. 选席目标与约束（沿用 frozen_rule v1.7.0 step6_md_batch）\n")
    A("字典序四阶段：① max 分层覆盖下限 t（3 轴 × 4 层 = "
      f"{n_bins} 格，参照资格池自分布）→ ② 在 t=t\\* 内 max CNN_VS 跨度 → "
      "③ 仍并列时 max Σ CNN_VS → ④ 规范化并列打破。每阶段 status=0、mip_gap=0，"
      "目标值等于名册实算值，且同输入重解一致。\n")
    A(f"约束：|S| = {N_SLOTS}、成对 ECFP4 Tc < {TC_CAP}、骨架 ≤1、化学型 ≤{CHEM_CAP}、"
      f"AM1-BCC 风险档 ≤{RISK_CAP}、≥ {MIN_LOW} 席来自 CNN_VS 下半区"
      f"（中位参照 = 合并候选池 {len(cand)} 的中位 {r4(med)}）。\n")
    A(f"配额口径说明：化学型帽取冻结登记值 {CHEM_CAP}；规则自带的按席位缩放式 "
      f"round(8N/20) 在 N={N_SLOTS + N_CTRL} 下会给出 {round_half(8 * (N_SLOTS + N_CTRL) / 20)}，"
      f"本层取更严的一侧。风险档两种取法都得 {RISK_CAP}。\n")
    A("## 3. 名册（按 CNN_VS 降序，对照固定占末两位）\n")
    A("| 卡号 | cid | 来源 | CNN_VS | 分位覆盖 | 化学型 | 席位角色 | 可合成性 |\n"
      "| --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in roster.itertuples():
        A(f"| {r.md_card_id} | {r.cid} | {r.source_short} | "
          f"{'' if pd.isna(r.CNN_VS) else f'{float(r.CNN_VS):.3f}'} | "
          f"{'' if not isinstance(r.coverage_bins, str) else r.coverage_bins} | "
          f"{'' if not isinstance(r.chemotype, str) else r.chemotype} | "
          f"{'' if not isinstance(r.selection_reason, str) else r.selection_reason} | "
          f"{'' if pd.isna(r.sa_score) else f'{float(r.sa_score):.3f} '}"
          f"{'' if not isinstance(r.triage_band, str) else r.triage_band} |")
    A("")
    A("## 4. 席位角色（留一法 + 反事实重解）\n")
    A(f"- 留一法：覆盖下限 {role_counts[REASON_COVER]}、跨度端点 {role_counts[REASON_SPAN]}、"
      f"质量兜底 {role_counts[REASON_FILL]}。\n"
      f"- 反事实重解（禁掉该席后仍可达标 ⇒ 可替代）：覆盖不可替代 {cf_cover}、"
      f"跨度不可替代 {cf_span}、可替代 {cf_replaceable}。\n"
      f"- 对照实现（v1.6.0 贪心）：t={ev_h['t']}、跨度 {r4(ev_h['span'])}、"
      f"Σ CNN_VS {r4(ev_h['sum'])}；与最优解重叠 {N_SLOTS - n_swaps} 席、互换 {n_swaps} 席。\n")
    A("## 5. 跨来源 CNN_VS 尺度敏感性对照（只作留痕，不改主名册）\n")
    A(f"把 CNN_VS 换成「该分子在**自己来源池**内的百分位」重解一次：与主名册重叠 "
      f"**{overlap}/{N_SLOTS} 席**、互换 {N_SLOTS - overlap} 席"
      f"（敏感性名册来源构成：{sens['source_counts']}）。\n")
    A("这一对照是必需的：SOP_5 的池是目录池按 CNN_VS 截断出的 top-4,990（下限 6.554），"
      "而 ML 池是生成库全范围（0.963–8.536）。阶段 3 的 Σ CNN_VS 会系统性偏向目录来源，"
      "这是选池口径差异而非化学优劣。重叠度越低，说明主名册的来源构成越受这个尺度影响。\n")
    A("## 6. 必须随结果转述的边界\n")
    for c in rep["caveats"]:
        A(f"- {c}")
    A("")
    A("## 7. 复现\n")
    A("```bash\ncd <本包根目录>\npython 10_final_md/scripts/17_final_md_merge.py\n"
      "python 10_final_md/scripts/18_final_md_figures.py\n```\n")
    A(f"产物：`tables/final_md_20.csv`（{len(roster)} 行）、`tables/final_md_excluded.csv`、"
      f"`reports/final_md.json`、`figures/fig_final_md_cards.png`、"
      f"`figures/fig_final_md_batch.png`。\n")
    write_text(OUT / "reports/FINAL_MD_20.md", "\n".join(lines))

    log.info(f"名册来源构成：{sc}；可合成性分诊：{bc.triage_band.value_counts().to_dict()}")
    log.info(f"写出 tables/final_md_20.csv（{len(roster)} 行）、tables/final_md_excluded.csv"
             f"（{len(excl) + (len(pool) - len(cids))} 行）、reports/final_md.json、"
             f"reports/FINAL_MD_20.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
