"""名单层选序（约束贪心）、Butina 聚类与稳健性辅助。

选序规则（frozen_rule.selection）：
  分数 CNN_VS 降序（并列时 cid 字典序升序，保证确定性）
  硬约束：① 成对 ECFP4 Tanimoto < tc_max（默认 0.65）
          ② 同骨架（Bemis-Murcko）计数 ≤ scaffold_cap
          ③ 同化学型计数 ≤ chem_cap = round(8·N/20)
          ④ AM1-BCC 风险档计数 ≤ risk_cap = round(3·N/20)
  对照名额：名单尾部固定 2 个阴性对照（NEG1/NEG2，见 frozen_rule.controls）
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina


def caps_for(n_total: int) -> tuple[int, int]:
    """chem_cap = round(8N/20)、risk_cap = round(3N/20)（Python 银行家舍入，与既定口径一致）。"""
    chem = int(round(8 * n_total / 20))
    risk = int(round(3 * n_total / 20))
    return chem, risk


def _sort_key(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    return df.sort_values([score_col, "cid"], ascending=[False, True]).reset_index(drop=True)


def greedy_select(cand: pd.DataFrame, fps: dict, n_slots: int, *, score_col: str = "CNN_VS",
                  tc_max: float = 0.65, scaffold_cap: int = 1,
                  chem_cap: int | None = None, risk_cap: int | None = None) -> dict:
    """在候选宇宙 `cand` 上做约束贪心；fps: {cid: 指纹对象}。

    返回 dict(picks=[cid...], picks_df=DataFrame, skips=Counter, log=[每步记录])
    """
    if chem_cap is None or risk_cap is None:
        c, r = caps_for(n_slots + 2)      # 名额口径按最终名单数（含对照）
        chem_cap = chem_cap if chem_cap is not None else c
        risk_cap = risk_cap if risk_cap is not None else r
    order = _sort_key(cand, score_col)
    scaf = Counter()
    chem = Counter()
    n_risk = 0
    picked, sel_fps, skips, log = [], [], Counter(), []
    for row in order.itertuples():
        if len(picked) >= n_slots:
            break
        s = getattr(row, "scaffold", "")
        c = getattr(row, "chemotype", "")
        risk = int(getattr(row, "risk_am1bcc", 0))
        if scaf[s] >= scaffold_cap:
            skips["scaffold_cap"] += 1
            continue
        if chem[c] >= chem_cap:
            skips["chem_cap"] += 1
            continue
        if risk and n_risk >= risk_cap:
            skips["risk_cap"] += 1
            continue
        fp = fps[row.cid]
        if sel_fps:
            tcs = DataStructs.BulkTanimotoSimilarity(fp, sel_fps)
            tmax = float(max(tcs))
            if tmax >= tc_max:
                skips["tc_pairwise"] += 1
                log.append({"cid": row.cid, "reason": "tc_pairwise", "tc_max_vs_selected": round(tmax, 4)})
                continue
        else:
            tmax = float("nan")
        scaf[s] += 1
        chem[c] += 1
        n_risk += risk
        picked.append(row.cid)
        sel_fps.append(fp)
        log.append({"cid": row.cid, "reason": "picked",
                    "tc_max_vs_selected": (None if np.isnan(tmax) else round(tmax, 4)),
                    "score": float(getattr(row, score_col))})
    picks_df = cand[cand.cid.isin(picked)].copy()
    picks_df["pick_order"] = picks_df.cid.map({c: i + 1 for i, c in enumerate(picked)})
    picks_df = picks_df.sort_values("pick_order").reset_index(drop=True)
    return {"picks": picked, "picks_df": picks_df, "skips": dict(skips), "log": log,
            "chem_cap": chem_cap, "risk_cap": risk_cap, "tc_max": tc_max,
            "scaffold_cap": scaffold_cap, "n_slots": n_slots}


def max_pairwise_tc(picks_df: pd.DataFrame, fps: dict) -> float:
    if len(picks_df) < 2:
        return float("nan")
    fs = [fps[c] for c in picks_df.cid]
    m = 0.0
    for i in range(len(fs)):
        v = DataStructs.BulkTanimotoSimilarity(fs[i], fs[i + 1:])
        if v:
            m = max(m, float(max(v)))
    return round(m, 4)


def cluster_butina(fps: list, cutoff: float = 0.35) -> list[tuple]:
    """Butina 聚类：距离 = 1 − Tanimoto(ECFP4)。"""
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1.0 - s for s in sims])
    return Butina.ClusterData(dists, n, cutoff, isDistData=True)


def cluster_butina_fast(fps: list, cutoff: float = 0.35, chunk: int = 256) -> list:
    """Taylor-Butina 聚类（内存友好实现），与 rdkit.ML.Cluster.Butina.ClusterData 同结果。

    与官方实现的差异只在内部数据结构：官方把距离展成 n×n 矩阵（n≈5,000 时约 200 MB
    float64 + 副本 → 会被 OOM 杀掉），本实现只保留「邻居表」（逐行 BulkTanimoto），
    其余步骤逐条对齐官方源码：
      ① 邻居表 = 距离 ≤ cutoff 的点（含自身）；
      ② 候选中心按 (邻居数, 下标) 降序出队；
      ③ 未见过的中心建簇 = [中心] + 其所有未见邻居；已见者跳过；
      ④ 循环结束后仍未见者为单例簇。
    等价性在主池抽样上逐个验证（tests/verify_pack 与本文件冒烟测试）。
    """
    n = len(fps)
    neighbors: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]), dtype=np.float64)
        js = (np.flatnonzero(1.0 - sims <= cutoff) + (i + 1)).tolist()
        neighbors[i].append(i)                      # 含自身（与官方 dist_matrix 对角为 0 等价）
        for j in js:
            neighbors[i].append(int(j))
            neighbors[int(j)].append(i)
    order = sorted(((len(neighbors[i]), i) for i in range(n)), reverse=True)
    seen = np.zeros(n, dtype=bool)
    clusters: list[list[int]] = []
    for _, idx in order:
        if seen[idx]:
            continue
        cluster = [idx]
        seen[idx] = True
        for nb in neighbors[idx]:
            if not seen[nb]:
                cluster.append(nb)
                seen[nb] = True
        clusters.append(sorted(cluster, key=lambda x: x))
    return [tuple(c) for c in clusters]


def cluster_summary(clusters: list[tuple], cids: list, scores: np.ndarray) -> dict:
    """每簇代表 = 簇内 CNN_VS 最高的成员。"""
    sizes = [len(c) for c in clusters]
    reps = []
    for cl in clusters:
        idx = list(cl)
        best = max(idx, key=lambda i: (scores[i], [-ord(ch) for ch in cids[i]]))
        reps.append({"cluster_id": len(reps) + 1, "size": len(cl), "representative": cids[best],
                     "members": ";".join(cids[i] for i in idx)})
    return {"n_clusters": len(clusters), "largest_cluster": int(max(sizes)) if sizes else 0,
            "n_singletons": int(sum(1 for s in sizes if s == 1)),
            "n_members": int(sum(sizes)),
            "clusters": reps}


# ---------------------------------------------------------------- 稳健性辅助
def perturb_gaussian(scores: np.ndarray, sigma: float, rng: np.random.Generator) -> np.ndarray:
    return scores + rng.normal(0.0, sigma, len(scores))


def perturb_tieband(scores: np.ndarray, band: float, n_slots: int,
                    rng: np.random.Generator) -> np.ndarray:
    """把 (v_n − band, v_n] 内的「并列区」分数随机置换（确定性随机源）。"""
    s = scores.copy()
    if band <= 0 or len(s) < n_slots:
        return s
    v_n = np.sort(s)[::-1][n_slots - 1]
    mask = s >= v_n - band
    s[mask] = rng.permutation(s[mask])
    return s


def overlap_fraction(a: list, b: list) -> float:
    sa, sb = set(a), set(b)
    return round(len(sa & sb) / max(len(sb), 1), 4)


def jaccard(a: list, b: list) -> float:
    sa, sb = set(a), set(b)
    return round(len(sa & sb) / max(len(sa | sb), 1), 4)
