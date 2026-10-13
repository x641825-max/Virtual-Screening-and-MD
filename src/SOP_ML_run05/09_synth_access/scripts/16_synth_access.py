#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 16（派生分析层 · 零侵入）：run05 名单可合成性评估。

定位
    这是主链之外的分析层，与 08_visuals 同类：**不改门、不改名单、不改规则版本**，
    产物只落 09_synth_access/。目录不在 12_determinism 的 INCLUDE_DIRS 内，也不重生成
    PACK_MANIFEST.sha256 / 交付 MANIFEST.sha256，因此不影响主链确定性基线。

评估对象
    主口径名单（04_selection/tables/list_main_QED0.50_anchor4.csv，run05 实测 34 行 =
    32 候选 + 2 阴性对照）。对照是 DUD-E 诱饵（PubChem CID，已有实体的已知化合物），
    其 SA 分只作参照，不进分诊带统计。行数、MD 批次席位数均由输入实测得出，不写死。

三套互补口径（不合并成一个「可合成性真值」）
    ① SAScore  : RDKit Contrib/SA_Score 的 sascorer（Ertl & Schuffenhauer 2009）
                  —— 复杂度–片段贡献启发式分，1（易）→ 10（难）；官方未给分档。
    ② 复杂旗标 : 8 条确定性结构规则（见 FLAG_DEFS），其中 4 条记作硬旗标。
    ③ 来源标注 : cid 前缀 → ZINC / MCULE（可按需合成目录）、PubChem CID（实体化合物）
                  或生成模型输出（SCALE- 前缀，run05 新增分支：**无目录来源**）。
                  **来源标注不等于现货可得**（本层无法离线核验库存与报价）。

分诊带（本层自定，非论文判据；定义随报告落盘；阈值与 SOP_5 步骤 16 完全一致）
    橙 优先人工评估 : SA > 6.0 或 命中任一硬旗标 或 软旗标 ≥ 3
    黄 需确认路线   : 4.0 < SA ≤ 6.0 或 软旗标 1–2
    绿 常规合成     : SA ≤ 4.0 且 无旗标

输出（09_synth_access/）
    tables/synth_access_list.csv      逐分子读数（含全部中间量）
    reports/synth_access.json         口径定义 + 统计 + 输入/工具哈希（无时间戳，可复现）
    reports/SYNTH_ACCESS.md           Markdown 小结
    figures/fig_synth_access.png/svg  三联图：SA 分布 / 旗标频次 / CNN_VS–SA 权衡

用法：python3 09_synth_access/scripts/16_synth_access.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
ROOT = _p
sys.path.insert(0, str(ROOT))

from rdkit import Chem, RDLogger, __version__ as RDKIT_VERSION  # noqa: E402
from rdkit.Chem import RDConfig, rdMolDescriptors  # noqa: E402

sys.path.insert(0, os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

from vendor import plots  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402  （后端已由 vendor.plots 设为 Agg）

RDLogger.DisableLog("rdApp.*")

OUT = ROOT / "09_synth_access"
LIST_CSV = ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv"
BATCH_CSV = ROOT / "07_md_batch/tables/md_batch20.csv"

SA_BANDS = [(0.0, 4.0, "易（SA ≤ 4.0）"), (4.0, 6.0, "中（4.0 < SA ≤ 6.0）"),
            (6.0, 10.01, "难（SA > 6.0）")]
BAND_COLOR = {"易（SA ≤ 4.0）": plots.GREEN, "中（4.0 < SA ≤ 6.0）": plots.AMBER,
              "难（SA > 6.0）": plots.RED}
TRIAGE_COLOR = {"绿 常规合成": plots.GREEN, "黄 需确认路线": plots.AMBER, "橙 优先人工评估": plots.RED}
ALLOWED_ELEMENTS = {"C", "N", "O", "S", "F", "Cl", "Br", "I", "P"}

# 8 条难点线；hard=True 者一旦命中直接进橙带（与 SOP_5 步骤 16 同一套定义）
FLAG_DEFS = {
    "stereo_ge4": ("立体中心 ≥ 4（含未指定手性中心）", True),
    "bridgehead": ("存在桥环原子（CalcNumBridgeheadAtoms > 0）", True),
    "macrocycle": ("存在环大小 ≥ 12 的大环", True),
    "rare_element": ("含 C/N/O/S/F/Cl/Br/I/P 之外的元素", True),
    "stereo_2_3": ("立体中心 2–3（手性合成/拆分负担）", False),
    "fused_tricyclic": ("最大稠合环系 ≥ 3 环（三环及以上稠合）", False),
    "spiro": ("螺环原子 ≥ 2", False),
    "ring_rich": ("环总数 ≥ 5", False),
}
HARD_FLAGS = tuple(k for k, (_, h) in FLAG_DEFS.items() if h)


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _f(x) -> float:
    """NaN 安全的 float（run05 名单里阴性对照的 CNN_VS/MW/QED 等列为空）。"""
    return float(x) if pd.notna(x) else float("nan")


def sa_band(v: float) -> str:
    for lo, hi, name in SA_BANDS:
        if v <= hi and (v > lo or lo == 0.0):
            return name
    return SA_BANDS[-1][2]


def triage(sa: float, flags: list[str]) -> str:
    soft = [f for f in flags if f not in HARD_FLAGS]
    if sa > 6.0 or any(f in HARD_FLAGS for f in flags) or len(soft) >= 3:
        return "橙 优先人工评估"
    if sa > 4.0 or soft:
        return "黄 需确认路线"
    return "绿 常规合成"


def catalog_source(cid: str) -> str:
    if cid.startswith("ZINC"):
        return "ZINC（可按需合成目录）"
    if cid.startswith("MCULE-"):
        return "MCULE（可按需合成目录）"
    if cid.startswith("C") and cid[1:].isdigit():
        return "PubChem CID（实体已知化合物）"
    if cid.startswith("SCALE-"):
        return "生成模型输出（无目录来源）"
    return "其他/未知"


def max_fused_ring_system(mol) -> int:
    """最大稠合环系里的环数（共享原子的环做并查集）。"""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if not rings:
        return 0
    parent = list(range(len(rings)))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if rings[i] & rings[j]:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
    groups: dict[int, int] = {}
    for i in range(len(rings)):
        groups[find(i)] = groups.get(find(i), 0) + 1
    return max(groups.values())


def mol_metrics(mol) -> dict:
    rings = [len(r) for r in mol.GetRingInfo().AtomRings()]
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    n_spec = sum(1 for _, tag in centers if tag in ("R", "S"))
    n_unspec = sum(1 for _, tag in centers if tag == "?")
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    return {
        "n_heavy": mol.GetNumHeavyAtoms(),
        "n_stereo_specified": n_spec,
        "n_stereo_unspecified": n_unspec,
        "n_stereo_total": n_spec + n_unspec,
        "n_bridgehead": int(rdMolDescriptors.CalcNumBridgeheadAtoms(mol)),
        "n_spiro": int(rdMolDescriptors.CalcNumSpiroAtoms(mol)),
        "n_rings": int(rdMolDescriptors.CalcNumRings(mol)),
        "n_aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "max_ring_size": int(max(rings)) if rings else 0,
        "max_fused_ring_system": max_fused_ring_system(mol),
        "rare_elements": ";".join(sorted(elements - ALLOWED_ELEMENTS)),
    }


def flags_of(m: dict) -> list[str]:
    f = []
    if m["n_stereo_total"] >= 4:
        f.append("stereo_ge4")
    if m["n_bridgehead"] > 0:
        f.append("bridgehead")
    if m["max_ring_size"] >= 12:
        f.append("macrocycle")
    if m["rare_elements"]:
        f.append("rare_element")
    if 2 <= m["n_stereo_total"] <= 3:
        f.append("stereo_2_3")
    if m["max_fused_ring_system"] >= 3:
        f.append("fused_tricyclic")
    if m["n_spiro"] >= 2:
        f.append("spiro")
    if m["n_rings"] >= 5:
        f.append("ring_rich")
    return f


def main() -> int:
    font = plots.setup_cjk()
    lst = pd.read_csv(LIST_CSV)
    batch = pd.read_csv(BATCH_CSV)
    batch_cids = set(batch.cid.astype(str))
    rows = []
    for r in lst.itertuples():
        mol = Chem.MolFromSmiles(r.smiles)
        if mol is None:
            raise RuntimeError(f"SMILES 解析失败：{r.cid}")
        try:
            sa = float(sascorer.calculateScore(mol))
        except Exception as exc:                                   # pragma: no cover
            raise RuntimeError(f"SAScore 计算失败：{r.cid} / {exc}") from exc
        m = mol_metrics(mol)
        fl = flags_of(m)
        is_ctrl = str(getattr(r, "candidate_type", "")) == "negative_control"
        rows.append({
            "cid": r.cid,
            "pick_order": int(r.pick_order),
            "candidate_type": r.candidate_type,
            "control_role": "" if pd.isna(getattr(r, "control_role", np.nan)) else r.control_role,
            "in_md_batch": int(r.cid in batch_cids),
            "CNN_VS": _f(r.CNN_VS),
            "cns_mpo": _f(r.cns_mpo),
            "MW": _f(r.MW), "qed": _f(r.qed),
            "smiles": r.smiles,
            "catalog_source": catalog_source(str(r.cid)),
            "sa_score": round(sa, 3),
            "sa_band": sa_band(sa),
            **m,
            "n_flags": len(fl),
            "n_soft_flags": len([x for x in fl if x not in HARD_FLAGS]),
            "complex_flags": ";".join(fl),
            "has_hard_flag": int(any(x in HARD_FLAGS for x in fl)),
            "triage_band": "（对照不评）" if is_ctrl else triage(sa, fl),
        })
    df = pd.DataFrame(rows)
    cand = df[df.candidate_type == "candidate"].copy()
    ctrl = df[df.candidate_type != "candidate"].copy()

    # ---------------- 统计
    tri_counts = {k: int((cand.triage_band == k).sum()) for k in TRIAGE_COLOR}
    band_counts = {name: int((cand.sa_band == name).sum()) for _, _, name in SA_BANDS}
    flag_counts = {k: int(cand.complex_flags.str.contains(k).sum()) for k in FLAG_DEFS}
    md = cand[cand.in_md_batch == 1]
    n_cand = int(len(cand))
    n_ge4 = int((cand.sa_score > 4.0).sum())
    n_ge6 = int((cand.sa_score > 6.0).sum())
    hardest = cand.nlargest(min(5, n_cand), "sa_score")[
        ["cid", "sa_score", "sa_band", "complex_flags", "triage_band", "CNN_VS"]]
    # 内部一致性补算（SOP_5 步骤 16 亦做过的同名补算，本层落盘；不是第四套口径）：
    # 名单内 SA 与结构描述符的 Spearman ρ，用来判断这个窄区间里的 SA 排序有没有可解释性。
    from scipy.stats import spearmanr
    corr_cols = {"n_heavy": cand.n_heavy, "n_rings": cand.n_rings, "MW": cand.MW,
                 "CNN_VS": cand.CNN_VS, "qed": cand.qed, "cns_mpo": cand.cns_mpo}
    sa_corr = {}
    for k, v in corr_cols.items():
        vv = pd.to_numeric(v, errors="coerce")
        ok = vv.notna() & cand.sa_score.notna()
        if ok.sum() >= 5 and vv[ok].nunique() > 1:
            rho = float(spearmanr(cand.sa_score[ok], vv[ok]).statistic)
            sa_corr[k] = round(rho, 3)
        else:
            sa_corr[k] = None
    stats = {
        "n_rows": int(len(df)), "n_candidates": n_cand, "n_controls": int(len(ctrl)),
        "sa_candidates": {
            "mean": round(float(cand.sa_score.mean()), 3),
            "median": round(float(cand.sa_score.median()), 3),
            "min": round(float(cand.sa_score.min()), 3),
            "max": round(float(cand.sa_score.max()), 3),
            "std": round(float(cand.sa_score.std()), 3),
            "p90": round(float(cand.sa_score.quantile(0.90)), 3),
        },
        "sa_band_counts": band_counts,
        "sa_band_degenerate": bool(band_counts[SA_BANDS[0][2]] == n_cand),
        "n_candidates_sa_gt4": n_ge4,
        "n_candidates_sa_gt6": n_ge6,
        "triage_counts": tri_counts,
        "flag_counts": flag_counts,
        "n_flags_hit": int(sum(1 for v in flag_counts.values() if v > 0)),
        "n_candidates_with_any_flag": int((cand.n_flags > 0).sum()),
        "n_candidates_with_hard_flag": int((cand.has_hard_flag == 1).sum()),
        "md_batch": {
            "n": int(len(md)),
            "n_rows_in_file": int(len(batch)),
            "sa_mean": round(float(md.sa_score.mean()), 3) if len(md) else None,
            "sa_max": round(float(md.sa_score.max()), 3) if len(md) else None,
            "n_flagged": int((md.n_flags > 0).sum()),
            "triage_counts": {k: int((md.triage_band == k).sum()) for k in TRIAGE_COLOR},
        },
        "controls": [{"cid": r.cid, "sa_score": r.sa_score, "sa_band": r.sa_band,
                      "complex_flags": r.complex_flags, "catalog_source": r.catalog_source}
                     for r in ctrl.itertuples()],
        "hardest5": hardest.to_dict(orient="records"),
        "sa_spearman_within_list": sa_corr,
        "observation": {
            "sa_top2_cids": [r.cid for r in cand.nlargest(min(2, n_cand), "sa_score").itertuples()],
            "sa_top2_have_flags": bool((cand.nlargest(min(2, n_cand), "sa_score").n_flags > 0).any()),
            "flagged_sa_range": ([round(float(cand[cand.n_flags > 0].sa_score.min()), 3),
                                  round(float(cand[cand.n_flags > 0].sa_score.max()), 3)]
                                 if (cand.n_flags > 0).any() else None),
            "controls_above_all_candidates": [
                r.cid for r in ctrl.itertuples() if len(cand) and r.sa_score > cand.sa_score.max()],
            "note": "SA 与结构旗标两套口径不重合（run05 实测：旗标全 0 命中，旗标口径无信息；"
                    "SA 口径给出排序，但见 sa_spearman_within_list 与对照位置——该窄区间内"
                    "SA 排序不具可合成性含义）",
        },
        "catalog_source_counts": {k: int(v) for k, v in df.catalog_source.value_counts().items()},
    }
    n_scope = f"{n_cand} 人名单 = {n_cand} 候选 + {len(ctrl)} 阴性对照"
    report = {
        "step": "16_synth_access",
        "layer": "派生分析层（零侵入：不改门、不改名单、不改规则版本；不在 12_determinism 的 INCLUDE_DIRS 内）",
        "scope": f"主口径 {n_scope}（04_selection/tables/list_main_QED0.50_anchor4.csv，run05 实测）",
        "caliber_note": "口径与 SOP_5 步骤 16 逐条一致（SA 分带 4.0/6.0、8 条难点线、分诊带阈值均未改）；"
                        "随数据变化的只有实测统计量。相对 SOP_5 的唯一新增是 catalog_source 的 "
                        "SCALE- 分支（run05 名单 cid 为生成模型输出，无目录来源可标）。",
        "method": {
            "sascore": {
                "impl": "RDKit Contrib/SA_Score/sascorer.py（随本机 RDKit 分发）",
                "reference": "Ertl P, Schuffenhauer A. J Cheminform 1:8 (2009)；"
                             "本地 README：<RDKit>/Contrib/SA_Score/README",
                "range": "1（易）→ 10（难），连续复杂度分；官方未给分档",
            },
            "sa_bands": [{"lo": lo, "hi": None if hi > 10 else hi, "name": name}
                         for lo, hi, name in SA_BANDS],
            "flag_defs": {k: {"desc": v[0], "hard": v[1]} for k, v in FLAG_DEFS.items()},
            "hard_flags": list(HARD_FLAGS),
            "triage_rule": {
                "橙 优先人工评估": "SA > 6.0 或 命中任一硬旗标 或 软旗标 ≥ 3",
                "黄 需确认路线": "4.0 < SA ≤ 6.0 或 软旗标 1–2",
                "绿 常规合成": "SA ≤ 4.0 且 无旗标",
                "note": "本层自定分诊带（经验），不是论文判据；定义随本报告落盘以便复核",
            },
            "catalog_note": "cid 前缀只表示来源：ZINC / MCULE 为可按需合成目录，PubChem CID 为实体化合物，"
                            "SCALE- 为生成模型输出（run05：**无目录来源可标**）；"
                            "**不等于现货可得**——本层无法离线核验库存或报价",
        },
        "toolchain": {
            "rdkit": RDKIT_VERSION, "font": font,
            "sascorer_py_sha256": sha256_file(Path(os.path.join(RDConfig.RDContribDir, "SA_Score", "sascorer.py"))),
            "fpscores_sha256": sha256_file(Path(os.path.join(RDConfig.RDContribDir, "SA_Score", "fpscores.pkl.gz"))),
            "input_list_sha256": sha256_file(LIST_CSV),
            "input_batch_sha256": sha256_file(BATCH_CSV),
            "timestamp_policy": "不写时间戳，保证同一份输入两次运行输出逐字节一致",
        },
        "stats": stats,
        "caveats": [
            "SA 是启发式复杂度分：不含路线、原料可得性、保护基策略与手性拆分代价；SA ≤ 4 不等于「立刻可买」",
            "旗标是确定性结构规则，命中的分子并非「不可合成」，只表示需要人工看路线",
            f"8 条难点线中 {8 - stats['n_flags_hit']} 条在本名单上零命中——这本身是读数（该名单未触及那些难点），不是漏检",
            "目录来源不等于现货可得；本名单为 ML 生成分子，**无目录来源可标**（SCALE- 前缀），"
            "2 个阴性对照是 DUD-E 诱饵（PubChem CID，已有实体），SA 分只作参照",
            "本层不进交付 MANIFEST 与确定性哈希范围；重跑整链不会改变本目录产物，反之亦然",
        ],
    }
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    (OUT / "reports").mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "tables/synth_access_list.csv", index=False, encoding="utf-8")
    (OUT / "reports/synth_access.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")

    # ---------------- 三联图（坐标范围随实测数据取，不写死）
    fig, axes = plt.subplots(1, 3, figsize=(7.9, 3.05), dpi=300,
                             gridspec_kw={"width_ratios": [1.05, 1.25, 1.05], "wspace": 0.52})
    # A：SA 分布（按分档着色）+ 对照位置
    ax = axes[0]
    step = 0.25
    lo_e = np.floor(float(cand.sa_score.min()) / step) * step
    hi_e = np.ceil(float(cand.sa_score.max()) / step) * step
    edges = np.arange(lo_e, hi_e + step / 2, step)
    cnt, _ = np.histogram(cand.sa_score, bins=edges)
    for lo, hi, name in SA_BANDS:
        idx = [i for i, c in enumerate(edges[:-1], start=0)
               if lo < c + step / 2 <= hi or (lo == 0.0 and c + step / 2 <= hi)]
        vals = [cnt[i] for i in idx if cnt[i] > 0]
        if vals:
            ax.bar([edges[i] + step / 2 for i in idx if cnt[i] > 0], vals, width=step * 0.84,
                   color=BAND_COLOR[name], alpha=0.9)
    ax.axvline(4.0, color=plots.GREY, ls="--", lw=0.9)
    for r in ctrl.itertuples():
        ax.axvline(r.sa_score, color=plots.GREY, ls=":", lw=0.9)
        ax.annotate("对照", (r.sa_score, ax.get_ylim()[1] * 0.985), rotation=90,
                    fontsize=6.2, color=plots.GREY, ha="right", va="top")
    xlo = min(lo_e, float(ctrl.sa_score.min()) - 0.15) if len(ctrl) else lo_e
    xhi = max(hi_e, float(ctrl.sa_score.max()) + 0.15) if len(ctrl) else hi_e
    ax.set_xlim(xlo, xhi)
    ax.set_xlabel("SAScore（1 易 → 10 难）", fontsize=8.5)
    ax.set_ylabel("候选分子数", fontsize=8.5)
    band_tail = "全部落 ≤ 4.0 的「易」档" if stats["sa_band_degenerate"] else f"> 4.0 的 {n_ge4} 个"
    ax.set_title(f"A {n_cand} 个候选的 SA 分布\n中位 {stats['sa_candidates']['median']} · "
                 f"范围 {stats['sa_candidates']['min']}–{stats['sa_candidates']['max']} · "
                 f"{band_tail}", fontsize=9, pad=6)
    # B：旗标命中频次（含零命中，说明规则覆盖面）
    ax = axes[1]
    names = list(FLAG_DEFS)
    vals = [flag_counts[k] for k in names]
    order = list(np.argsort(vals))
    cols = [plots.RED if names[i] in HARD_FLAGS else plots.BLUE for i in order]
    ax.barh(range(len(order)), [vals[i] for i in order], color=cols, height=0.6)
    for j, i in enumerate(order):
        ax.text(max(vals[i], 0) + 0.25, j, str(vals[i]), va="center", fontsize=7.5,
                color=plots.GREY if vals[i] == 0 else "black")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([names[i] for i in order], fontsize=7.4)
    ax.set_xlim(0, max(vals) * 1.3 + 1.5)
    ax.set_xlabel(f"命中候选数（{n_cand} 个候选）", fontsize=8.5)
    ax.set_title(f"B 8 条难点线命中频次\n红 = 硬旗标（本名单 {stats['n_candidates_with_hard_flag']} 命中）",
                 fontsize=9, pad=6)
    ax.grid(axis="y", visible=False)
    # C：CNN_VS – SA 权衡
    ax = axes[2]
    for band, col in TRIAGE_COLOR.items():
        s = cand[cand.triage_band == band]
        if len(s):
            ax.scatter(s.CNN_VS, s.sa_score, s=26, color=col, alpha=0.9, label=band, zorder=3)
    mdb = cand[cand.in_md_batch == 1]
    if len(mdb):
        ax.scatter(mdb.CNN_VS, mdb.sa_score, s=88, facecolors="none", edgecolors=plots.GREY,
                   linewidths=0.9, label=f"进 MD 批次（{stats['md_batch']['n']} 席）", zorder=4)
    ax.axhline(4.0, color=plots.GREY, ls="--", lw=0.8)
    top2 = cand.nlargest(min(2, n_cand), "sa_score")
    ax.scatter(top2.CNN_VS, top2.sa_score, s=110, facecolors="none", edgecolors=plots.RED,
               linewidths=1.0, zorder=5)
    ax.text(0.015, 0.92, "SA 最高两个：\n" + "\n".join(
        f"{r.cid} {r.sa_score}" for r in top2.itertuples()),
        transform=ax.transAxes, fontsize=6.2, color=plots.RED, va="top", ha="left")
    ax.set_xlabel("CNN_VS（对接打分）", fontsize=8.5)
    ax.set_ylabel("SAScore", fontsize=8.5)
    ax.set_title("C 分数–可合成性权衡", fontsize=9, pad=6)
    ax.legend(fontsize=6.4, frameon=False, loc="lower right", handletextpad=0.4)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.845, bottom=0.145)
    plots.save(fig, OUT / "figures/fig_synth_access.png", dpi=300)

    # ---------------- Markdown 小结
    obs = stats["observation"]
    obs_cids = "、".join(obs["sa_top2_cids"])
    obs_max = " / ".join(f"{r.sa_score}" for r in hardest.head(2).itertuples())
    obs_n = stats["n_candidates_with_any_flag"]
    obs_range = (f"{obs['flagged_sa_range'][0]}–{obs['flagged_sa_range'][1]}"
                 if obs["flagged_sa_range"] else "—")
    flag_rows = "\n".join(
        f"| {k} | {FLAG_DEFS[k][0]} | {flag_counts[k]} | {'硬' if k in HARD_FLAGS else '软'} |"
        for k in FLAG_DEFS)
    corr_txt = "；".join(f"{k} {('+' if v > 0 else '')}{v}" if v is not None else f"{k} —"
                         for k, v in sa_corr.items())
    band_tail_md = f"{band_counts[SA_BANDS[0][2]]}/{n_cand} 落「易」"
    flagged_block = "\n".join(
        f"| {r.cid} | {r.sa_score} | {r.complex_flags} | {r.triage_band} | "
        f"{'是' if r.in_md_batch else '否'} | {r.CNN_VS:.4f} |"
        for r in cand[cand.n_flags > 0].sort_values("sa_score", ascending=False).itertuples())
    if not flagged_block:
        flagged_block = "| — | — | — | — | — | — |"
    cross_caliber = (
        f"**两套口径不重合**：SA 最高的两个候选（{obs_cids}，SA {obs_max}）**无任何结构旗标**；"
        f"被标记的 {obs_n} 个候选 SA 落在 {obs_range} 之间，全部低于它们。"
        if obs_n else
        f"**旗标口径在 run05 名单上是 0 信息**：{n_cand} 个候选无一命中 8 条难点线，"
        f"所以没有任何「SA 与旗标不重合」的对照可比。SA 口径仍给出排序，"
        f"最高两个是 {obs_cids}（SA {obs_max}）；但这个窄区间里的排序**不具可合成性含义**——"
        f"见下条对照。")
    ctrl_note = ("；".join(f"{c['cid']} SA {c['sa_score']}（{c['sa_band'].split('（')[0]}档，"
                           f"{c['complex_flags'] or '无旗标'}）"
                           for c in stats["controls"]) or "（无对照）")
    md_txt = f"""# SOP-05 · run05 名单可合成性评估（步骤 16 · 派生分析层）

**结论一句话**：{n_cand} 个候选的 SAScore 中位 **{stats['sa_candidates']['median']}**（范围
{stats['sa_candidates']['min']}–{stats['sa_candidates']['max']}，std {stats['sa_candidates']['std']}，
{'全部落 ≤ 4.0 的「易」分档' if stats['sa_band_degenerate'] else f'其中 {n_ge4} 个 > 4.0'}），
分诊结果 **绿 {tri_counts['绿 常规合成']} / 黄 {tri_counts['黄 需确认路线']} / 橙 {tri_counts['橙 优先人工评估']}**；
8 条难点线中命中 {stats['n_flags_hit']} 条，共 {obs_n} 个候选被标记，
硬旗标（桥环 / 大环 / 多手性 / 罕见元素）**{stats['n_candidates_with_hard_flag']} 命中**。
进 MD 批次的 {stats['md_batch']['n']} 席 SA 均值 {stats['md_batch']['sa_mean']}、
最大 {stats['md_batch']['sa_max']}，其中 {stats['md_batch']['n_flagged']} 个被标记，分诊
绿 {stats['md_batch']['triage_counts']['绿 常规合成']} / 黄 {stats['md_batch']['triage_counts']['黄 需确认路线']} / 橙 {stats['md_batch']['triage_counts']['橙 优先人工评估']}。
**本层在 run05 名单上无区分度**（全绿、零旗标），因此它是红旗扫描的「未报警」，不是「确认好合成」。

## 1. 三套口径

| 口径 | 实现 | 说明 |
| --- | --- | --- |
| SAScore | RDKit `Contrib/SA_Score/sascorer.py` | Ertl & Schuffenhauer 2009 的复杂度–片段贡献启发式分，1（易）→ 10（难）；官方**未给分档**，本层用 4.0 / 6.0 作经验分带 |
| 复杂旗标 | 8 条确定性结构规则 | 见下表；桥环 / 大环 / 多手性 / 罕见元素为**硬旗标** |
| 来源标注 | cid 前缀 | ZINC / MCULE（可按需合成目录）、PubChem CID（实体化合物）、**SCALE-（生成模型输出，无目录来源）**；**不等于现货可得** |

分诊带（本层自定，随报告落盘，阈值与 SOP_5 步骤 16 逐条一致）：橙 = SA > 6.0 或命中硬旗标或软旗标 ≥ 3；黄 = 4.0 < SA ≤ 6.0 或软旗标 1–2；绿 = SA ≤ 4.0 且无旗标。

## 2. 读数

| SA 分档 | 候选数 |
| --- | --- |
| {SA_BANDS[0][2]} | {band_counts[SA_BANDS[0][2]]} |
| {SA_BANDS[1][2]} | {band_counts[SA_BANDS[1][2]]} |
| {SA_BANDS[2][2]} | {band_counts[SA_BANDS[2][2]]} |

**SA 分档在本名单上不具区分度**（{band_tail_md}）：该分是复杂度–片段贡献的启发式分，
本批分子骨架常规、分子量集中（重原子 {int(cand.n_heavy.min())}–{int(cand.n_heavy.max())}），
所以分化本该来自结构旗标——而旗标在本名单上全为 0。

| 难点线 | 定义 | 命中候选数 | 类型 |
| --- | --- | --- | --- |
{flag_rows}

被标记的候选：

| cid | SAScore | 旗标 | 分诊带 | 进 MD 批次 | CNN_VS |
| --- | --- | --- | --- | --- | --- |
{flagged_block}

{cross_caliber}

2 个阴性对照（DUD-E 诱饵，已有实体，SA 分只作参照、不进分诊统计）：{ctrl_note}。
{'**对照位置是重要读数**：' + '、'.join(obs['controls_above_all_candidates']) + f" 的 SA 高于全部 {n_cand} 个候选——一个已被合成出来的真实化合物在 SA 上比所有候选都「难」，说明本区间内 SA 排序不具备可合成性含义。" if obs['controls_above_all_candidates'] else '两个对照的 SA 均落在候选区间内。'}

**内部一致性补算**（SOP_5 步骤 16 亦做过的同名补算，本层落盘；不是第四套口径）：
名单内 SA 与结构描述符的 Spearman ρ = {corr_txt}。

## 3. 如实说明（写进材料前请保留）

1. **SA 是启发式分**，不含路线设计、原料可得性、保护基策略与手性拆分代价；SA ≤ 4 不等于「立刻可买」。
2. **旗标命中 ≠ 不可合成**：只表示该分子有需要人工看路线的结构特征；run05 名单上 8 条难点线**全部零命中**，
   即这批 ML 生成分子未触及桥环 / 大环 / ≥4 手性中心 / 罕见元素 / 螺环 / 稠合三环 / 环数 ≥ 5 / 手性 2–3
   任何一条——这是读数（名单未触及这些难点），不是漏检。
3. **绿带是无信息，不是正面结论**：{tri_counts['绿 常规合成']} 个候选落绿带只说明「启发式未报警」。
4. **目录来源不等于现货可得**：run05 名单为 ML 生成分子，**没有目录来源可标**（SCALE- 前缀），
   这一点与 SOP_5 名单（ZINC / MCULE / PubChem CID）本质不同——SOP_5 那套「目录来源」读数在本批上不可用。
5. **本层不改名单、不改门、不改规则版本**：run05 主名单（{n_cand} 候选 + {len(ctrl)} 对照）与 MD 批次
   （{stats['md_batch']['n']} 席，文件 {stats['md_batch']['n_rows_in_file']} 行含 2 对照）一字未动，分诊只作标注与人工复核线索。
6. 本目录不在 `12_determinism` 的 `INCLUDE_DIRS` 内，也不进 `PACK_MANIFEST.sha256` / 交付 `MANIFEST.sha256`；
   重跑整链不会改变本目录产物，反之本层重跑也不影响主链确定性判据。

## 4. 出处

- 评估对象：`04_selection/tables/list_main_QED0.50_anchor4.csv`（{stats['n_rows']} 行 = {n_cand} 候选 + {len(ctrl)} 对照）
- MD 批次标记：`07_md_batch/tables/md_batch20.csv`（{stats['md_batch']['n_rows_in_file']} 行，实为 {stats['md_batch']['n']} 候选 + 2 对照）
- 打分实现：RDKit `Contrib/SA_Score/sascorer.py` + `fpscores.pkl.gz`（报告内登记 sha256）
- 口径与方法定义：`09_synth_access/reports/synth_access.json`
- 逐分子读数：`09_synth_access/tables/synth_access_list.csv`
- 图：`09_synth_access/figures/fig_synth_access.png` / `.svg`
- 文献：Ertl P, Schuffenhauer A. *J Cheminform* 1:8 (2009)，http://www.jcheminf.com/content/1/1/8
- 同口径基线：`<包根>/src/SOP_5/09_synth_access/`（SOP_5 名单 48 候选：中位 2.726、范围 2.197–3.547、绿 40 / 黄 8 / 橙 0）
"""
    (OUT / "reports/SYNTH_ACCESS.md").write_text(md_txt, encoding="utf-8")

    print(f"[16_synth_access] 字体 {font} · 行数 {len(df)}（候选 {n_cand} / 对照 {len(ctrl)}）")
    print(f"[16_synth_access] SA 中位 {stats['sa_candidates']['median']} · "
          f"范围 {stats['sa_candidates']['min']}–{stats['sa_candidates']['max']} · 分档 {band_counts}")
    print(f"[16_synth_access] 分诊 {tri_counts} · 旗标命中 {stats['n_candidates_with_any_flag']}"
          f"（硬旗标 {stats['n_candidates_with_hard_flag']}）· 频次 {flag_counts}")
    print(f"[16_synth_access] MD {stats['md_batch']['n']} 席：SA 均值 {stats['md_batch']['sa_mean']} / 最大 "
          f"{stats['md_batch']['sa_max']} · 标记 {stats['md_batch']['n_flagged']} · "
          f"{stats['md_batch']['triage_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
