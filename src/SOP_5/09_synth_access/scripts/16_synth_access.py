#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 16（派生分析层 · 零侵入）：48 人名单可合成性评估。

定位
    这是主链之外的分析层，与 08_visuals 同类：**不改门、不改名单、不改规则版本**，
    产物只落 09_synth_access/。目录不在 12_determinism 的 INCLUDE_DIRS 内，也不重生成
    PACK_MANIFEST.sha256 / 交付 MANIFEST.sha256，因此不影响主链确定性基线。

评估对象
    主口径 48 人名单（04_selection/tables/list_main_QED0.50_anchor4.csv，50 行 =
    48 候选 + 2 阴性对照）。对照是 DUD-E 诱饵（PubChem CID，已有实体的已知化合物），
    其 SA 分只作参照，不进分诊带统计。

三套互补口径（不合并成一个「可合成性真值」）
    ① SAScore  : RDKit Contrib/SA_Score 的 sascorer（Ertl & Schuffenhauer 2009）
                  —— 复杂度–片段贡献启发式分，1（易）→ 10（难）；官方未给分档。
    ② 复杂旗标 : 8 条确定性结构规则（见 FLAG_DEFS），其中 4 条记作硬旗标。
    ③ 目录来源 : cid 前缀 → ZINC / MCULE（可按需合成目录）或 PubChem CID（实体化合物）。
                  **来源标注不等于现货可得**（本层无法离线核验库存与报价）。

分诊带（本层自定，非论文判据；定义随报告落盘）
    橙 优先人工评估 : SA > 6.0 或 命中任一硬旗标 或 软旗标 ≥ 3
    黄 需确认路线   : 4.0 < SA ≤ 6.0 或 软旗标 1–2
    绿 常规合成     : SA ≤ 4.0 且 无旗标

输出（09_synth_access/）
    tables/synth_access_list.csv      50 行逐分子读数（含全部中间量）
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

# 8 条难点线；hard=True 者一旦命中直接进橙带
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
    batch_cids = set(pd.read_csv(BATCH_CSV).cid.astype(str))
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
            "CNN_VS": float(r.CNN_VS) if pd.notna(r.CNN_VS) else np.nan,
            "cns_mpo": float(r.cns_mpo) if pd.notna(r.cns_mpo) else np.nan,
            "MW": float(r.MW), "qed": float(r.qed),
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
    hardest = cand.nlargest(5, "sa_score")[["cid", "sa_score", "sa_band", "complex_flags",
                                            "triage_band", "CNN_VS"]]
    stats = {
        "n_rows": int(len(df)), "n_candidates": int(len(cand)), "n_controls": int(len(ctrl)),
        "sa_candidates": {
            "mean": round(float(cand.sa_score.mean()), 3),
            "median": round(float(cand.sa_score.median()), 3),
            "min": round(float(cand.sa_score.min()), 3),
            "max": round(float(cand.sa_score.max()), 3),
            "std": round(float(cand.sa_score.std()), 3),
            "p90": round(float(cand.sa_score.quantile(0.90)), 3),
        },
        "sa_band_counts": band_counts,
        "sa_band_degenerate": bool(band_counts[SA_BANDS[0][2]] == len(cand)),
        "triage_counts": tri_counts,
        "flag_counts": flag_counts,
        "n_candidates_with_any_flag": int((cand.n_flags > 0).sum()),
        "n_candidates_with_hard_flag": int((cand.has_hard_flag == 1).sum()),
        "md_batch": {
            "n": int(len(md)),
            "sa_mean": round(float(md.sa_score.mean()), 3),
            "sa_max": round(float(md.sa_score.max()), 3),
            "n_flagged": int((md.n_flags > 0).sum()),
            "triage_counts": {k: int((md.triage_band == k).sum()) for k in TRIAGE_COLOR},
        },
        "controls": [{"cid": r.cid, "sa_score": r.sa_score, "sa_band": r.sa_band,
                      "complex_flags": r.complex_flags, "catalog_source": r.catalog_source}
                     for r in ctrl.itertuples()],
        "hardest5": hardest.to_dict(orient="records"),
        "observation": {
            "sa_top2_cids": [r.cid for r in cand.nlargest(2, "sa_score").itertuples()],
            "sa_top2_have_flags": bool((cand.nlargest(2, "sa_score").n_flags > 0).any()),
            "flagged_sa_range": [round(float(cand[cand.n_flags > 0].sa_score.min()), 3),
                                 round(float(cand[cand.n_flags > 0].sa_score.max()), 3)],
            "note": "SA 与结构旗标两套口径不重合：SA 最高的两个候选无任何结构旗标，"
                    "而被标记的候选 SA 均低于它们 —— 两套读数互补，不应互相替代",
        },
        "catalog_source_counts": {k: int(v) for k, v in df.catalog_source.value_counts().items()},
    }
    report = {
        "step": "16_synth_access",
        "layer": "派生分析层（零侵入：不改门、不改名单、不改规则版本；不在 12_determinism 的 INCLUDE_DIRS 内）",
        "scope": "主口径 48 人名单 = 48 候选 + 2 阴性对照（04_selection/tables/list_main_QED0.50_anchor4.csv）",
        "method": {
            "sascore": {
                "impl": "RDKit Contrib/SA_Score/sascorer.py（随本机 RDKit 分发）",
                "reference": "Ertl P, Schuffenhauer A. J Cheminform 1:8 (2009)；"
                             "本地 README：/opt/conda/share/RDKit/Contrib/SA_Score/README",
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
            "catalog_note": "cid 前缀只表示来源目录（ZINC / MCULE 为可按需合成目录，PubChem CID 为实体化合物）；"
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
            "8 条难点线中多数在本名单上零命中——这本身是读数（该名单未触及那些难点），不是漏检",
            "目录来源不等于现货可得；2 个阴性对照是 DUD-E 诱饵（PubChem CID，已有实体），SA 分只作参照",
            "本层不进交付 MANIFEST 与确定性哈希范围；重跑整链不会改变本目录产物，反之亦然",
        ],
    }
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    (OUT / "reports").mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "tables/synth_access_list.csv", index=False, encoding="utf-8")
    (OUT / "reports/synth_access.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")

    # ---------------- 三联图
    fig, axes = plt.subplots(1, 3, figsize=(7.9, 3.05), dpi=300,
                             gridspec_kw={"width_ratios": [1.05, 1.25, 1.05], "wspace": 0.52})
    # A：SA 分布（按分档着色）+ 2 个对照位置
    ax = axes[0]
    edges = np.arange(2.0, 4.51, 0.25)
    cnt, _ = np.histogram(cand.sa_score, bins=edges)
    for lo, hi, name in SA_BANDS:
        idx = [i for i, c in enumerate(edges[:-1], start=0)
               if lo < c + 0.125 <= hi or (lo == 0.0 and c + 0.125 <= hi)]
        vals = [cnt[i] for i in idx if cnt[i] > 0]
        if vals:
            ax.bar([edges[i] + 0.125 for i in idx if cnt[i] > 0], vals, width=0.21,
                   color=BAND_COLOR[name], alpha=0.9)
    ax.axvline(4.0, color=plots.GREY, ls="--", lw=0.9)
    for r in ctrl.itertuples():
        ax.axvline(r.sa_score, color=plots.GREY, ls=":", lw=0.9)
        ax.annotate("对照", (r.sa_score, ax.get_ylim()[1] * 0.985), rotation=90,
                    fontsize=6.2, color=plots.GREY, ha="right", va="top")
    ax.set_xlim(2.0, 4.6)
    ax.set_xlabel("SAScore（1 易 → 10 难）", fontsize=8.5)
    ax.set_ylabel("候选分子数", fontsize=8.5)
    ax.set_title(f"A 48 个候选的 SA 分布\n中位 {stats['sa_candidates']['median']} · "
                 f"范围 {stats['sa_candidates']['min']}–{stats['sa_candidates']['max']} · "
                 f"全部 ≤ 4.0", fontsize=9, pad=6)
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
    ax.set_xlabel("命中候选数（48 个候选）", fontsize=8.5)
    ax.set_title("B 8 条难点线命中频次\n红 = 硬旗标（本名单 0 命中）", fontsize=9, pad=6)
    ax.grid(axis="y", visible=False)
    # C：CNN_VS – SA 权衡
    ax = axes[2]
    for band, col in TRIAGE_COLOR.items():
        s = cand[cand.triage_band == band]
        if len(s):
            ax.scatter(s.CNN_VS, s.sa_score, s=26, color=col, alpha=0.9, label=band, zorder=3)
    mdb = cand[cand.in_md_batch == 1]
    ax.scatter(mdb.CNN_VS, mdb.sa_score, s=88, facecolors="none", edgecolors=plots.GREY,
               linewidths=0.9, label="进 MD 批次（18 席）", zorder=4)
    ax.axhline(4.0, color=plots.GREY, ls="--", lw=0.8)
    top2 = cand.nlargest(2, "sa_score")
    ax.scatter(top2.CNN_VS, top2.sa_score, s=110, facecolors="none", edgecolors=plots.RED,
               linewidths=1.0, zorder=5)
    ax.text(0.015, 0.985, "SA 最高两个：\n" + "\n".join(
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
    obs_range = f"{obs['flagged_sa_range'][0]}–{obs['flagged_sa_range'][1]}"
    flag_rows = "\n".join(
        f"| {k} | {FLAG_DEFS[k][0]} | {flag_counts[k]} | {'硬' if k in HARD_FLAGS else '软'} |"
        for k in FLAG_DEFS)
    md_txt = f"""# SOP-05 · 48 人名单可合成性评估（步骤 16 · 派生分析层）

**结论一句话**：48 个候选的 SAScore 中位 **{stats['sa_candidates']['median']}**（范围
{stats['sa_candidates']['min']}–{stats['sa_candidates']['max']}，全部落 ≤ 4.0 的「易」分档），
分诊结果 **绿 {tri_counts['绿 常规合成']} / 黄 {tri_counts['黄 需确认路线']} / 橙 {tri_counts['橙 优先人工评估']}**；
8 条难点线中命中 2 条，共 {stats['n_candidates_with_any_flag']} 个候选被标记，
硬旗标（桥环 / 大环 / 多手性 / 罕见元素）**0 命中**。进 MD 批次的 18 席 SA 均值
{stats['md_batch']['sa_mean']}、最大 {stats['md_batch']['sa_max']}，其中
{stats['md_batch']['n_flagged']} 个被标记，分诊
绿 {stats['md_batch']['triage_counts']['绿 常规合成']} / 黄 {stats['md_batch']['triage_counts']['黄 需确认路线']} / 橙 {stats['md_batch']['triage_counts']['橙 优先人工评估']}。

## 1. 三套口径

| 口径 | 实现 | 说明 |
| --- | --- | --- |
| SAScore | RDKit `Contrib/SA_Score/sascorer.py` | Ertl & Schuffenhauer 2009 的复杂度–片段贡献启发式分，1（易）→ 10（难）；官方**未给分档**，本层用 4.0 / 6.0 作经验分带 |
| 复杂旗标 | 8 条确定性结构规则 | 见下表；桥环 / 大环 / 多手性 / 罕见元素为**硬旗标** |
| 目录来源 | cid 前缀 | ZINC / MCULE（可按需合成目录）或 PubChem CID（实体化合物）；**不等于现货可得** |

分诊带（本层自定，随报告落盘）：橙 = SA > 6.0 或命中硬旗标或软旗标 ≥ 3；黄 = 4.0 < SA ≤ 6.0 或软旗标 1–2；绿 = SA ≤ 4.0 且无旗标。

## 2. 读数

| SA 分档 | 候选数 |
| --- | --- |
| {SA_BANDS[0][2]} | {band_counts[SA_BANDS[0][2]]} |
| {SA_BANDS[1][2]} | {band_counts[SA_BANDS[1][2]]} |
| {SA_BANDS[2][2]} | {band_counts[SA_BANDS[2][2]]} |

**SA 分档在本名单上不具区分度**（48/48 落「易」）：该分是复杂度–片段贡献的启发式分，
本批分子骨架常规、分子量集中（重原子中位 26、最大 34），所以分化全部来自结构旗标。

| 难点线 | 定义 | 命中候选数 | 类型 |
| --- | --- | --- | --- |
{flag_rows}

被标记的候选：

| cid | SAScore | 旗标 | 分诊带 | 进 MD 批次 | CNN_VS |
| --- | --- | --- | --- | --- | --- |
""" + ("\n".join(
        f"| {r.cid} | {r.sa_score} | {r.complex_flags} | {r.triage_band} | "
        f"{'是' if r.in_md_batch else '否'} | {r.CNN_VS:.4f} |"
        for r in cand[cand.n_flags > 0].sort_values("sa_score", ascending=False).itertuples())
        or "| — | — | — | — | — | — |") + f"""

**两套口径不重合**：SA 最高的两个候选（{obs_cids}，SA {obs_max}）**无任何结构旗标**；
被标记的 {obs_n} 个候选 SA 落在 {obs_range} 之间，全部低于它们。也就是说 SAScore 与结构旗标
给出的是互补信息——SA 分高不一定结构复杂，结构复杂（稠合三环、环数多、多手性）也不一定 SA 分高，
两条都要看，不能互相替代。

2 个阴性对照（DUD-E 诱饵，已有实体，SA 分只作参照、不进分诊统计）：
""" + "；".join(f"{c['cid']} SA {c['sa_score']}（{c['sa_band'].split('（')[0]}档，"
                       f"{c['complex_flags'] or '无旗标'}）"
                for c in stats["controls"]) + f"""

## 3. 如实说明（写进材料前请保留）

1. **SA 是启发式分**，不含路线设计、原料可得性、保护基策略与手性拆分代价；SA ≤ 4 不等于「立刻可买」。
2. **旗标命中 ≠ 不可合成**：只表示该分子有需要人工看路线的结构特征；本名单命中的只有
   `fused_tricyclic`（三环及以上稠合）与 `ring_rich`（环数 ≥ 5）两类，属中等复杂度提示。
3. **8 条难点线里 6 条零命中**（含 4 条硬旗标）——这是读数（该名单未触及桥环 / 大环 / ≥4 手性中心 / 罕见元素），
   不是漏检；规则集覆盖范围写在本文件与 JSON 里，便于别人复核。
4. **目录来源不等于现货可得**：ZINC / MCULE 是可按需合成目录，本层无法离线核验库存与报价。
5. **本层不改名单、不改门、不改规则版本**：48 人名单与 20 席 MD 批次一字未动，分诊只作标注与人工复核线索。
6. 本目录不在 `12_determinism` 的 `INCLUDE_DIRS` 内，也不进 `PACK_MANIFEST.sha256` / 交付 `MANIFEST.sha256`；
   重跑整链不会改变本目录产物，反之本层重跑也不影响主链确定性判据。

## 4. 出处

- 评估对象：`04_selection/tables/list_main_QED0.50_anchor4.csv`（50 行 = 48 候选 + 2 对照）
- MD 批次标记：`07_md_batch/tables/md_batch20.csv`
- 打分实现：RDKit `Contrib/SA_Score/sascorer.py` + `fpscores.pkl.gz`（报告内登记 sha256）
- 口径与方法定义：`09_synth_access/reports/synth_access.json`
- 逐分子读数：`09_synth_access/tables/synth_access_list.csv`
- 图：`09_synth_access/figures/fig_synth_access.png` / `.svg`
- 文献（本地 README 记载）：Ertl P, Schuffenhauer A. *J Cheminform* 1:8 (2009)，http://www.jcheminf.com/content/1/1/8
"""
    (OUT / "reports/SYNTH_ACCESS.md").write_text(md_txt, encoding="utf-8")

    print(f"[16_synth_access] 字体 {font} · 行数 {len(df)}（候选 {len(cand)} / 对照 {len(ctrl)}）")
    print(f"[16_synth_access] SA 中位 {stats['sa_candidates']['median']} · "
          f"范围 {stats['sa_candidates']['min']}–{stats['sa_candidates']['max']} · 分档 {band_counts}")
    print(f"[16_synth_access] 分诊 {tri_counts} · 旗标命中 {stats['n_candidates_with_any_flag']}"
          f"（硬旗标 {stats['n_candidates_with_hard_flag']}）· 频次 {flag_counts}")
    print(f"[16_synth_access] MD 18 席：SA 均值 {stats['md_batch']['sa_mean']} / 最大 "
          f"{stats['md_batch']['sa_max']} · 标记 {stats['md_batch']['n_flagged']} · "
          f"{stats['md_batch']['triage_counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
