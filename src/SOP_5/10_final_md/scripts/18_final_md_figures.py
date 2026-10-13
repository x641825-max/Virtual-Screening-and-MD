#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 17 的图：最终 MD 候选结构卡（22 张）+ 批次构成图。

输入：tables/final_md_20.csv、reports/final_md.json、两份上游 MD 批次表（取资格池分数）
输出：figures/fig_final_md_cards.png(+.svg)、figures/fig_final_md_batch.png(+.svg)

卡面沿用 SOP_5 08_visuals/figures/fig01_md_structure_cards 的版式（RDKit 描画、4 列、
CNS/Peripheral 标注、hERG 第二模型标记），新增两项：**来源标注**（目录 / 生成）与
**可合成性分诊带**。编号 A2AR-MD-01…22。
"""
from __future__ import annotations

import io
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

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402

from vendor import plots  # noqa: E402
from vendor.io_utils import read_json  # noqa: E402

RDLogger.DisableLog("rdApp.*")
ML_ROOT = Path(os.environ.get("SOPML_ROOT", str(ROOT.parent / "SOP_ML_run05")))
OUT = ROOT / "10_final_md"
FIG = OUT / "figures"

BLUE, GREY, RED, GREEN, AMBER = plots.BLUE, plots.GREY, plots.RED, plots.GREEN, plots.AMBER
INK, MUTED = "#0B0C0E", "#6B7079"
SRC_DIR = "#1F7A5C"      # 目录来源
SRC_GEN = "#8A4FBE"      # ML 生成来源
WARN_TX, OK_TX = "#B07908", "#15915F"


def mol_png(smiles: str, w: int = 640, h: int = 508) -> bytes:
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)
    dr = rdMolDraw2D.MolDraw2DCairo(w, h)
    o = dr.drawOptions()
    o.bondLineWidth = 2
    o.padding = 0.05
    o.fixedBondLength = 34
    rdMolDraw2D.PrepareAndDrawMolecule(dr, m)
    dr.FinishDrawing()
    return dr.GetDrawingText()


# ---------------------------------------------------------------- 图 1：结构卡
def fig_cards(roster: pd.DataFrame) -> dict:
    b = roster.sort_values("batch_order").reset_index(drop=True)
    n = len(b)
    rows, cols = 6, 4
    fig, axes = plt.subplots(rows, cols, figsize=(13.8, 19.8))
    fig.subplots_adjust(left=0.015, right=0.985, top=0.945, bottom=0.010,
                        hspace=0.06, wspace=0.04)
    is_cand = b.candidate_type == "candidate"
    n_cns = int(np.sum(is_cand & (b.cns_mpo >= 4).to_numpy()))
    n_per = int(is_cand.sum()) - n_cns
    n_dir = int((b.source_short == "目录").sum())
    n_gen = int((b.source_short == "生成").sum())
    n_hl = int(b.herg_cross_model_flag.fillna("").eq("cross_model_fail").sum())
    n_yellow = int(b.triage_band.fillna("").str.startswith("黄").sum())

    fig.text(0.015, 0.988, "SOP-05 · 最终 MD 候选 20 席结构卡（合并两分支批次）",
             fontsize=16.5, color=INK, ha="left", va="top")
    fig.text(0.015, 0.967,
             f"共 {n} 格 = {int(is_cand.sum())} 候选（目录 {n_dir} + ML 生成 {n_gen}） + "
             f"2 共用阴性对照｜按 CNN_VS 降序（对照固定占末两位）｜"
             f"CNS = CNS MPO ≥ 4：{n_cns} 个；Peripheral < 4：{n_per} 个",
             fontsize=10, color=MUTED, ha="left", va="top")
    fig.text(0.015, 0.9505,
             f"黄框 = 可合成性黄带「需确认路线」：{n_yellow} 个（绿带是「启发式未报警」，"
             f"不是「确认好合成」）｜△ = hERG 第二模型（ChEMBL 自建集成）分歧：保留并标注，"
             f"共 {n_hl} 个",
             fontsize=9.2, color=MUTED, ha="left", va="top")
    fig.text(0.985, 0.988, "来源：■ 目录 ZINC/MCULE　■ ML 生成 SCALE-",
             fontsize=9.2, color=MUTED, ha="right", va="top")

    for i in range(rows * cols):
        ax = axes.ravel()[i]
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        if i >= n:
            continue
        r = b.iloc[i]
        flagged = str(r.herg_cross_model_flag or "") == "cross_model_fail"
        yellow = isinstance(r.triage_band, str) and r.triage_band.startswith("黄")
        if isinstance(r.smiles, str) and r.smiles:
            ax.imshow(plt.imread(io.BytesIO(mol_png(r.smiles)), format="png"),
                      extent=[0.02, 0.98, 0.27, 1.0], aspect="auto")
        src_col = {"目录": SRC_DIR, "生成": SRC_GEN}.get(r.source_short, GREY)
        if r.candidate_type == "candidate":
            cls = "CNS" if float(r.cns_mpo) >= 4 else "Peripheral"
        else:
            cls = "对照 " + str(r.control_role or "")
        ax.text(0.5, 0.225, f"{r.md_card_id} | {cls}", ha="center", va="top",
                fontsize=10.5, color=WARN_TX if flagged else "#2347C4")
        if flagged:
            ax.text(0.985, 0.235, "△", ha="right", va="top", fontsize=11, color=WARN_TX)
        ax.text(0.5, 0.135, str(r.cid), ha="center", va="top", fontsize=9.5,
                color=INK, family="monospace")
        if r.candidate_type == "candidate":
            ax.text(0.5, 0.058,
                    f"来源 {r.source_short}　CNN_VS {float(r.CNN_VS):.3f}　"
                    f"SA {float(r.sa_score):.2f}",
                    ha="center", va="top", fontsize=8.6, color=src_col)
            ax.text(0.5, -0.012, "黄 需确认路线" if yellow else "绿 常规合成",
                    ha="center", va="top", fontsize=8.6,
                    color=WARN_TX if yellow else OK_TX)
        else:
            ax.text(0.5, 0.058, f"来源 两分支共用　SA {float(r.sa_score):.2f}",
                    ha="center", va="top", fontsize=8.6, color=GREY)
            ax.text(0.5, -0.012, "对照不评可合成性", ha="center", va="top",
                    fontsize=8.6, color=GREY)
        if yellow:
            ax.add_patch(FancyBboxPatch((0.04, -0.065), 0.92, 0.325,
                                        boxstyle="round,pad=0,rounding_size=0.03",
                                        linewidth=1.2, edgecolor=AMBER,
                                        facecolor="none", transform=ax.transAxes))
    plots.save(fig, FIG / "fig_final_md_cards.png", dpi=200)
    return {"n_cards": n, "n_candidates": int(is_cand.sum()), "n_dir": n_dir, "n_gen": n_gen,
            "n_cns": n_cns, "n_peripheral": n_per, "n_herg_cross_flag": n_hl,
            "n_synth_yellow": n_yellow}


# ---------------------------------------------------------------- 图 2：批次构成
def _pool_table() -> pd.DataFrame:
    """资格池逐分子分数：直接读两份上游 MD 批次表（与选席脚本同源）。"""
    ups = []
    for p in [ROOT / "07_md_batch/tables/md_batch20.csv",
              ML_ROOT / "07_md_batch/tables/md_batch20.csv"]:
        d = pd.read_csv(p)
        d["src"] = np.where(d.cid.astype(str).str.startswith(("ZINC", "MCULE")), "目录", "生成")
        ups.append(d[d.candidate_type == "candidate"])
    return pd.concat(ups, ignore_index=True)


def fig_batch(roster: pd.DataFrame, rep: dict) -> dict:
    b = roster[roster.candidate_type == "candidate"].copy()
    src_counts = rep["roster"]["source_counts"]
    pool_src = rep["eligibility"]["pool_source_counts"]
    bins = rep["optimality"]["bin_counts"]
    src5, srcm = list(rep["merge"]["source_counts"].keys())
    sens = rep["sensitivity"]

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.6))
    fig.subplots_adjust(left=0.055, right=0.985, top=0.78, bottom=0.17, wspace=0.30)

    # A 来源构成：资格池 vs 本层名册 vs 敏感性对照名册
    ax = axes[0]
    labels = ["资格池（去重后 33）", "本层名册（20 席）", "敏感性对照（20 席）"]
    dirs = [pool_src.get(src5, 0), src_counts.get(src5, 0), sens["source_counts"].get(src5, 0)]
    gens = [pool_src.get(srcm, 0), src_counts.get(srcm, 0), sens["source_counts"].get(srcm, 0)]
    y = np.arange(3)[::-1]
    ax.barh(y, dirs, color=SRC_DIR, height=0.5, label="目录 ZINC/MCULE")
    ax.barh(y, gens, left=dirs, color=SRC_GEN, height=0.5, label="ML 生成 SCALE-")
    for yy, dv, gv in zip(y, dirs, gens):
        ax.text(dv / 2, yy, str(dv), ha="center", va="center", fontsize=9.5, color="white")
        ax.text(dv + gv / 2, yy, str(gv), ha="center", va="center", fontsize=9.5, color="white")
        ax.text(dv + gv + 0.6, yy, f"{dv} / {gv}", ha="left", va="center",
                fontsize=8.6, color=MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.8)
    ax.set_xlim(0, max(d + g for d, g in zip(dirs, gens)) + 7)
    ax.set_xlabel("席位数", fontsize=9)
    ax.set_title("A 来源构成：合并池 → 名册", fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="lower right")

    # B 三轴分位覆盖：12 格各 5 席（理论最大 = 20×3/12）
    ax = axes[1]
    names = [f"CNN_VS|Q{i}" for i in range(1, 5)] + \
            [f"cLogP|Q{i}" for i in range(1, 5)] + \
            [f"CNS MPO|Q{i}" for i in range(1, 5)]
    vals = [bins.get(nm, 0) for nm in names]
    cols = [BLUE] * 4 + [GREEN] * 4 + [AMBER] * 4
    xs = np.arange(12)
    ax.bar(xs, vals, color=cols, width=0.68)
    ax.axhline(5, ls="--", lw=1.0, color=RED)
    ax.text(-0.45, 5.18, "理论最大 20×3/12 = 5", ha="left", va="bottom",
            fontsize=8.6, color=RED)
    for x, v in zip(xs, vals):
        ax.text(x, v - 0.42, str(v), ha="center", va="top", fontsize=9, color="white")
    ax.set_xticks(xs)
    ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4"] * 3, fontsize=8.5)
    for nm, xc, col in [("CNN_VS", 1.5, BLUE), ("cLogP", 5.5, GREEN), ("CNS MPO", 9.5, AMBER)]:
        ax.text(xc, -0.145, nm, ha="center", va="top", fontsize=9.2, color=col,
                transform=ax.get_xaxis_transform())
    ax.set_ylim(0, 7.4)
    ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
    ax.set_ylabel("席位数", fontsize=9)
    ax.set_title(f"B 三轴分位覆盖：t* = {rep['optimality']['t_star']}，12/12 格铺满", fontsize=11)

    # C 名册在 CNN_VS 轴上的铺开
    ax = axes[2]
    up = _pool_table()
    sel = set(b.cid)
    up["in_roster"] = up.cid.isin(sel)
    out = up[~up.in_roster]
    ax.scatter(out.CNN_VS, np.full(len(out), 0.12), s=26, color="#C9CED6",
               edgecolor="#9AA2AC", lw=0.5, label=f"资格池未入席（{len(out)}）", zorder=2)
    for src, col, yy in [("目录", SRC_DIR, 0.38), ("生成", SRC_GEN, 0.62)]:
        s = up[up.in_roster & (up.src == src)]
        ax.scatter(s.CNN_VS, np.full(len(s), yy), s=52, color=col,
                   label=f"入席·{src}（{len(s)}）", zorder=3)
    zmin, zmax = float(b.CNN_VS.min()), float(b.CNN_VS.max())
    ax.plot([zmin, zmax], [0.86, 0.86], color=RED, lw=1.4, zorder=4)
    ax.text((zmin + zmax) / 2, 0.90, f"名册跨度 {zmax - zmin:.3f}", ha="center", va="bottom",
            fontsize=8.8, color=RED)
    cid_max = up.loc[up.CNN_VS.idxmax(), "cid"]
    cid_min = up.loc[up.CNN_VS.idxmin(), "cid"]
    ax.annotate(f"跨度端点 {cid_max}\n{zmax:.3f}", xy=(zmax, 0.62), xytext=(zmax - 0.62, 1.07),
                fontsize=7.8, color=MUTED, ha="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
    ax.annotate(f"跨度端点 {cid_min}\n{zmin:.3f}", xy=(zmin, 0.62), xytext=(zmin + 0.72, 1.07),
                fontsize=7.8, color=MUTED, ha="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
    ax.set_xlim(float(out.CNN_VS.min()) - 0.42, float(out.CNN_VS.max()) + 0.42)
    ax.set_yticks([0.12, 0.38, 0.62, 0.86])
    ax.set_yticklabels(["资格池\n未入席", "入席\n目录", "入席\n生成", "名册\n跨度"], fontsize=8)
    ax.set_ylim(0, 1.32)
    ax.set_xlabel("CNN_VS（对接打分）", fontsize=9)
    ax.set_title("C 名册在 CNN_VS 轴上的铺开", fontsize=11)
    ax.legend(fontsize=7.6, frameon=False, loc="lower right")

    fig.text(0.015, 0.955,
             f"合并候选池 {rep['merge']['n_candidates']} → 资格池 {rep['eligibility']['n_pool']}"
             f"（通用骨架去重剔除 {rep['eligibility']['n_deduped']}） → 名册 "
             f"{rep['roster']['n']} 行 = 20 候选 + 2 共用对照",
             fontsize=10.5, color=INK, ha="left", va="top")
    fig.text(0.015, 0.898,
             "字典序最优（frozen_rule v1.7.0 step6_md_batch）：① 分层覆盖下限 t → ② CNN_VS 跨度 "
             "→ ③ Σ CNN_VS → ④ 规范并列打破；四阶段 status=0、mip_gap=0",
             fontsize=9, color=MUTED, ha="left", va="top")
    fig.text(0.015, 0.850,
             f"边界：跨来源 CNN_VS 不可直接比（目录池是 top-4,990 截断，生成池是全范围）；"
             f"敏感性口径（来源池内百分位）与主名册重叠 {sens['overlap_with_main']}/20 席。"
             f"跨度与 Σ 的绝对量程含合并带来的量程差异，不作化学优劣解读。",
             fontsize=8.6, color=WARN_TX, ha="left", va="top")
    plots.save(fig, FIG / "fig_final_md_batch.png", dpi=200)
    return {"t_star": rep["optimality"]["t_star"], "span": rep["optimality"]["span"],
            "sum": rep["optimality"]["sum_CNN_VS"],
            "sensitivity_overlap": sens["overlap_with_main"], "source_counts": src_counts}


def main() -> int:
    font = plots.setup_cjk()
    roster = pd.read_csv(OUT / "tables/final_md_20.csv")
    rep = read_json(OUT / "reports/final_md.json")
    print(f"[18_final_md_figures] 字体 {font} · 名册 {len(roster)} 行")
    s1 = fig_cards(roster)
    print(f"[18_final_md_figures] fig_final_md_cards：{s1}")
    s2 = fig_batch(roster, rep)
    print(f"[18_final_md_figures] fig_final_md_batch：{s2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
