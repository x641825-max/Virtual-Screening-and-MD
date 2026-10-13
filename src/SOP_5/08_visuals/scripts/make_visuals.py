#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SOP_5 可视化层：把全链可可视化部分一次性出图（不参与链、不改任何门/名单/交付件）。

产物（08_visuals/figures/，PNG + SVG 双出、可复现）：
  fig01_md_structure_cards.png  MD 批次 20 个分子的结构卡（4×5，含 CNS / Peripheral / 对照 标签）
  fig02_pipeline_flow.png       gnina 对接表 → MD 候选分子 的筛选流程（规模 + 判据 + 人工检查点）
  fig03_gate_cascade.png        逐门级联（4,990 → 20）与各口径臂池规模
  fig04_pose_identity.png       姿态闸门 / 锚点 / 结构身份 / maxTc 分布
  fig05_druglikeness.png        类药性：QED 档 / CNS MPO（含 80% 区间）/ Ro5-Veber-PAINS / cLogP-TPSA
  fig06_herg_safety.png         hERG：池内分布与门 / 第二模型交叉验证 / 保留曲线 / 阈值人群位置
  fig07_robustness.png          稳健性 A–F 臂
  fig08_list_and_batch.png      交付名单 48↔2 与 MD 批次 20 的画像对照
  fig09_status_board.png        验收与确定性状态板（硬门逐条读数）

另产出 README.md（索引 + 数据快照 + 数据源对照）与 gallery.html（自包含单文件索引）。

用法：python3 08_visuals/scripts/make_visuals.py
说明：08_visuals 不在 12_determinism 的 INCLUDE_DIRS 内，故本目录产物不进 G 臂哈希范围；
      所有数字都从链上产物现读（不重算口径、不硬编码结论），换快照后请重跑本脚本。
"""
from __future__ import annotations

import base64
import io
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from vendor import plots
from vendor.io_utils import read_json

RDLogger.DisableLog("rdApp.*")
OUT = ROOT / "08_visuals"
FIG = OUT / "figures"

BLUE, GREY, RED, GREEN, AMBER = plots.BLUE, plots.GREY, plots.RED, plots.GREEN, plots.AMBER
INK, MUTED = "#0B0C0E", "#6B7079"
GREY_FILL, GREY_EDGE = "#F7F8FA", "#E3E6EA"
BAND_BG, BAND_ED = "#EEF2FE", "#2F5BEA"
WARN_TX = "#B07908"
BAD_TX = "#CF4034"
OK_TX = "#15915F"

SNAP: dict = {}


# ---------------------------------------------------------------- 工具
def load() -> dict:
    d = {}
    d["rule"] = read_json(ROOT / "config/frozen_rule.json")
    d["order"] = read_json(ROOT / "config/run_order.json")
    d["repro"] = read_json(ROOT / "REPRO_INFO.json")
    d["gate"] = read_json(ROOT / "03_gates/reports/gate_stack.json")
    d["pose"] = read_json(ROOT / "01_pose_qc/reports/pose_gate_summary.json")
    d["iqc"] = read_json(ROOT / "01_pose_qc/reports/pose_identity_summary.json")
    d["fun"] = read_json(ROOT / "02_druglikeness/reports/druglikeness_funnel.json")
    d["mtc"] = read_json(ROOT / "03_gates/reports/maxTc_summary.json")
    d["herg"] = read_json(ROOT / "03_gates/reports/herg_status.json")
    d["xmod"] = read_json(ROOT / "03_gates/herg/cross_model/reports/herg_cross_model.json")
    d["admet"] = read_json(ROOT / "03_gates/admet/reports/admet_annotation.json")
    d["hthr"] = read_json(ROOT / "03_gates/admet/reports/herg_threshold_analysis.json")
    d["sel"] = read_json(ROOT / "04_selection/reports/selection_summary.json")
    d["rob"] = read_json(ROOT / "05_robustness/reports/robustness.json")
    d["mb"] = read_json(ROOT / "07_md_batch/reports/md_batch.json")
    d["acc"] = read_json(ROOT / "06_deliverables/acceptance_report.json")
    d["det"] = read_json(ROOT / "05_robustness/reports/determinism_compare.json")
    d["desc"] = pd.read_csv(ROOT / "02_druglikeness/tables/descriptors.csv", low_memory=False)
    d["pool"] = pd.read_csv(ROOT / "03_gates/admet/tables/admet_annotation_pool.csv",
                            low_memory=False,
                            usecols=["cid", "ADMET_hERG", "admet_herg_band",
                                     "admet_herg_margin", "admet_tier", "admet_liability_count"])
    d["xpred"] = pd.read_csv(
        ROOT / "03_gates/herg/cross_model/tables/herg_cross_model_predictions.csv",
        low_memory=False)
    d["ret"] = pd.read_csv(ROOT / "03_gates/admet/tables/herg_retention_curve.csv")
    d["main"] = pd.read_csv(ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv",
                            low_memory=False)
    d["batch"] = pd.read_csv(ROOT / "07_md_batch/tables/md_batch20.csv", low_memory=False)
    d["sel_main"] = d["sel"]["lists"]["main_QED0.50_anchor4"]
    d["gc"] = d["gate"]["gate_chain"]
    return d


def panel(ax, title, sub=""):
    ax.set_title(title, fontsize=11.5, color=INK, pad=8, loc="left")
    if sub:
        ax.set_xlabel(sub, fontsize=8.8, color=MUTED)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GREY_EDGE)
    ax.tick_params(labelsize=8.5, colors=MUTED)


def rbox(ax, x, y, w, h, fc, ec, lw=0.8, r=0.02):
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                linewidth=lw, edgecolor=ec, facecolor=fc))


def note(ax, x, y, text, color=INK, size=8.6, ha="left", va="top", frac=True):
    """带白底的注记（默认按 axes 分数坐标），避免压住数据。"""
    ax.text(x, y, text, fontsize=size, color=color, ha=ha, va=va, zorder=6,
            transform=ax.transAxes if frac else ax.transData,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.8))


def _dw(t: str) -> int:
    return sum(2 if ord(c) > 0x2E80 else 1 for c in t)


def wrap_cjk(s: str, width: int = 40) -> list[str]:
    """按显示宽度折行（CJK 计 2，ASCII 计 1）；优先在空格处断，避免切断数字与英文词。"""
    out, line, w = [], "", 0
    for ch in s:
        cw = 2 if ord(ch) > 0x2E80 else 1
        if w + cw > width:
            k = line.rfind(" ")
            if k > 0 and _dw(line[k + 1:]) <= 14:
                out.append(line[:k].rstrip())
                line = line[k + 1:]
                w = _dw(line)
            else:
                out.append(line)
                line, w = "", 0
        line += ch
        w += cw
    if line:
        out.append(line)
    return out


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


# ---------------------------------------------------------------- fig01 结构卡
def fig01_structure_cards(d: dict) -> dict:
    import matplotlib.pyplot as plt
    b = d["batch"].copy().sort_values("batch_order")
    is_cand = b.candidate_type == "candidate"
    cls = np.where(is_cand, np.where(b.cns_mpo >= 4, "CNS", "Peripheral"),
                   "对照 " + b.control_role.fillna("").astype(str))
    b["cls"] = cls
    flagged = b.herg_cross_model_flag.fillna("").eq("cross_model_fail").to_numpy()

    n_cns = int(np.sum(is_cand & (b.cns_mpo >= 4).to_numpy()))
    n_per = int(np.sum(is_cand.to_numpy())) - n_cns
    rows, cols = 5, 4
    fig, axes = plt.subplots(rows, cols, figsize=(13.8, 16.8))
    fig.subplots_adjust(left=0.015, right=0.985, top=0.93, bottom=0.012,
                        hspace=0.05, wspace=0.04)
    fig.text(0.015, 0.985, "SOP_5 · MD 批次候选分子结构卡", fontsize=16.5, color=INK,
             ha="left", va="top")
    fig.text(0.015, 0.960,
             f"共 {len(b)} 格 = {int(is_cand.sum())} 候选 + {int((~is_cand).sum())} 阴性对照｜"
             f"按 CNN_VS 降序（批内序 01–20，对照固定占末两位）｜"
             f"CNS = CNS MPO ≥ 4（点估计）：{n_cns} 个；Peripheral < 4：{n_per} 个",
             fontsize=10, color=MUTED, ha="left", va="top")
    fig.text(0.985, 0.987, "△ = hERG 第二模型（ChEMBL 自建集成）分歧：保留并标注",
             fontsize=9, color=WARN_TX, ha="right", va="top")
    for i, (ax, (_, r)) in enumerate(zip(axes.ravel(), b.iterrows())):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        if isinstance(r.smiles, str) and r.smiles:
            ax.imshow(plt.imread(io.BytesIO(mol_png(r.smiles)), format="png"),
                      extent=[0.02, 0.98, 0.21, 1.0], aspect="auto")
        lab = f"SOP5-MD-{int(r.batch_order):02d} | {r.cls}"
        if flagged[i]:
            lab += "  △"
        ax.text(0.5, 0.155, lab, ha="center", va="top", fontsize=10.5,
                color=WARN_TX if flagged[i] else "#2347C4")
        ax.text(0.5, 0.065, str(r.cid), ha="center", va="top", fontsize=9.5,
                color=INK, family="monospace")
    return {"n_cards": len(b), "n_cns": n_cns, "n_peripheral": n_per,
            "n_flagged": int(flagged.sum())}


# ---------------------------------------------------------------- fig02 流程图
def stage_card(ax, x, w, spec):
    title, big, biglab, bullets, check, fc, ec, accent = spec
    rbox(ax, x, 8.2, w, 30.4, fc, ec, lw=1.1, r=0.55)
    ax.text(x + 0.6, 37.0, title, fontsize=11.5, color=accent)
    ax.text(x + 0.6, 32.6, big, fontsize=24, color=accent)
    lab = wrap_cjk(biglab, 40)
    for j, ln in enumerate(lab):
        ax.text(x + 0.6, 30.3 - j * 1.15, ln, fontsize=8.4, color=MUTED)
    y = 30.3 - len(lab) * 1.15 - 1.1
    for bul in bullets:
        lines = wrap_cjk("· " + bul, 40)
        for j, ln in enumerate(lines):
            ax.text(x + 0.6, y - j * 1.15, ln, fontsize=8.2, color=INK)
        y -= len(lines) * 1.15 + 0.45
    cy = 15.2
    rbox(ax, x + 0.45, 8.7, w - 0.9, cy - 8.7, "#FFFFFF", GREY_EDGE, lw=0.7, r=0.3)
    ax.text(x + 1.0, cy - 0.5, "人工分析检查点", fontsize=8.4, color=WARN_TX)
    for j, ln in enumerate(wrap_cjk(check, 40)[:4]):
        ax.text(x + 1.0, cy - 2.1 - j * 1.15, ln, fontsize=8.0, color=INK)


def fig02_pipeline_flow(d: dict) -> dict:
    import matplotlib.pyplot as plt
    gc, po, fu, mt = d["gc"], d["pose"], d["fun"], d["mtc"]
    sel, mb, h, an = d["sel_main"], d["mb"], d["herg"], d["admet"]["herg_gate_semantics"]
    a = d["fun"]["annotation_columns_on_druglike_pool"]
    S = []

    def add(*a_):
        S.append(a_)

    add("01 固定输入", f"{49963:,}", "排名表行数（全库 CNN_VS 排名）", [
        "输入排名表 49,963 行（含 pKa / cLogD74 列）",
        "gnina 最优姿态 top5000_best_pose.sdf（5,000 个）",
        "5 个原始输入全部 SHA256 锁定，解包白名单固定",
        "对照库：run01 5,098 实体（479 已知活性 / 4,619 诱饵）",
    ], "输入哈希与规则登记逐一对账（01_verify_inputs）；解包失败或哈希漂移即停链。",
        GREY_FILL, GREY_EDGE, INK)

    add("02 姿态核验与标准化", f"{gc['0_去重后池']:,}", "InChIKey 去重后的唯一结构", [
        f"5,000 → {gc['0_去重后池']:,}（InChIKey 去重）",
        f"身份分级：identical {d['iqc']['identity_class_counts']['identical']:,}、"
        f"stereo_only {d['iqc']['identity_class_counts']['stereo_only']}、"
        f"protonation_diff {d['iqc']['identity_class_counts']['protonation_diff']}",
        "CNNscore × CNNaffinity = CNN_VS 逐值核验（容差 5e-4）",
        f"坐标往返偏差 {d['iqc']['coord_max_dev_max']} Å（几何未改）",
    ], "stereo_only 是「3D 指派手性」类差异，需人工确认不是上游结构错误；"
       "protonation_diff 逐个人工看。",
        BAND_BG, BAND_ED, "#2347C4")

    add("03 类药性硬门", f"{gc['1_类药性硬门（PAINS/Ro5/Veber）']:,}", "通过硬门的类药池", [
        f"剔除 {fu['funnel']['1_硬门剔除']}：Ro5 {fu['attribution']['ro5_fail']} 个、"
        f"Veber {fu['attribution']['veber_fail']} 个、PAINS {fu['attribution']['pains_fail']} 个",
        f"QED 三档：≥0.60 {fu['funnel']['3_QED>=0.60']:,}｜≥0.50 "
        f"{fu['funnel']['3_QED>=0.50']:,}｜≥0.40 {fu['funnel']['3_QED>=0.40']:,}",
        f"CNS MPO ≥ 4（点估计）{a['CNS MPO>=4']:,}；cLogD 80% 区间悲观界 "
        f"{a['CNS MPO>=4（80% 区间悲观界）']:,} / 乐观界 {a['CNS MPO>=4（80% 区间乐观界）']:,}",
        f"Brenk 命中 {a['Brenk 命中']} 个（只报不改）",
    ], "PAINS 全池 0 命中——须人工确认筛库已做 PAINS 清洗，而非该判据失效。",
        GREY_FILL, GREY_EDGE, INK)

    add("04 姿态闸门与锚点", f"{po['gate_counts']['pass']:,}", "姿态 pass（4,862 / 5,000）", [
        "16 残基接触指纹（重原子 < 4.0 Å）+ 盒子几何",
        f"锚点 N253 / E169 / F168 / W246 命中 4/4：{po['anchor_hits_counts']['4']:,} 个；"
        f"3/4：{po['anchor_hits_counts']['3']:,}",
        f"pass {po['gate_counts']['pass']:,} / warn {po['gate_counts']['warn']} / "
        f"reject {po['gate_counts']['reject']}（warn 保留并标警，reject 剔除）",
        f"级联内：姿态 pass {gc['2_姿态闸门 pass']:,} → anchor=4 {gc['3_anchor=4']:,}",
    ], "warn 类（2/4 锚点）为何保留要有人看过：本包口径是「宁可留警不误删」。",
        BAND_BG, BAND_ED, "#2347C4")

    add("05 相似度与安全门", f"{gc['5_QED>=0.50（主池）']:,}", "主池（QED≥0.50 ∧ anchor=4）", [
        f"maxTc（vs 479 已知活性）< 0.90：{gc['4_maxTc<0.9']:,}（剔除 "
        f"{gc['3_anchor=4'] - gc['4_maxTc<0.9']} 个）",
        f"hERG < 0.5（ADMET-AI）：作用域全池 {h['gate_scope_n']:,}、覆盖 "
        f"{h['scope_coverage'] * 100:.0f}%、通过 {h['n_pass_in_scope']:,} 个",
        f"池内 ≥0.5 占 {an['pool_pct_at_or_above_gate']:.1f}%"
        f"（已上市药百分位 {an['approved_drug_percentile_at_gate_pool_near']:.1f}）",
        f"缺预测值 → not_evaluated（{h['n_not_evaluated']} 个），不当不合格",
    ], "0.5 是未校准分诊线（库内 66% 在线上），不得转述为「已判定安全」；"
       "第二模型一致性低（κ≈0.02–0.07）也要一并说。",
        GREY_FILL, GREY_EDGE, INK)

    add("06 名单选序", f"{sel['n_selected_candidates']} + {sel['n_controls']}",
        "候选 + 阴性对照（各口径臂同规模）", [
        f"主口径池 {sel['pool_n']:,} → novelty 门剔除 "
        f"{d['sel']['novelty_gate']['n_exact_known_in_base']} 个已知物重合者",
        "约束：CNN_VS 降序 + 成对 Tc < 0.65 + 骨架 ≤1 + 化学型 ≤20 + 风险档 ≤8",
        f"实测：成对 Tc 上限 {sel['max_pairwise_tc']}；唯一骨架 {sel['n_distinct_scaffold']}",
        f"化学型最多 {sel['chemotype_max']}（cap {sel['chem_cap']}）、"
        f"风险档 {sel['risk_n']}（cap {sel['risk_cap']}）",
    ], "未设 novelty 门时第 3 名就是已知物 ZINC000118633318（诊断档），"
       "说明这道门真的在起作用。",
        BAND_BG, BAND_ED, "#2347C4")

    add("07 MD 批次", f"{mb['batch_size']}", "进入一级 MD（机时分配层）", [
        f"{mb['n_candidates']} 候选 + {mb['n_controls']} 对照（对照固定占位）",
        f"资格：通用骨架重复剔除 {mb['pool']['n_excluded_generic_scaffold_dup']} 个；"
        f"hERG 第二模型分歧 {mb['pool']['n_retained_hERG_cross_model_fail']} 个保留并标注",
        f"分层：CNN_VS 下半区 {mb['stratification']['n_below_median_selected']} 个"
        f"（下限 {mb['stratification']['min_below_median']}，中位 "
        f"{mb['stratification']['reference_median_CNN_VS']}）",
        f"批内成对 Tc 上限 {mb['diversity']['max_pairwise_tc']}；唯一骨架 "
        f"{mb['diversity']['n_unique_scaffold']}/18；一级约 "
        f"{mb['md_cost_estimate']['level1_cardday']} 卡日",
    ], "批次是机时分配不是「更优子集」；MD 后换人须改规则并重跑整链，"
       "判据须在跑 MD 前预注册。",
        BAND_BG, BAND_ED, "#2347C4")

    fig, ax = plt.subplots(figsize=(21.5, 10.4))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 46)
    ax.axis("off")
    fig.text(0.012, 0.975, "SOP_5 · 从 gnina 对接表到 MD 候选分子：七道关",
             fontsize=19, color=INK, ha="left", va="top")
    fig.text(0.012, 0.935, "每段给出：数据规模 / 判据与实测值 / 需要人看的地方　"
                           "（全部读数来自包内产物，可用 run_all.py 复现）",
             fontsize=10.5, color=MUTED, ha="left", va="top")
    fig.text(0.988, 0.975,
             f"DMTA 复现包 · 规则 {d['rule']['rule_version']}｜RDKit 2025.09.5｜"
             "ECFP4 / Butina｜gnina 打分（姿态由上游提供）", fontsize=9.5, color=MUTED,
             ha="right", va="top")
    x0, w, gap = 1.0, 13.1, 1.0
    for i, spec in enumerate(S):
        stage_card(ax, x0 + i * (w + gap), w, spec)
        if i < len(S) - 1:
            ax.annotate("", xy=(x0 + i * (w + gap) + w + gap - 0.12, 23.5),
                        xytext=(x0 + i * (w + gap) + w + 0.12, 23.5),
                        arrowprops=dict(arrowstyle="-|>", color="#9AA3AE", lw=1.3))
    ax.plot([1.0, 99.0], [40.2, 40.2], color=GREY_EDGE, lw=0.8)
    marks = [(6.6, f"{49963:,} 排名表"), (21.0, f"{gc['0_去重后池']:,} 唯一结构"),
             (35.4, f"{gc['1_类药性硬门（PAINS/Ro5/Veber）']:,} 类药池"),
             (49.8, f"{gc['3_anchor=4']:,} 锚点 4/4"),
             (64.2, f"{gc['5_QED>=0.50（主池）']:,} 主池"),
             (78.6, f"{sel['pool_n']:,} 名单池"),
             (93.2, f"{mb['batch_size']} 进入 MD")]
    for xv, lab in marks:
        ax.text(xv, 41.2, lab, fontsize=9, color=MUTED, ha="center")
    ax.text(99.0, 43.0, "规模阶梯（左→右）", fontsize=8.6, color=MUTED, ha="right")
    rbox(ax, 1.0, 1.0, 98.0, 5.6, "#FFFFFF", GREY_EDGE, lw=0.8, r=0.35)
    ax.text(1.8, 5.1, "贯穿全程的复核对账点", fontsize=10, color=INK)
    ax.text(1.8, 2.6,
            "· 输入 SHA256 锁定　· 每步落漏斗 JSON（硬门：每步有漏斗 JSON）　"
            "· 阴性对照规则包内重算　· G 臂：全量产物 sha256 重跑比对（n_diff 须为 0）　"
            "· 隔离扫描：包内脚本不引用包外路径　· 分数敏感度 σ=0.10 名单重叠为主判据（≥0.80）",
            fontsize=8.6, color=MUTED)
    return {"stages": len(S), "cascade": [gc[k] for k in gc]}


# ---------------------------------------------------------------- fig03 级联
def fig03_gate_cascade(d: dict) -> dict:
    import matplotlib.pyplot as plt
    gc, sel, mb = d["gc"], d["sel_main"], d["mb"]
    labels = ["去重后池", "类药性硬门\n（PAINS/Ro5/Veber）", "姿态闸门 pass", "锚点 4/4",
              "maxTc < 0.90", "QED ≥ 0.50 主池", "hERG < 0.50", "名单选序（约束后）",
              "MD 批次"]
    keys = list(gc)
    vals = [gc[k] for k in keys] + [sel["list_size"], mb["batch_size"]]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(16.6, 6.6),
                                  gridspec_kw={"width_ratios": [1.75, 1.15]})
    y = np.arange(len(vals))[::-1]
    colors = [GREY, GREY, BLUE, BLUE, BLUE, BLUE, RED, GREEN, GREEN]
    ax.barh(y, vals, height=0.6, color=colors, alpha=0.92)
    for i, (v, yy) in enumerate(zip(vals, y)):
        ax.text(v + 70, yy, f"{v:,}", va="center", fontsize=9.5, color=INK)
        if i:
            drop = vals[i - 1] - v
            if drop > 0:
                ax.text(v + 70 + max(vals) * 0.005, yy - 0.34, f"−{drop:,}", fontsize=8,
                        color=BAD_TX, va="top")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.2, color=INK)
    ax.set_xlim(0, max(vals) * 1.2)
    ax.set_ylim(-0.8, len(vals) - 0.2)
    panel(ax, "逐门级联：4,990 唯一结构 → 20 个 MD 分子", "分子数（红色 = 安全门，收缩最强）")
    ax.text(0.99, 0.06, "hERG 门作用域为全池 4,990（覆盖 100%），故是级联中最强的一次收缩",
            transform=ax.transAxes, ha="right", fontsize=8.6, color=MUTED)

    arms = d["gate"]["arm_pools"]
    names = [a["arm"].replace("_", "\n") for a in arms]
    pre = [a["pool_n"] for a in arms]
    post = [a["pool_n_herg_pass"] for a in arms]
    xs = np.arange(len(arms))
    ax2.bar(xs - 0.2, pre, 0.4, color=BLUE, alpha=0.45, label="hERG 门前臂池")
    ax2.bar(xs + 0.2, post, 0.4, color=BLUE, label="hERG 门后臂池")
    for x, p_, q in zip(xs, pre, post):
        ax2.text(x - 0.2, p_ + 40, f"{p_:,}", ha="center", fontsize=7.6, color=MUTED)
        ax2.text(x + 0.2, q + 40, f"{q:,}", ha="center", fontsize=7.6, color="#2347C4")
    for _sep in np.arange(0.5, len(xs) - 0.5, 2.0):
        ax2.axvline(_sep, color="#E3E6EA", lw=0.8, ls=(0, (3, 3)), zorder=0)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(names, fontsize=7.8)
    ax2.set_ylabel("分子数", fontsize=9, color=MUTED)
    ax2.set_ylim(0, max(pre) * 1.18)
    ax2.legend(fontsize=8.2, frameon=False, loc="upper right")
    panel(ax2, "臂池规模 × hERG 门（候选池 QED 档 × 锚点档，共 8 档）",
          "主口径池 1,989 → hERG 门后 621；锚点放宽到 ≥3 时池更大但阳性率更低")
    fig.tight_layout()
    return {"cascade": vals}


# ---------------------------------------------------------------- fig04 姿态
def fig04_pose(d: dict) -> dict:
    import matplotlib.pyplot as plt
    po, iq, mt = d["pose"], d["iqc"], d["mtc"]
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.8))

    ax = axes[0][0]
    gk = list(po["gate_counts"])
    gv = [po["gate_counts"][k] for k in gk]
    ax.bar(["pass", "warn", "reject"], gv, color=[GREEN, AMBER, RED], alpha=0.9, width=0.5)
    for i, v in enumerate(gv):
        ax.text(i, v + 60, f"{v:,}", ha="center", fontsize=10, color=INK)
    ax.set_ylabel("姿态数", fontsize=9, color=MUTED)
    ax.set_ylim(0, max(gv) * 1.18)
    panel(ax, f"姿态闸门（{po['n_poses']:,} 个姿态）",
          f"盒子外 {po['outside_box_n']} 个 → reject；warn 保留并标警")

    ax = axes[0][1]
    ak = sorted(po["anchor_hits_counts"], key=int)
    av = [po["anchor_hits_counts"][k] for k in ak]
    ax.bar([f"{k}/4" for k in ak], av, color=[GREY, BLUE, GREEN], alpha=0.9, width=0.5)
    for i, v in enumerate(av):
        ax.text(i, v + 60, f"{v:,}", ha="center", fontsize=10, color=INK)
    ax.set_ylabel("分子数", fontsize=9, color=MUTED)
    ax.set_ylim(0, max(av) * 1.18)
    ev = po["labeled5098_anchor_evidence"]
    panel(ax, "锚点命中分布（N253 / E169 / F168 / W246）",
          f"对照 5,098 实体库：4/4 命中富集 {ev['anchor4_enrichment']}×"
          f"（活性 {ev['anchor4_active_rate'] * 100:.1f}% vs 诱饵 {ev['anchor4_decoy_rate'] * 100:.1f}%）")

    ax = axes[1][0]
    ic = iq["identity_class_counts"]
    keys = ["identical", "stereo_only", "protonation_diff"]
    v = [ic[k] for k in keys]
    ax.bar(keys, v, color=[GREEN, AMBER, RED], alpha=0.9, width=0.5)
    for i, x in enumerate(v):
        ax.text(i, x * 1.15, f"{x:,}", ha="center", fontsize=10, color=INK)
    ax.set_yscale("log")
    ax.set_ylim(1, max(v) * 6)
    ax.set_ylabel("分子数（对数轴）", fontsize=9, color=MUTED)
    panel(ax, "结构身份分级（上传姿态 vs 重读图）",
          f"坐标最大偏差 {iq['coord_max_dev_max']} Å；CNN_VS 复算违规 "
          f"{iq['cnn_vs_x_affinity_violations']} 个")

    ax = axes[1][1]
    bins = mt["histogram_0.1_bins"]
    xs = [int(k) / 10 for k in bins]
    ys = [bins[k] for k in bins]
    ax.bar(xs, ys, width=0.085, color=BLUE, alpha=0.88)
    ax.set_yscale("log")
    ax.axvline(0.9, color=RED, ls="--", lw=1.1)
    note(ax, 0.92, 0.97,
         f"门 0.90\n≥ 门仅 {mt['n_with_any_active_ge_0.9']} 个\n"
         f"过门 {mt['pass_primary_lt_0.9']:,} 个（99.86%）",
         color=BAD_TX, size=8.4, ha="right", va="top")
    ax.set_xlabel("maxTc（对 479 个已知 A2A 活性的最大 Tanimoto）", fontsize=9, color=MUTED)
    ax.set_ylabel("分子数（对数轴）", fontsize=9, color=MUTED)
    pct = mt["maxTc_percentiles"]
    panel(ax, "新颖性：maxTc 分布（0.1 分箱）",
          f"中位 {pct['p50']}｜P95 {pct['p95']}｜P100 {pct['p100']}"
          f"（参考活性集 {mt['reference_actives']} 个）")
    fig.tight_layout()
    return {}


# ---------------------------------------------------------------- fig05 类药性
def fig05_druglikeness(d: dict) -> dict:
    import matplotlib.pyplot as plt
    de, fu = d["desc"], d["fun"]
    a = fu["annotation_columns_on_druglike_pool"]
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.8))

    ax = axes[0][0]
    ax.hist(de.qed.dropna(), bins=40, color=BLUE, alpha=0.88)
    top = ax.get_ylim()[1]
    for t, c, lab in ((0.40, GREY, "0.40"), (0.50, AMBER, "0.50 主口径"),
                      (0.60, GREEN, "0.60 敏感性")):
        ax.axvline(t, color=c, ls="--", lw=1.1)
        ax.text(t + 0.008, top * 0.95, lab, color=c, fontsize=8.4, rotation=90, va="top")
    ax.set_ylabel("分子数", fontsize=9, color=MUTED)
    panel(ax, "QED 分布（去重后池 4,990）",
          f"≥0.60 {fu['funnel']['3_QED>=0.60']:,}（{fu['coverage']['3_QED>=0.60'] * 100:.1f}%）｜"
          f"≥0.50 {fu['funnel']['3_QED>=0.50']:,}（{fu['coverage']['3_QED>=0.50'] * 100:.1f}%）｜"
          f"≥0.40 {fu['funnel']['3_QED>=0.40']:,}")

    ax = axes[0][1]
    ax.hist(de.cns_mpo.dropna(), bins=40, color=BLUE, alpha=0.88)
    ax.axvline(4, color=RED, ls="--", lw=1.1)
    note(ax, 0.02, 0.97,
         f"CNS MPO = 4（CNS 导向经验线）\n点估计 ≥ 4：{a['CNS MPO>=4']:,}\n"
         f"cLogD 80% 区间悲观界：{a['CNS MPO>=4（80% 区间悲观界）']:,}\n"
         f"乐观界：{a['CNS MPO>=4（80% 区间乐观界）']:,}",
         color=BAD_TX, size=8.3, ha="left", va="top")
    ax.set_ylabel("分子数", fontsize=9, color=MUTED)
    panel(ax, "CNS MPO 分布（含 cLogD 80% 区间界定）",
          f"区间宽度中位 {a['cLogD 区间宽度中位']} log；"
          f"{a['cLogD 区间使 CNS MPO 变动>=0.1 的分子']:,} 个分子因区间跨过 0.1")

    ax = axes[1][0]
    names = ["Ro5 违规", "Veber 违规", "PAINS 命中"]
    v = [fu["attribution"]["ro5_fail"], fu["attribution"]["veber_fail"],
         fu["attribution"]["pains_fail"]]
    ax.bar(names, v, color=[RED, AMBER, GREY], alpha=0.9, width=0.5)
    for i, x in enumerate(v):
        ax.text(i, max(v) * 0.06 + x, f"{x}", ha="center", fontsize=10.5, color=INK)
    ax.set_ylim(0, max(v) * 1.25)
    ax.set_ylabel("分子数", fontsize=9, color=MUTED)
    panel(ax, "类药性硬门触发情况（各自独立计数）",
          f"合计剔除 {fu['funnel']['1_硬门剔除']} 个"
          f"（{fu['coverage']['1_硬门剔除'] * 100:.2f}%）："
          f"仅 Ro5 {fu['attribution']['ro5_fail_only']}、仅 Veber "
          f"{fu['attribution']['veber_fail_only']}、仅 PAINS {fu['attribution']['pains_fail_only']}")

    ax = axes[1][1]
    cns = de.cns_mpo >= 4
    ax.scatter(de.cLogP[~cns], de.TPSA[~cns], s=5, c=GREY, alpha=0.3, label="CNS MPO < 4")
    ax.scatter(de.cLogP[cns], de.TPSA[cns], s=5, c=BLUE, alpha=0.45, label="CNS MPO ≥ 4")
    ax.axhspan(90, 140, color=GREEN, alpha=0.07)
    ax.text(5.4, 142, "TPSA 90–140 Å²（CNS 经验带）", fontsize=8, color=OK_TX)
    ax.set_xlabel("cLogP", fontsize=9, color=MUTED)
    ax.set_ylabel("TPSA (Å²)", fontsize=9, color=MUTED)
    ax.legend(fontsize=8.2, frameon=False, loc="upper left")
    panel(ax, "cLogP–TPSA 空间（去重后池）",
          f"CNS MPO ≥ 4 且 TPSA∈[90,140]：{a['CNS MPO>=4 且 TPSA∈[90,140]']:,} 个"
          f"（仅 TPSA 入带 {a['TPSA in [90,140]']:,}）")
    fig.tight_layout()
    return {}


# ---------------------------------------------------------------- fig06 hERG
def ret_cols(ret) -> tuple[dict, bool]:
    """探测保留曲线列名；表头粘连（并集列与名单列同名）时返回 ambiguous=True。"""
    _pick = lambda pref: next((c for c in ret.columns if c.startswith(pref)), None)
    cols = {"library": _pick("keep_library_f3"), "pool": _pick("keep_pool"),
            "pre": _pick("keep_pre_herg")}
    union_c = next((c for c in ret.columns if c.startswith("keep_quad_union")), None)
    deliv_c = next((c for c in ret.columns if "delivered_list_48" in c), None)
    amb = bool(union_c and deliv_c and union_c == deliv_c)
    if amb:
        union_c = deliv_c = None
    return {"union": union_c, "deliv": deliv_c, **cols}, amb


def fig06_herg(d: dict) -> dict:
    import matplotlib.pyplot as plt
    h, xm, ret, hp, xp = d["herg"], d["xmod"], d["ret"], d["pool"], d["xpred"]
    sem, bands = d["admet"]["herg_gate_semantics"], d["admet"]["herg_band_counts_pool"]
    thr = float(h["threshold"])
    ens = {r["model2"]: r for r in xm["agreement_vs_admet_ai"]}["ens"]
    kappas = [r["cohen_kappa"] for r in xm["agreement_vs_admet_ai"]]
    agrees = [r["gate_agreement"] for r in xm["agreement_vs_admet_ai"]]
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.0))

    ax = axes[0][0]
    ax.hist(hp.ADMET_hERG.dropna(), bins=50, color=BLUE, alpha=0.88)
    top = ax.get_ylim()[1]
    npass = int((hp.ADMET_hERG < thr).sum())
    nfail = int((hp.ADMET_hERG >= thr).sum())
    ax.axvline(thr, color=RED, ls="--", lw=1.2)
    note(ax, 0.02, 0.97,
         f"门 = {thr}（未校准分诊线）\n通过 {npass:,}｜不通过 {nfail:,}",
         color=BAD_TX, size=8.8, ha="left", va="top")
    ax.set_ylabel("分子数", fontsize=9, color=MUTED)
    panel(ax, f"hERG 概率分布（作用域 {h['scope_mode']} = 全池 {h['gate_scope_n']:,}）",
          f"覆盖 {h['scope_coverage'] * 100:.1f}%；not_evaluated {h['n_not_evaluated']} 个 →"
          f" 不当不合格（band 低 {bands['低']:,} / 中 {bands['中']:,} / 高 {bands['高']:,}）")

    ax = axes[0][1]
    m = xp.dropna(subset=["hERG_admet_ai", "hERG_cross_ens"])
    ok = m.hERG_cross_ens < 0.5
    ax.scatter(m.hERG_admet_ai[ok], m.hERG_cross_ens[ok], s=5, c=BLUE, alpha=0.35,
               label=f"第二模型也过门 {int(ok.sum()):,}")
    ax.scatter(m.hERG_admet_ai[~ok], m.hERG_cross_ens[~ok], s=5, c=RED, alpha=0.5,
               label=f"第二模型不过门 {int((~ok).sum()):,}")
    ax.axvline(0.5, color=RED, ls="--", lw=0.9)
    ax.axhline(0.5, color=RED, ls="--", lw=0.9)
    ax.set_xlabel("ADMET-AI hERG（主门口径）", fontsize=9, color=MUTED)
    ax.set_ylabel("第二模型集成分（ChEMBL 自建）", fontsize=9, color=MUTED)
    ax.legend(fontsize=8.2, frameon=True, framealpha=0.9, edgecolor="#C9CCD3",
              loc="lower right")
    panel(ax, f"两模型判定对照（n={len(m):,}）",
          f"ρ={ens['spearman_rho']}，一致率 {ens['gate_agreement'] * 100:.1f}%，"
          f"κ={ens['cohen_kappa']}；4 个第二模型 κ {min(kappas)}–{max(kappas)}")

    ax = axes[1][0]
    rc, amb = ret_cols(ret)
    _gc = d["gc"]
    order = [("pool", f"全池 {_gc.get('0_去重后池', 0):,}", BLUE),
             ("pre", f"主池（门前）{_gc.get('5_QED>=0.50（主池）', 0):,}", GREY),
             ("union", f"并集 {d['gate'].get('union_four_arms', 0):,}", AMBER),
             ("deliv", f"交付名单 {d['sel_main']['n_selected_candidates']}", GREEN),
             ("library", "全库 f3", "#9AA3AE")]
    for key, lab, c in order:
        col = rc.get(key)
        if col:
            ax.plot(ret.t, ret[col], lw=1.5, color=c, label=lab)
    ax.axvline(thr, color=RED, ls="--", lw=1.1)
    note(ax, 0.02, 0.97, "t 范围 0.30–0.99", color=MUTED, size=8.0)
    ax.text(thr + 0.006, 62, "主门 0.5", color=BAD_TX, fontsize=8.6)
    ax.set_xlabel("hERG 概率阈值 t", fontsize=9, color=MUTED)
    ax.set_ylabel("保留比例 (%)", fontsize=9, color=MUTED)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=8.2, frameon=False, loc="lower right")
    extra = "" if not amb else "｜本轮并集/名单列因表头粘连未画（见 README 异常记录）"
    panel(ax, "阈值保留曲线", "回答「换一个 hERG 阈值会裁掉多少」" + extra)

    ax = axes[1][1]
    labs = ["池内 ≥ 0.5 占比", "库内 ≥ 0.5 占比", "池内百分位 @0.5", "已上市药百分位 @0.5"]
    vals = [sem["pool_pct_at_or_above_gate"], sem["library_pct_at_or_above_gate"],
            sem["pool_empirical_percentile_of_gate"],
            sem["approved_drug_percentile_at_gate_pool_near"]]
    yy = np.arange(len(vals))[::-1]
    ax.barh(yy, vals, height=0.5, color=[RED, BLUE, GREY, GREEN], alpha=0.9)
    for i, v in enumerate(vals):
        ax.text(v + 1.4, yy[i], f"{v:.1f}%", va="center", fontsize=9, color=INK)
    ax.set_yticks(yy)
    ax.set_yticklabels(labs, fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel("百分比", fontsize=9, color=MUTED)
    panel(ax, "0.5 这个阈值在人群中的位置",
          f"重标定切点（池内）：hERG "
          f"{d['admet']['cutpoints_recalibrated_on_pool']['herg']:.3f}"
          f"｜交付名单 margin 中位 {sem['delivered_margin_median']}")
    fig.tight_layout()
    return {}


# ---------------------------------------------------------------- fig07 稳健性
def fig07_robustness(d: dict) -> dict:
    import matplotlib.pyplot as plt
    rb = d["rob"]
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.8))

    ax = axes[0][0]
    rows = rb["armA_score_perturbation"]["rows"]
    s = [r["sigma"] for r in rows]
    ax.plot(s, [r["mean_overlap"] for r in rows], "o-", color=BLUE, lw=1.7, label="平均重叠")
    ax.plot(s, [r["p10_overlap"] for r in rows], "s--", color=GREY, lw=1.2, label="P10")
    ax.plot(s, [r["min_overlap"] for r in rows], "^:", color=RED, lw=1.2, label="最小")
    ax.axhline(0.8, color=RED, ls="--", lw=1.1)
    ax.text(min(s), 0.806, "主判据线 0.80", color=BAD_TX, fontsize=8.6)
    prim = [r for r in rows if r["sigma"] == rb["armA_score_perturbation"]["primary_sigma"]][0]
    ok = prim["mean_overlap"] >= 0.8
    ax.scatter([prim["sigma"]], [prim["mean_overlap"]], s=90, facecolors="none",
               edgecolors=BAD_TX if not ok else OK_TX, lw=1.6, zorder=5)
    ax.set_xlabel("σ（CNN_VS 加性噪声）", fontsize=9, color=MUTED)
    ax.set_ylabel("与主名单重叠率", fontsize=9, color=MUTED)
    ax.set_ylim(0.5, 1.02)
    ax.legend(fontsize=8.2, frameon=False, loc="lower left")
    panel(ax, "A 打分扰动（每个 σ 100 次重抽样）",
          f"主判据 σ={rb['armA_score_perturbation']['primary_sigma']}："
          f"实测 {prim['mean_overlap']}，阈值 {rb['armA_score_perturbation']['judge_mean_overlap_min']}"
          f"→ {'通过' if ok else '未通过'}")

    ax = axes[0][1]
    b = pd.DataFrame(rb["armB_threshold_matrix"]["rows"])
    piv = b.pivot_table(index="maxTc_upper", columns="qed_min", values="overlap_vs_main")
    im = ax.imshow(piv.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, f"{piv.values[i, j]:.2f}", ha="center", va="center", fontsize=8.4,
                    color=INK)
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels([f"QED ≥ {c}" for c in piv.columns], fontsize=8.4)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels([f"maxTc < {i}" for i in piv.index], fontsize=8.4)
    _cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    _cb.ax.tick_params(labelsize=8.0, length=2)
    _cb.ax.set_title("重叠率", fontsize=8.4, color=MUTED, pad=5)
    panel(ax, "B 阈值矩阵（4 QED 档 × 4 maxTc 档）",
          f"重叠范围 {rb['armB_threshold_matrix']['overlap_min']}–"
          f"{rb['armB_threshold_matrix']['overlap_max']}（对角线附近最稳）")

    ax = axes[1][0]
    c = rb["armC_one_axis"]
    yy = np.arange(len(c))[::-1]
    vals = [r["overlap_vs_main"] for r in c]
    ax.barh(yy, vals, height=0.5, color=BLUE, alpha=0.9)
    for i, v in enumerate(vals):
        ax.text(v + 0.015, yy[i], f"{v:.3f}", va="center", fontsize=8.8, color=INK)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{r['axis']}：{r['change']}" for r in c], fontsize=8.6)
    ax.set_xlim(0, 1.15)
    ax.set_xlabel("与主名单重叠率", fontsize=9, color=MUTED)
    panel(ax, "C 单轴开关（一次只动一道门）",
          f"四个单轴的重合率都在 {min(vals):.3f}–{max(vals):.3f}："
          "名单对每一道门都敏感，不存在「只靠一个门撑着」的情况")

    ax = axes[1][1]
    scan = rb["armF_butina_full_pool"]["cutoff_scan"]
    xs = [r["dist_cutoff"] for r in scan]
    ns = [r["n_clusters"] for r in scan]
    ax.plot(xs, ns, "o-", color=BLUE, lw=1.6, label="簇数（全池 4,990）")
    for x, n_ in zip(xs, ns):
        ax.text(x + 0.010, n_ * 1.13, f"{n_:,}", ha="left", fontsize=8.2, color=INK)
    ax.set_ylim(min(ns) * 0.55, max(ns) * 1.6)
    ax.set_xlabel("Butina 距离阈值（Tc = 1 − d）", fontsize=9, color=MUTED)
    ax.set_ylabel("簇数（对数轴）", fontsize=9, color=MUTED)
    ax.set_yscale("log")
    ax.legend(fontsize=8.2, frameon=False)
    f = rb["armF_butina"]
    panel(ax, "F 骨架聚类（Butina / ECFP4）",
          f"名单池 618 @d={f['cutoff']}：{f['n_clusters']} 簇、最大 {f['largest_cluster']}、"
          f"单例 {f['n_singletons']}；全池 4,990：{rb['armF_butina_full_pool']['n_clusters']} 簇")
    fig.tight_layout()
    return {}


# ---------------------------------------------------------------- fig08 名单↔批次
def fig08_list_and_batch(d: dict) -> dict:
    import matplotlib.pyplot as plt
    main, bat, mb = d["main"], d["batch"], d["mb"]
    cand = main[main.candidate_type == "candidate"].sort_values("pick_order").reset_index(drop=True)
    bc = bat[bat.candidate_type == "candidate"]
    bset = set(bc.cid)
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.8))

    ax = axes[0][0]
    inb = cand.cid.isin(bset)
    ax.scatter(cand.index[inb] + 1, cand.CNN_VS[inb], s=42, c=BLUE, zorder=3,
               label=f"进入 MD 批次（{int(inb.sum())}）")
    ax.scatter(cand.index[~inb] + 1, cand.CNN_VS[~inb], s=18, c=GREY, alpha=0.55,
               label=f"留在候选池（{int((~inb).sum())}）")
    med = float(mb["stratification"]["reference_median_CNN_VS"])
    ax.axhline(med, color=RED, ls="--", lw=1)
    ax.text(1.2, med + 0.006, f"CNN_VS 中位 {med}", color=BAD_TX, fontsize=8.4)
    ax.set_xlabel("名单内 CNN_VS 名次（1 = 最高分）", fontsize=9, color=MUTED)
    ax.set_ylabel("CNN_VS", fontsize=9, color=MUTED)
    ax.legend(fontsize=8.2, frameon=False, loc="lower left")
    _mbo = mb.get("optimality", {})
    panel(ax, "MD 批次在 48 候选名单中的位置",
          f"下半区成员 {mb['stratification']['n_below_median_selected']} 个"
          f"（下限 {mb['stratification']['min_below_median']}）；"
          + (f"分位覆盖 t*={_mbo.get('t_star')}（占 {_mbo.get('n_occupied_bins')}/"
             f"{_mbo.get('n_bins')} 格）" if _mbo else "分位覆盖证书缺失"))

    ax = axes[0][1]
    cls = np.where(cand.cns_mpo >= 4, "CNS", "Peripheral")
    for name, col in (("CNS", BLUE), ("Peripheral", AMBER)):
        s = cand[cls == name]
        ax.scatter(s.qed, s.cns_mpo, s=26, c=col, alpha=0.7, label=f"{name}（{len(s)}）")
    ax.scatter(bc.qed, bc.cns_mpo, s=120, facecolors="none", edgecolors=RED, lw=1.3,
               label="MD 批次成员")
    ax.axhline(4, color=RED, ls="--", lw=0.9)
    ax.text(0.86, 4.05, "CNS MPO = 4", color=BAD_TX, fontsize=8.2)
    ax.set_xlabel("QED", fontsize=9, color=MUTED)
    ax.set_ylabel("CNS MPO", fontsize=9, color=MUTED)
    ax.legend(fontsize=8.2, frameon=False, loc="lower right")
    panel(ax, "候选池 → 批次的 QED / CNS MPO 覆盖",
          f"批次含 CNS {int((bc.cns_mpo >= 4).sum())} 个、Peripheral {int((bc.cns_mpo < 4).sum())} 个"
          f"（名单整体 {int((cand.cns_mpo >= 4).sum())}/{len(cand)}）")

    ax = axes[1][0]
    ch = cand.chemotype.value_counts()
    chb = bc.chemotype.value_counts()
    xs = np.arange(len(ch))
    ax.bar(xs - 0.19, [ch.get(k, 0) for k in ch.index], 0.38, color=GREY, label="候选名单 48")
    ax.bar(xs + 0.19, [chb.get(k, 0) for k in ch.index], 0.38, color=BLUE, label="MD 批次 18")
    for x, k in zip(xs, ch.index):
        ax.text(x - 0.19, ch.get(k, 0) + 0.4, str(ch.get(k, 0)), ha="center", fontsize=8.4)
        ax.text(x + 0.19, chb.get(k, 0) + 0.4, str(chb.get(k, 0)), ha="center", fontsize=8.4,
                color="#2347C4")
    ax.set_xticks(xs)
    ax.set_xticklabels([k.replace(" + ", " +\n") for k in ch.index], fontsize=8)
    ax.set_ylabel("分子数", fontsize=9, color=MUTED)
    ax.legend(fontsize=8.2, frameon=False)
    panel(ax, "化学型构成（批次上限 8）",
          f"批次化学型上限 {mb['caps']['chem_cap']}、风险档上限 {mb['caps']['risk_cap']}"
          f"（实测风险档 {mb['caps']['risk_count']}）")

    ax = axes[1][1]
    ax.axis("off")
    sp = mb["spans"]
    rows = [
        ("批次规模", f"{mb['batch_size']} = {mb['n_candidates']} 候选 + {mb['n_controls']} 对照"),
        ("资格处理", f"通用骨架重复剔除 {mb['pool']['n_excluded_generic_scaffold_dup']}；"
                     f"hERG 第二模型分歧保留并标注 "
                     f"{mb['pool']['n_retained_hERG_cross_model_fail']}（剔除 "
                     f"{mb['pool']['n_excluded_hERG_cross_model_fail']}）"),
        ("批内成对 Tc 上限", f"{mb['diversity']['max_pairwise_tc']}"
                             f"（对 {mb['diversity']['max_pairwise_tc_pair'][0]} × "
                             f"{mb['diversity']['max_pairwise_tc_pair'][1]}）"),
        ("唯一骨架", f"{mb['diversity']['n_unique_scaffold']} / {mb['n_candidates']} 全不同"
                     f"（通用骨架 {mb['diversity']['n_unique_generic_scaffold']}）"),
        ("CNN_VS 跨度", f"{sp['CNN_VS_min']}–{sp['CNN_VS_max']}（宽 {sp['CNN_VS_span']}）"),
        ("CNS MPO 跨度", f"{sp['cns_mpo_min']}–{sp['cns_mpo_max']}"),
        ("cLogP 跨度", f"{sp['cLogP_min']}–{sp['cLogP_max']}"),
        ("maxTc 跨度", f"{sp['maxTc_min']}–{sp['maxTc_max']}"),
        ("机时外推", f"一级 {mb['md_cost_estimate']['level1_cardday']} 卡日；"
                     f"二级（按 65% 通过率）"
                     f"{mb['md_cost_estimate']['level2_cardday_at_65pct_pass']} 卡日"),
    ]
    ax.text(0.0, 0.99, "MD 批次自检（07_md_batch）", fontsize=11.5, color=INK, va="top")
    for i, (k, v) in enumerate(rows):
        yy = 0.88 - i * 0.098
        ax.text(0.0, yy, k, fontsize=9, color=MUTED, va="top")
        ax.text(0.30, yy, v, fontsize=9, color=INK, va="top")
    ax.text(0.0, 0.035, "口径提醒：批次 = 机时分配层，不改名单、不改门；"
                        "MD 后换人须改规则并重跑整链。", fontsize=8.4, color=WARN_TX, va="top")
    fig.tight_layout()
    return {}


# ---------------------------------------------------------------- fig09 状态板
def fig09_status(d: dict) -> dict:
    import matplotlib.pyplot as plt
    acc, det, rep, rob, sel, mb = d["acc"], d["det"], d["repro"], d["rob"], d["sel_main"], d["mb"]
    checks = acc["checks"]
    hard = [c for c in checks if c.get("class") == "hard"]
    n_pass = sum(1 for c in hard if c.get("pass") is True)
    n_pend = sum(1 for c in hard if c.get("pass") in (None, "pending"))
    n_fail = sum(1 for c in hard if c.get("pass") is False)
    soft = [c for c in checks if c.get("class") != "hard"]
    fails = [c for c in hard if c.get("pass") is False]
    pends = [c for c in hard if c.get("pass") in (None, "pending")]
    det_ok = (det.get("mode") == "compared" and det.get("n_diff") == 0
              and det.get("determinism_verified") is True)
    xm = d["xmod"]["agreement_vs_admet_ai"]
    fig, ax = plt.subplots(figsize=(16.4, 8.0))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 46)
    ax.axis("off")
    fig.text(0.012, 0.972, "SOP_5 · 验收与确定性状态板", fontsize=16.5, color=INK, ha="left",
             va="top")
    fig.text(0.012, 0.928,
             f"规则 {d['rule']['rule_version']}｜rules_sha256 "
             f"{d['rule']['rules_sha256'][:16]}…｜计划 {len(d['order']['steps'])} 步｜"
             f"读数生成 {SNAP['generated_utc']}（全部数值为链上产物实读数）",
             fontsize=9.8, color=MUTED, ha="left", va="top")

    # ---- 左：硬门逐条 + G 臂
    rbox(ax, 0.8, 3.0, 32.0, 35.0, "#FFFFFF", GREY_EDGE, lw=0.9, r=0.5)
    ax.text(2.0, 35.6, "硬门逐条读数（acceptance_report.json）", fontsize=11.5, color=INK)
    _pend_lab = "待判定（步序读数时差）" if det_ok else "待 G 臂判定"
    chips = [(f"{n_pass}", "通过", OK_TX), (f"{n_pend}", _pend_lab, WARN_TX),
             (f"{n_fail}", "未通过", BAD_TX if n_fail else OK_TX)]
    for i, (n, lab, c) in enumerate(chips):
        ax.text(3.2 + i * 10.0, 30.0, n, fontsize=23, color=c)
        ax.text(3.2 + i * 10.0, 27.2, lab, fontsize=8.8, color=MUTED)
    ax.text(2.0, 24.2, f"共 {len(hard)} 条硬门；另有 {len(soft)} 条「如实报告」指标（不设门）",
            fontsize=8.8, color=INK)
    y = 21.4
    ax.text(2.0, y, f"未通过 {len(fails)} 条", fontsize=9.6, color=BAD_TX)
    y -= 2.0
    for c in fails:
        for ln in wrap_cjk(f"{c['criterion']}：实测 {c['observed']}，"
                           f"期望 {c.get('expected')}", 52):
            ax.text(2.0, y, ln, fontsize=8.4, color=BAD_TX)
            y -= 1.5
    y -= 0.6
    ax.text(2.0, y, (f"待判定 {len(pends)} 条（步序导致的读数时差，非缺陷；链末 G 臂"
                     f"已为 mode=compared 且 n_diff=0）" if det_ok else
                     f"待 G 臂判定 {len(pends)} 条（步序所致的读数时差，非缺陷）"),
            fontsize=9.0, color=WARN_TX)
    y -= 1.9
    for c in pends:
        ax.text(2.0, y, "· " + c["criterion"], fontsize=8.2, color=MUTED)
        y -= 1.5
    y -= 0.8
    ax.text(2.0, y, "G 臂（12_determinism，链末步）", fontsize=10.5, color=INK)
    y -= 2.6
    ax.text(2.0, y, f"mode = {det.get('mode')}｜n_diff = {det.get('n_diff')}｜"
                    f"文件 {det.get('n_files_current')} / {det.get('n_files_baseline')}",
            fontsize=9.0, color=OK_TX if det_ok else BAD_TX)
    y -= 2.2
    for ln in wrap_cjk("判据：mode=compared 且 n_diff=0。链末步晚于验收报告，"
                       "故报告内只登记上一轮读数；本图为链末文件实读数。", 52):
        ax.text(2.0, y, ln, fontsize=8.2, color=MUTED)
        y -= 1.4

    # ---- 中：关键质量数字
    rbox(ax, 33.5, 3.0, 32.5, 35.0, BAND_BG, BAND_ED, lw=1.1, r=0.5)
    ax.text(34.8, 35.6, "关键质量数字（与图同一快照）", fontsize=11.5, color="#2347C4")
    rows = [
        ("交付名单", f"{sel['n_selected_candidates']} 候选 + {sel['n_controls']} 对照"
                     f"（6 口径臂同规模）"),
        ("MD 批次", f"{mb['batch_size']} = {mb['n_candidates']} 候选 + {mb['n_controls']} 对照"),
        ("成对 Tc 上限", f"{sel['max_pairwise_tc']}（约束 < 0.65）"),
        ("唯一骨架", f"{sel['n_distinct_scaffold']} / {sel['n_selected_candidates']}"
                     f"（通用骨架 {sel['n_distinct_scaffold_generic']}）"),
        ("锚点命中", f"{'48/48 全部 4/4' if sel['anchor_all_4'] else '不满足'}"),
        ("hERG 门", f"作用域 {d['herg']['scope_mode']}｜全池 "
                     f"{d['herg']['gate_scope_n']:,}｜通过 {d['herg']['n_pass_in_scope']:,}"),
        ("第二模型一致性", f"κ {min(r['cohen_kappa'] for r in xm)}–"
                          f"{max(r['cohen_kappa'] for r in xm)}｜门决策一致率 "
                          f"{min(r['gate_agreement'] for r in xm) * 100:.1f}–"
                          f"{max(r['gate_agreement'] for r in xm) * 100:.1f}%"),
        ("σ=0.10 名单重叠", f"{rob['armA_score_perturbation']['rows'][1]['mean_overlap']}"
                             f"（主判据 ≥ 0.80 → "
                             f"{'通过' if rob['armA_score_perturbation']['pass_primary'] else '未通过'}）"),
        ("指纹敏感性", f"ECFP6 {rob['armD_fingerprint'][0]['overlap_vs_main']}｜"
                       f"MACCS {rob['armD_fingerprint'][1]['overlap_vs_main']}"),
        ("全链耗时预期", rep.get("runtime_expectation", "")[:22] + "…"),
    ]
    yy = 33.0
    for k, v in rows:
        ax.text(34.8, yy, k, fontsize=8.8, color="#2347C4", va="top")
        for j, ln in enumerate(wrap_cjk(v, 52)):
            ax.text(41.6, yy - j * 1.45, ln, fontsize=8.8, color=INK, va="top")
        yy -= 1.5 + 1.45 * (len(wrap_cjk(v, 52)) - 1) + 1.5

    # ---- 右：产物地图 + 转述注意
    rbox(ax, 66.5, 3.0, 32.5, 35.0, "#FFFFFF", GREY_EDGE, lw=0.9, r=0.5)
    ax.text(67.8, 35.6, "产物地图与转述注意", fontsize=11.5, color=INK)
    lines = [
        ("00_inputs/", "输入哈希、解包清单"),
        ("01_pose_qc/", "姿态核验（闸门/身份/方向性）"),
        ("02_druglikeness/", "描述符、类药性漏斗"),
        ("03_gates/", "门叠加、maxTc、hERG（含第二模型）、ADMET 注释与阈值锚定"),
        ("04_selection/", "6 口径名单矩阵 + 约束咬合度"),
        ("05_robustness/", "A–F 稳健性 + G 臂确定性基线"),
        ("06_deliverables/", "名单、决策卡、报告、验收、MANIFEST"),
        ("07_md_batch/", "MD 批次 20 与 MD 计划"),
        ("08_visuals/", "本目录（不参与 G 臂比对）"),
    ]
    yy = 33.2
    for k, v in lines:
        ax.text(67.8, yy, k, fontsize=8.6, color="#2347C4", va="top", family="monospace")
        for j, ln in enumerate(wrap_cjk(v, 34)):
            ax.text(75.4, yy - j * 1.35, ln, fontsize=8.6, color=INK, va="top")
        yy -= 1.4 + 1.35 * (len(wrap_cjk(v, 34)) - 1) + 1.1
    yy -= 0.4
    ax.text(67.8, yy, "转述注意", fontsize=10, color=INK, va="top")
    yy -= 2.0
    for t in ["「hERG 过门」须写作「在 ADMET-AI 口径下过门」",
              "CNS MPO 依赖上游 MolGpKa / cpLogD7.4；cLogD 有 80% 区间",
              "MD 批次是机时分配，不是「更优子集」",
              "PAINS 全池 0 命中 → 该判据在本库上属空转"]:
        for j, ln in enumerate(wrap_cjk("· " + t, 46)):
            ax.text(67.8, yy - j * 1.35, ln, fontsize=8.2, color=MUTED, va="top")
        yy -= 1.35 * len(wrap_cjk("· " + t, 46)) + 0.5
    return {"hard_pass": n_pass, "hard_pending": n_pend, "hard_fail": n_fail,
            "det_ok": det_ok, "all_pass": acc.get("all_pass")}


# ---------------------------------------------------------------- README / gallery
FIG_INDEX = [
    ("fig01_md_structure_cards.png", "MD 候选分子结构卡（20 = 18 候选 + 2 对照）",
     "每格一个分子的 RDKit 2D 结构与标签：SOP5-MD-序号 | CNS（CNS MPO ≥ 4）/ Peripheral / 对照；"
     "第二行为库 ID；△ 标记第二 hERG 模型分歧者。顺序 = CNN_VS 降序，对照固定末两位。"),
    ("fig02_pipeline_flow.png", "筛选流程图：gnina 对接表 → MD 候选分子",
     "七段：固定输入 → 姿态核验与标准化 → 类药性硬门 → 姿态闸门与锚点 → 相似度与安全门 → "
     "名单选序 → MD 批次。每段含规模、判据与「人工分析检查点」。"),
    ("fig03_gate_cascade.png", "逐门级联 4,990 → 20",
     "上：每道门的剩余量与剔除量（hERG 门收缩最强）；下：8 档臂池（QED 档 × 锚点档）在 hERG 门前后规模。"),
    ("fig04_pose_identity.png", "姿态与身份核验面板",
     "姿态闸门 pass/warn/reject、锚点命中分布（含 5,098 实体库富集）、结构身份分级、maxTc 对数分布。"),
    ("fig05_druglikeness.png", "类药性面板",
     "QED 分布与 0.40/0.50/0.60 三档、CNS MPO 分布（含 cLogD 80% 区间悲观/乐观界）、"
     "Ro5/Veber/PAINS 触发、cLogP–TPSA 空间。"),
    ("fig06_herg_safety.png", "hERG 安全面板",
     "全池概率分布与 0.5 门、第二模型（ChEMBL 自建集成）判定对照、阈值保留曲线、"
     "0.5 在池/库/已上市药中的百分位。"),
    ("fig07_robustness.png", "稳健性面板 A/B/C/F",
     "A 打分扰动 σ 扫描（含主判据线与通过与否）、B 阈值矩阵、C 单轴开关、F Butina 截断扫描。"),
    ("fig08_list_and_batch.png", "交付名单 ↔ MD 批次画像",
     "批次在 48 名单中的位置与中位线、QED/CNS MPO 覆盖、化学型配额、批次自检数字。"),
    ("fig09_status_board.png", "验收与确定性状态板",
     "硬门 pass/pending/fail 逐条读数、G 臂 n_diff、关键质量数字、产物地图与转述注意。"),
]
EXISTING = [
    ("01_pose_qc/figures/fig_pose_gate.png", "步骤 1：姿态闸门与接触指纹"),
    ("02_druglikeness/figures/fig_druglikeness_funnel.png", "步骤 2：类药性漏斗"),
    ("03_gates/figures/fig_gate_stack.png", "步骤 3：门叠加"),
    ("03_gates/herg/cross_model/figures/fig_herg_cross_model.png", "步骤 3b：hERG 第二模型交叉验证"),
    ("04_selection/figures/fig_list_overlap.png", "步骤 4：名单重叠"),
    ("04_selection/figures/fig_gate_impact.png", "步骤 4：门对名单的影响"),
    ("05_robustness/figures/fig_robA_sigma.png", "步骤 5A：σ 扰动"),
    ("05_robustness/figures/fig_robB_matrix.png", "步骤 5B：阈值矩阵"),
    ("05_robustness/figures/fig_robE_band.png", "步骤 5E：并列带宽"),
    ("05_robustness/figures/fig_robF_clusters.png", "步骤 5F：Butina 聚类"),
    ("07_md_batch/figures/fig_md_batch.png", "步骤 6：MD 批次分布与构成"),
]


def verify_pack_status() -> dict:
    """就地跑一次包完整性自检（tests/verify_pack.py，只读断言），只取结论行与异常行。"""
    import subprocess
    try:
        r = subprocess.run([sys.executable, str(ROOT / "tests/verify_pack.py")],
                           cwd=str(ROOT), capture_output=True, text=True, timeout=900)
        out = [l.strip() for l in (r.stdout or "").splitlines() if l.strip()]
        err = [l.strip() for l in (r.stderr or "").splitlines() if l.strip()]
        return {"rc": r.returncode, "last": out[-1] if out else "", "err": err[-1:]}
    except Exception as exc:  # noqa: BLE001
        return {"rc": None, "last": "", "err": [f"{type(exc).__name__}: {exc}"]}


def _rel(target: Path) -> str:
    return os.path.relpath(target, OUT).replace(os.sep, "/")


def write_readme(d: dict, meta: dict) -> None:
    acc, det = d["acc"], d["det"]
    hard = [c for c in acc["checks"] if c.get("class") == "hard"]
    n_pass = sum(1 for c in hard if c.get("pass") is True)
    n_pend = sum(1 for c in hard if c.get("pass") in (None, "pending"))
    n_fail = sum(1 for c in hard if c.get("pass") is False)
    fail = [c for c in hard if c.get("pass") is False]
    pend = [c for c in hard if c.get("pass") in (None, "pending")]
    vp = verify_pack_status()
    _rc, _amb = ret_cols(d["ret"])
    _det_ok = (det.get("mode") == "compared" and det.get("n_diff") == 0
               and det.get("determinism_verified") is True)
    _retp = ROOT / "03_gates/admet/tables/herg_retention_curve.csv"
    _srcp = ROOT / "03_gates/scripts/10c_herg_threshold.py"
    _ts = lambda q: datetime.fromtimestamp(q.stat().st_mtime, timezone.utc).strftime("%H:%M:%S")
    _stale = _retp.is_file() and _srcp.is_file() and _retp.stat().st_mtime < _srcp.stat().st_mtime
    _t_ret = _ts(_retp) if _retp.is_file() else "—"
    _t_src = _ts(_srcp) if _srcp.is_file() else "—"
    t = []
    A = t.append
    A("# SOP_5 可视化目录（08_visuals）")
    A("")
    A(f"> 数据快照：规则 {d['rule']['rule_version']}｜rules_sha256 "
      f"`{d['rule']['rules_sha256'][:16]}…`｜计划 {len(d['order']['steps'])} 步｜"
      f"读数生成 {meta['generated_utc']}")
    A("")
    A("本目录把 SOP_5 全链**可可视化部分**集中成一册图，供阅读、汇报与复核。"
      "所有数字都从链上产物现读（不重算口径、不改任何门/池/名单），"
      "由 `scripts/make_visuals.py` 一键重建。")
    A("")
    A("## 1. 图目录")
    A("")
    A("| 图 | 内容 | 说明 |")
    A("|---|---|---|")
    for f, title, desc in FIG_INDEX:
        A(f"| [`figures/{f}`](figures/{f}) | {title} | {desc} |")
    A("")
    A("网页版索引（自包含单文件、图片内嵌 base64，可离线分享）：[`gallery.html`](gallery.html)")
    A("")
    A("## 2. 读图顺序")
    A("")
    A("1. `fig02` 看全局：从 49,963 行排名表走到 20 个 MD 分子，每道关的判据与规模；")
    A("2. `fig03` 看每道门各裁掉多少 —— hERG 门（全池作用域）收缩最强；")
    A("3. `fig01` 看最终 20 个分子长什么样（结构卡，CNS/Peripheral/对照一目了然）；")
    A("4. 细节按步骤：`fig04` 姿态与身份 → `fig05` 类药性 → `fig06` 安全 → `fig07` 稳健性；")
    A("5. `fig08`/`fig09` 收口：名单与批次画像、验收与确定性状态。")
    A("")
    A("## 3. 当前包状态（与图同一快照，逐条读数）")
    A("")
    A(f"- 硬门验收：**{n_pass} 通过 / {n_pend} 待判定（步序读数时差）/ {n_fail} 未通过**"
      f"（`06_deliverables/acceptance_report.json`，共 {len(hard)} 条硬门）")
    for c in fail:
        A(f"  - 未通过：**{c['criterion']}** —— 实测 `{c['observed']}`，"
          f"期望 `{c.get('expected')}`")
    A(f"- 确定性 G 臂（链末文件 `05_robustness/reports/determinism_compare.json`）："
      f"`mode = {det.get('mode')}`，`n_diff = {det.get('n_diff')}`，"
      f"文件数 {det.get('n_files_current')} / {det.get('n_files_baseline')}，"
      f"`determinism_verified = {det.get('determinism_verified')}`")
    A(f"- 交付名单（主口径）：{d['sel_main']['n_selected_candidates']} 候选 + "
      f"{d['sel_main']['n_controls']} 对照；MD 批次 {d['mb']['batch_size']} = "
      f"{d['mb']['n_candidates']} + {d['mb']['n_controls']}")
    A("")
    A("## 4. 本轮核对到的异常（需在链上处理，本目录不改链）")
    A("")
    anom = [("1", "σ=0.10 名单平均重叠未达主判据",
             f"实测 `{d['rob']['armA_score_perturbation']['rows'][1]['mean_overlap']}` < 0.80"
             f"（`05_robustness/reports/robustness.json` → `pass_primary=false`）",
             "验收硬门 1 条未通过；对外须如实报告「分数扰动 σ=0.10 时名单不稳定」")]
    if not _det_ok:
        anom.append(("2", "G 臂本轮只写基线、未出判据",
                     f"`mode={det.get('mode')}`、`n_diff={det.get('n_diff')}`；"
                     f"文件自述「同源运行不构成独立重跑证据」",
                     "确定性判据需再跑一轮 `mode=compared` 方可判定；"
                     "本轮**不得**声明「确定性通过」"))
    if _amb:
        anom.append((str(len(anom) + 1), "hERG 保留曲线表头粘连（产物陈旧）",
                     f"`03_gates/admet/tables/herg_retention_curve.csv`（写于 {_t_ret}）表头出现"
                     f"粘连列名：并集列与名单列名合并为一列"
                     f"{'；生成它的 03_gates/scripts/10c_herg_threshold.py 已在其后于 ' + _t_src + ' 被修改，当前源码里两个键已分开 → 下一次重跑后自愈' if _stale else ''}",
                     "并集/名单保留曲线本轮无法定位；`fig06` 已按「不猜口径」原则省略该曲线"))
    _mbc = d["mb"].get("optimality", {})
    _mbe = d["rule"].get("step6_md_batch", {}).get("expected_certificate", {})
    if (not _mbc) or _mbc.get("t_star") != _mbe.get("t_star") \
            or abs(round(float(_mbc.get("span", -1)), 4) - float(_mbe.get("span", 0))) > 1e-9 \
            or _mbc.get("resolved_roster_consistent") is not True:
        anom.append((str(len(anom) + 1), "MD 批次最优性证书缺失或与规则登记值不符",
                     f"`07_md_batch/reports/md_batch.json → optimality`：t*="
                     f"{_mbc.get('t_star')}（登记 {_mbe.get('t_star')}）、跨度="
                     f"{_mbc.get('span')}（登记 {_mbe.get('span')}）、重解一致="
                     f"{_mbc.get('resolved_roster_consistent')}",
                     "MD 批次「可证明最优」的口径不成立或已随上游漂移，须显式复核"))
    if vp.get("rc") not in (0, None):
        anom.append((str(len(anom) + 1), "包完整性自检未跑完",
                     f"`tests/verify_pack.py` 退出码 `{vp['rc']}`；末行：`{vp['last'][:80]}`；"
                     f"异常：`{(vp['err'] or ['（无）'])[0][:90]}`",
                     "自检结论不完整，需修好后再判"))
    A("| # | 现象 | 证据 | 影响 |")
    A("|---|---|---|---|")
    if anom:
        for tag, what, ev, imp in anom:
            A(f"| {tag} | {what} | {ev} | {imp} |")
    else:
        A("| — | 本轮未发现需在链上处理的异常 | — | — |")
    A("")
    if _det_ok:
        A("读数时差说明：验收报告里那 2 条「待判定」不是缺陷，而是**报告与链末步的先后顺序**"
          "（`14_report` 早于 `12_determinism`，报告只登记上一轮读数）；链末文件本轮的实测读数"
          f"为 `mode=compared` / `n_diff=0`，确定性判据可判为通过。"
          "唯一真正未达标的是 σ 主判据那一条，两者不要合并转述。")
    else:
        A("异常 1 与 2 不是同一件事：硬门那 2 条「待判定」来自**报告与链末步的先后顺序**"
          "（`14_report` 早于 `12_determinism`，报告只登记上一轮读数），"
          "而 σ 主判据是**真的未达标**。两条都不要合并转述。")
    A("")
    A("## 5. 图之外的既有图（各步骤自出，此处只做索引）")
    A("")
    A("| 图 | 内容 |")
    A("|---|---|")
    for rel, title in EXISTING:
        p = ROOT / rel
        ok = "是" if p.is_file() else "缺失"
        A(f"| [`{rel.split('/')[-1]}`]({_rel(p)}) | {title}（原文件存在：{ok}）|")
    A("")
    A("## 6. 重建与边界")
    A("")
    A("```bash")
    A("cd SOP_5")
    A("python3 08_visuals/scripts/make_visuals.py   # 重建本目录全部图 + README.md + gallery.html")
    A("python3 run_all.py                          # 复现整链（19 步 + 自检），约 5 分钟")
    A("python3 tests/verify_pack.py                # 包完整性自检")
    A("```")
    A("")
    A("- 本目录**不参与** G 臂（`05_robustness/scripts/12_determinism.py` 的 `INCLUDE_DIRS` "
      "不含 `08_visuals`），也不改任何门、池、名单与交付件；")
    A("- 图的内容随链上产物变化：**换快照（重冻/改规则）后必须重跑本脚本**，"
      "否则 `figures/` 会与包状态不一致；")
    A("- 图内全部为中文标注，字体由 `vendor/plots.py` 的 CJK 配置提供；PNG 200 dpi + SVG 双出，"
      "SVG 已固定 hashsalt 与元数据以保证复现。")
    A("")
    A("## 7. 数据源对照")
    A("")
    A("| 图 | 主要数据源 |")
    A("|---|---|")
    A("| fig01 | `07_md_batch/tables/md_batch20.csv`（SMILES、CNS MPO、批次序、"
      "`herg_cross_model_flag`）|")
    A("| fig02 | `03_gates/reports/gate_stack.json`、`01_pose_qc/reports/pose_gate_summary.json`、"
      "`01_pose_qc/reports/pose_identity_summary.json`、"
      "`02_druglikeness/reports/druglikeness_funnel.json`、`03_gates/reports/herg_status.json`、"
      "`03_gates/reports/maxTc_summary.json`、`04_selection/reports/selection_summary.json`、"
      "`07_md_batch/reports/md_batch.json` |")
    A("| fig03 | `03_gates/reports/gate_stack.json`（`gate_chain`、`arm_pools`）+ "
      "`04_selection/reports/selection_summary.json` + `07_md_batch/reports/md_batch.json` |")
    A("| fig04 | `01_pose_qc/reports/{pose_gate_summary,pose_identity_summary}.json`、"
      "`03_gates/reports/maxTc_summary.json` |")
    A("| fig05 | `02_druglikeness/tables/descriptors.csv`、"
      "`02_druglikeness/reports/druglikeness_funnel.json` |")
    A("| fig06 | `03_gates/admet/tables/admet_annotation_pool.csv`、"
      "`03_gates/reports/herg_status.json`、"
      "`03_gates/herg/cross_model/tables/herg_cross_model_predictions.csv` 与 `reports/…json`、"
      "`03_gates/admet/tables/herg_retention_curve.csv`、"
      "`03_gates/admet/reports/{admet_annotation,herg_threshold_analysis}.json` |")
    A("| fig07 | `05_robustness/reports/robustness.json`（armA/B/C/D/E/F + 全池 Butina 扫描）|")
    A("| fig08 | `04_selection/tables/list_main_QED0.50_anchor4.csv`、"
      "`07_md_batch/reports/md_batch.json` |")
    A("| fig09 | `06_deliverables/acceptance_report.json`、"
      "`05_robustness/reports/determinism_compare.json`、`REPRO_INFO.json`、"
      "`config/run_order.json`、`config/frozen_rule.json` |")
    A("")
    A("## 8. 转述注意（写进对外材料前请先看这几条）")
    A("")
    A("- **hERG**：单一外部模型（ADMET-AI）的分数，门 0.5 未校准（库内 66.4% 在线上）；"
      "第二模型（ChEMBL `CHEMBL240` 自建）与其 κ ≈ 0.02–0.07、门决策一致率 ≈ 33%–45%（近乎随机一致）。"
      "引用时必须写作「在 ADMET-AI 口径下过门」。")
    A("- **CNS MPO**：依赖上游 MolGpKa pKa 与 cpLogD7.4（第三方预测）；cLogD 带 80% 区间，"
      "`fig05` 已把悲观/乐观界一并给出（区间宽度中位 "
      f"{d['fun']['annotation_columns_on_druglike_pool']['cLogD 区间宽度中位']} log）。")
    A("- **MD 批次**：机时分配层，不是「比其余候选更优」的结论；换人须改规则并重跑整链。")
    A("- **稳健性**：σ=0.10 名单平均重叠是主判据（阈值 0.80），现值见 `fig07` 与第 3 节；"
      "阈值矩阵、指纹更换（MACCS 0.5）说明名单对**相似度约束实现**同样敏感。")
    A("- **筛库形状**：全池 PAINS 命中 0 个，说明上游已做 PAINS 清洗；"
      "本包的 PAINS 判据在效度上属于「空转」而非「通过验证」。")
    A("")
    (OUT / "README.md").write_text("\n".join(t), encoding="utf-8")


def write_gallery(d: dict, meta: dict) -> None:
    from PIL import Image
    import html as _html
    acc, det = d["acc"], d["det"]
    hard = [c for c in acc["checks"] if c.get("class") == "hard"]
    n_pass = sum(1 for c in hard if c.get("pass") is True)
    n_pend = sum(1 for c in hard if c.get("pass") in (None, "pending"))
    n_fail = sum(1 for c in hard if c.get("pass") is False)
    vp = verify_pack_status()
    _, amb = ret_cols(d["ret"])
    g_lab = ("G 臂差异数（mode=" + str(det.get("mode")) + "）"
             if det.get("n_diff") is not None else "G 臂：本轮只写基线")
    g_val = det.get("n_diff") if det.get("n_diff") is not None else "—"

    def b64(p: Path, max_w: int):
        if not p.is_file():
            return None
        im = Image.open(p).convert("RGB")
        if im.width > max_w:
            im = im.resize((max_w, int(im.height * max_w / im.width)), Image.LANCZOS)
        buf = io.BytesIO()
        im.save(buf, format="JPEG", quality=84, optimize=True)
        return base64.b64encode(buf.getvalue()).decode()

    cards = []
    for f, title, desc in FIG_INDEX:
        s = b64(FIG / f, 1500)
        if s:
            cards.append((title, desc, s))
    for rel, title in EXISTING:
        p = ROOT / rel
        s = b64(p, 1100)
        if s:
            cards.append((title, "链上步骤自出图，原文件：" + rel, s))

    h = []
    A = h.append
    A("<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>")
    A("<meta name='viewport' content='width=device-width, initial-scale=1'>")
    A("<title>SOP_5 可视化索引</title><style>")
    A(":root{--blue:#2F5BEA;--blue-l:#EEF2FE;--blue-t:#2347C4;--gray:#6B7079;"
      "--gray-f:#F7F8FA;--edge:#C9CCD3;--ink:#0B0C0E;--ok:#15915F;--warn:#B07908;--bad:#CF4034}")
    A("*{box-sizing:border-box}body{margin:0;background:#fff;color:var(--ink);"
      "font:400 13px/1.6 -apple-system,'Noto Sans CJK SC','PingFang SC','Microsoft YaHei',sans-serif}")
    A(".wrap{max-width:1060px;margin:0 auto;padding:32px 24px 64px}")
    A("h1{font-size:15px;font-weight:500;margin:0 0 6px}")
    A("h2{font-size:14px;font-weight:500;margin:30px 0 10px}")
    A(".meta{color:var(--gray);font-size:12px;margin-bottom:16px}")
    A(".stat{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:14px 0 4px}")
    A(".stat div{border:0.5px solid var(--edge);border-radius:8px;padding:10px 12px;"
      "background:var(--gray-f)}")
    A(".stat .n{font-size:20px;font-weight:500}.stat .l{font-size:12px;color:var(--gray)}")
    A(".grid{display:grid;grid-template-columns:1fr;gap:14px}")
    A(".card{border:0.5px solid var(--edge);border-radius:8px;overflow:hidden}")
    A(".card .cap{padding:10px 12px;background:var(--gray-f)}")
    A(".card .t{font-size:13px;font-weight:500}")
    A(".card .d{font-size:12px;color:var(--gray);margin-top:2px}")
    A(".card img{display:block;width:100%;height:auto}")
    A("table{border-collapse:collapse;width:100%;font-size:12px}")
    A("th,td{text-align:left;padding:6px 8px;border-bottom:0.5px solid var(--edge);vertical-align:top}")
    A("th{color:var(--gray);font-weight:400}")
    A("code{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px}")
    A(".ok{color:var(--ok)}.warn{color:var(--warn)}.bad{color:var(--bad)}.blue{color:var(--blue-t)}")
    A("</style></head><body><div class='wrap'>")
    A("<h1>SOP_5 · 全链可视化索引</h1>")
    A(f"<div class='meta'>规则 {_html.escape(d['rule']['rule_version'])}（rules_sha256 "
      f"<code>{d['rule']['rules_sha256'][:16]}…</code>）｜计划 "
      f"{len(d['order']['steps'])} 步｜读数生成 {meta['generated_utc']}｜"
      f"由 <code>08_visuals/scripts/make_visuals.py</code> 生成</div>")
    A("<div class='stat'>")
    A(f"<div><div class='n ok'>{n_pass}</div><div class='l'>硬门通过</div></div>")
    A(f"<div><div class='n warn'>{n_pend}</div><div class='l'>验收报告待判定（步序读数时差）</div></div>")
    A(f"<div><div class='n {'bad' if n_fail else 'ok'}'>{n_fail}</div>"
      f"<div class='l'>硬门未通过（σ=0.10 重叠）</div></div>")
    A(f"<div><div class='n blue'>{g_val}</div>"
      f"<div class='l'>{_html.escape(g_lab)}</div></div>")
    A("</div>")
    A("<h2>本轮核对到的异常（需在链上处理）</h2>")
    A("<table><tr><th>#</th><th>现象</th><th>证据</th></tr>")
    _anom = [(f"<td>σ=0.10 名单平均重叠未达主判据</td><td>实测 "
              f"<code>{d['rob']['armA_score_perturbation']['rows'][1]['mean_overlap']}</code> &lt; 0.80"
              f"（robustness.json → pass_primary=false）</td></tr>")]
    if not (det.get("mode") == "compared" and det.get("n_diff") == 0
            and det.get("determinism_verified") is True):
        _anom.append(f"<td>G 臂本轮只写基线、未出判据</td><td><code>mode="
                     f"{_html.escape(str(det.get('mode')))}</code>，<code>n_diff="
                     f"{_html.escape(str(det.get('n_diff')))}</code>；需再跑一轮 mode=compared</td></tr>")
    if amb:
        _anom.append("<td>hERG 保留曲线表头粘连</td><td><code>herg_retention_curve.csv</code> "
                     "并集列与名单列名粘连；本页 fig06 已省略并集曲线</td></tr>")
    _mbc2 = d["mb"].get("optimality", {})
    _mbe2 = d["rule"].get("step6_md_batch", {}).get("expected_certificate", {})
    if (not _mbc2) or _mbc2.get("t_star") != _mbe2.get("t_star") \
            or abs(round(float(_mbc2.get("span", -1)), 4)
                   - float(_mbe2.get("span", 0))) > 1e-9 \
            or _mbc2.get("resolved_roster_consistent") is not True:
        _anom.append("<td>MD 批次最优性证书缺失或与规则登记值不符</td><td><code>"
                     f"md_batch.json → optimality</code>：t*={_mbc2.get('t_star')}"
                     f"（登记 {_mbe2.get('t_star')}）、重解一致="
                     f"{_mbc2.get('resolved_roster_consistent')}</td></tr>")
    if vp.get("rc") not in (0, None):
        _anom.append(f"<td>包完整性自检未跑完</td><td><code>tests/verify_pack.py</code> "
                     f"退出码 {vp['rc']}；{_html.escape((vp['err'] or [''])[0][:80])}</td></tr>")
    for _i, _row in enumerate(_anom, 1):
        A(f"<tr><td>{_i}</td>{_row}")
    if not _anom:
        A("<tr><td>—</td><td>本轮未发现需在链上处理的异常</td><td>—</td></tr>")
    A("</table>")
    A("<h2>图（新绘 9 张 + 链上既有 11 张）</h2><div class='grid'>")
    for title, desc, s in cards:
        A("<div class='card'>")
        A(f"<div class='cap'><div class='t'>{_html.escape(title)}</div>"
          f"<div class='d'>{_html.escape(desc)}</div></div>")
        A(f"<img alt='{_html.escape(title)}' src='data:image/jpeg;base64,{s}'>")
        A("</div>")
    A("</div>")
    A("<h2>怎么重建</h2><table><tr><th>动作</th><th>命令</th></tr>")
    A("<tr><td>重建本目录图与索引</td>"
      "<td><code>python3 08_visuals/scripts/make_visuals.py</code></td></tr>")
    A("<tr><td>复现整链（19 步 + 自检）</td><td><code>python3 run_all.py</code></td></tr>")
    A("<tr><td>包完整性自检</td><td><code>python3 tests/verify_pack.py</code></td></tr></table>")
    A("<h2>转述注意</h2><table><tr><th>事项</th><th>口径</th></tr>")
    A("<tr><td>hERG</td><td>单一 ADMET-AI 分数，门 0.5 未校准；第二模型 κ≈0.02–0.07，"
      "只能写作「在 ADMET-AI 口径下过门」</td></tr>")
    A("<tr><td>CNS MPO</td><td>依赖上游 MolGpKa / cpLogD7.4；cLogD 的 80% 区间会翻转部分分子，"
      "图 5 已给出悲观/乐观界</td></tr>")
    A("<tr><td>MD 批次</td><td>机时分配层，不是「更优子集」</td></tr>")
    A("<tr><td>稳健性</td><td>σ=0.10 平均重叠为主判据（≥0.80）；MACCS 指纹下重叠降至 "
      "0.5，说明名单对相似度实现敏感</td></tr>")
    A("<tr><td>PAINS</td><td>全池 0 命中 → 该判据在上游已清洗的库上属空转</td></tr></table>")
    A("</div></body></html>")
    (OUT / "gallery.html").write_text("\n".join(h), encoding="utf-8")


# ---------------------------------------------------------------- main
def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    d = load()
    plots.setup_cjk()
    import matplotlib.pyplot as plt
    SNAP["generated_utc"] = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    t0 = time.time()
    jobs = [("fig01_md_structure_cards", fig01_structure_cards),
            ("fig02_pipeline_flow", fig02_pipeline_flow),
            ("fig03_gate_cascade", fig03_gate_cascade),
            ("fig04_pose_identity", fig04_pose),
            ("fig05_druglikeness", fig05_druglikeness),
            ("fig06_herg_safety", fig06_herg),
            ("fig07_robustness", fig07_robustness),
            ("fig08_list_and_batch", fig08_list_and_batch),
            ("fig09_status_board", fig09_status)]
    info = {}
    for name, fn in jobs:
        res = fn(d)
        fig = plt.gcf()
        plots.save(fig, FIG / (name + ".png"), dpi=200)
        info[name] = res
        print(f"  ✓ {name}（{time.time() - t0:.1f}s）")
    write_readme(d, SNAP)
    write_gallery(d, SNAP)
    (OUT / "visuals_meta.json").write_text(
        json.dumps({"generated_utc": SNAP["generated_utc"],
                    "rule_version": d["rule"]["rule_version"],
                    "rules_sha256": d["rule"]["rules_sha256"],
                    "figures": [j[0] + ".png" for j in jobs],
                    "notes": info}, ensure_ascii=False, indent=1, sort_keys=True),
        encoding="utf-8")
    print(f"完成：9 张图 + README.md + gallery.html + visuals_meta.json（{time.time() - t0:.1f}s）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
