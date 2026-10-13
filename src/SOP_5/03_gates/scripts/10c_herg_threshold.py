#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3e（分析步 10c）：hERG 阈值三锚定 + 保留曲线 + 冲突消解（不改门）。

性质：**分析步** —— 主口径 hERG<0.5 不变；本步只把"阈值该怎么选"变成可复现、可审计、
可迁移到其他项目的流程证据（依据见 03_gates/admet/raw 下的 ADMET-AI 冻结表）。

三个锚定（全部在包内实算）：
  A1 库内相对  —— 包内库（f3 49,963）与包内池（4,990）的 hERG 分位；仅用于库内排序
  A2 已上市药参照 —— 用已上市药百分位列映射出 50/75/90 分位对应的 hERG 值；
                    75 分位 = 软线（「不差于 3/4 的已上市药」），50/90 分位 = 硬线合理性窗口
  A3 取舍驱动  —— 在"对接头部保留率 ≥ min_keep"约束下最严的阈值；两种头部定义分别给出：
                  ① 库级头部（f3 rank ≤ K，复现外部流程的 0.80）
                  ② 本包头部（主池 QED≥0.50 ∧ anchor=4，给出本包自洽的取舍线）

产出：
  - 03_gates/admet/tables/herg_retention_curve.csv
  - 03_gates/admet/tables/herg_threshold_anchors.csv
  - 03_gates/admet/reports/herg_threshold_analysis.json
  - 03_gates/admet/HERG_THRESHOLD_SOP.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.io_utils import (StepLog, counts_dict, r4, r6, read_json, sha256_file,
                             write_csv, write_json, write_text)

OUT = ROOT / "03_gates/admet"


def _strictest_t(values: np.ndarray, ks: list, grid: np.ndarray, min_keep: float) -> float | None:
    """在网格上取满足「保留率 ≥ min_keep」的最严（最小）阈值；网格按升序。"""
    ok = [float(t) for t in grid if float((values < t).mean()) >= min_keep]
    return min(ok) if ok else None


def _round_up(t: float | None, step: float) -> float | None:
    if t is None:
        return None
    if not step:
        return round(t, 6)
    return round(np.ceil(round(t / step, 8)) * step, 6)


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    spec = rule["step3_gates"]["herg_threshold_analysis"]
    al = rule["step3_gates"]["admet_layer"]
    with StepLog(ROOT, "10c_herg_threshold") as log:
        raw_p = ROOT / al["raw_table"]
        got = sha256_file(raw_p)
        if got != al["raw_sha256"]:
            raise RuntimeError(f"ADMET 冻结表哈希不符：{got}")
        adt = pd.read_csv(raw_p, usecols=["rank", "cid", "ADMET_hERG",
                                          "ADMET_hERG_drugbank_approved_percentile"],
                          low_memory=False)
        if adt.ADMET_hERG.isna().any():
            raise RuntimeError("库级 hERG 存在缺失值")
        ann = pd.read_csv(OUT / "tables/admet_annotation_pool.csv",
                          usecols=["cid", "ADMET_hERG", "admet_herg_band"], low_memory=False)
        univ = pd.read_csv(ROOT / "03_gates/tables/candidate_universe.csv",
                           usecols=["cid", "CNN_VS", "qed", "druglike_pass", "pose_gate",
                                    "anchor_hits_4A", "maxTc_actives", "known_exact_hit"],
                           low_memory=False)
        union = pd.read_csv(ROOT / "03_gates/herg/SOP-05_herg_input_union.csv",
                            usecols=["cid"], low_memory=False)
        ml = pd.read_csv(ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv",
                         low_memory=False)
        list48 = ml.loc[ml.candidate_type == "candidate", "cid"].astype(str)

        pool = univ.merge(ann[["cid", "ADMET_hERG"]], on="cid", how="left")
        base = pool[(pool.druglike_pass == 1) & (pool.pose_gate == "pass")
                    & (pool.maxTc_actives < 0.9) & (pool.known_exact_hit.fillna(0) == 0)]
        head_main = base[(base.qed >= 0.50) & (base.anchor_hits_4A == 4)]
        hd = {
            "library_f3": adt.ADMET_hERG.to_numpy(float),
            "pool_4990": pool.ADMET_hERG.to_numpy(float),
            "pre_herg_main_pool": head_main.ADMET_hERG.to_numpy(float),
            # 键名不再内嵌旧规模（原 quad_union_3236）；集合 = 口径重算后的四臂并集（审计件）；
            # 规模一律以数组实际长度为准。
            "quad_union_recomputed": pool[pool.cid.isin(
                set(union.cid.astype(str)))].ADMET_hERG.to_numpy(float),
            "delivered_list_48": pool[pool.cid.isin(set(list48))].ADMET_hERG.to_numpy(float),
        }
        for k in spec["alt_top_ks"]:
            sub = adt[adt["rank"] <= k]
            hd[f"library_top{k}"] = sub.ADMET_hERG.to_numpy(float)
        log.info("头部规模：" + "、".join(f"{k}={len(v)}" for k, v in hd.items()))

        grid = np.round(np.arange(spec["grid"][0], spec["grid"][1] + 1e-9, spec["grid"][2]), 2)
        step = float(spec["round_step"])
        mk = float(spec["min_keep"])

        # ---------- A1 库内相对
        a1 = {"library_f3": {str(q): r6(float(np.percentile(hd["library_f3"], q * 100)))
                             for q in (0.5, 0.75, 0.9, 0.95)},
              "pool_4990": {str(q): r6(float(np.percentile(hd["pool_4990"], q * 100)))
                            for q in (0.5, 0.75, 0.9, 0.95)},
              "pre_herg_main_pool": {str(q): r6(float(np.percentile(hd["pre_herg_main_pool"], q * 100)))
                                     for q in (0.5, 0.75, 0.9, 0.95)}}

        # ---------- A2 已上市药参照（百分位列 → hERG 值）
        pct = adt.ADMET_hERG_drugbank_approved_percentile.to_numpy(float)
        herg = adt.ADMET_hERG.to_numpy(float)
        o = np.argsort(pct, kind="mergesort")
        a2 = {str(q): r6(float(np.interp(float(q), pct[o], herg[o])))
              for q in spec["ref_percentiles"]}
        soft_line = a2["75.0"] if "75.0" in a2 else a2["75"]
        window = [min(float(v) for v in a2.values()), max(float(v) for v in a2.values())]

        # ---------- A3 取舍驱动（两种头部定义）
        heads = {"library_head@K=%d" % spec["primary_top_k"]: hd[f"library_top{spec['primary_top_k']}"],
                 "pack_main_pool": hd["pre_herg_main_pool"]}
        a3_rows, a3 = [], {}
        for name, vals in heads.items():
            t_raw = _strictest_t(vals, [], grid, mk)
            t_rnd = _round_up(t_raw, step)
            keep = r4(float((vals < t_rnd).mean()) * 100) if t_rnd is not None else None
            a3[name] = {"strictest_t_0.01_grid": r6(t_raw), "feasible_t": r6(t_rnd),
                        "retention_at_feasible_pct": keep, "head_n": int(len(vals))}
            a3_rows.append({"head": name, "min_keep": mk, "head_n": int(len(vals)),
                            "strictest_t": r6(t_raw), "feasible_t_rounded": r6(t_rnd),
                            "retention_pct": keep})
        hard_line = a3[f"library_head@K={spec['primary_top_k']}"]["feasible_t"]

        # 多 min_keep 档（供决策时查看代价）
        for keep in (0.75, 0.80, 0.90):
            for name, vals in heads.items():
                t_raw = _strictest_t(vals, [], grid, keep)
                a3_rows.append({"head": name, "min_keep": keep, "head_n": int(len(vals)),
                                "strictest_t": r6(t_raw),
                                "feasible_t_rounded": r6(_round_up(t_raw, step)),
                                "retention_pct": r4(float((vals < _round_up(t_raw, step)).mean()) * 100)
                                if t_raw is not None else None})
        write_csv(OUT / "tables/herg_threshold_anchors.csv", pd.DataFrame(a3_rows))

        # ---------- 冲突消解规则
        branch = "none（A3 落在合理性窗口内）"
        if hard_line is not None and hard_line > window[1]:
            branch = "A3 超过 A2@90 → 以 90 分位封顶"
        elif hard_line is not None and hard_line < window[0]:
            branch = "A3 低于 A2@50 → 改用 A2@75 作硬线"
        sanity_ok = bool(hard_line is not None and window[0] <= hard_line <= window[1])

        # ---------- 保留曲线
        rows = []
        for t in grid:
            r = {"t": float(t)}
            for k, v in hd.items():
                r[f"keep_{k}_pct"] = round(float((v < t).mean()) * 100, 4)
            rows.append(r)
        ret = pd.DataFrame(rows)
        write_csv(OUT / "tables/herg_retention_curve.csv", ret)

        # ---------- 名单标注档（低/中/高）
        bands = {"low": f"< {soft_line}", "mid": f"[{soft_line}, {hard_line})", "high": f">= {hard_line}"}
        band_counts = {
            "pool_4990": counts_dict(ann.admet_herg_band.fillna("covered")),
            "delivered_list_48": counts_dict(pool.loc[pool.cid.isin(set(list48)),
                                                      "ADMET_hERG"].map(
                lambda v: "高" if v >= hard_line else ("中" if v >= soft_line else "低"))),
        }

        # ---------- 门作用域缺口（非判据）：池内未落入落门表作用域的分子，其冻结表 hERG 状况
        gap = {}
        gp = ROOT / "03_gates/herg/herg_predictions.csv"
        if gp.is_file():
            gh = pd.read_csv(gp, low_memory=False)
            hcol = "hERG" if "hERG" in gh.columns else next(
                c for c in gh.columns if "herg" in c.lower())
            thr = float(rule["step3_gates"]["herg"]["threshold"])
            oos = pool.loc[~pool.cid.astype(str).isin(set(gh.cid.astype(str)))]
            j = gh.assign(cid=gh.cid.astype(str)).merge(
                pool.assign(cid=pool.cid.astype(str))[["cid", "ADMET_hERG"]], on="cid", how="inner")
            gap = {
                "gate_table": "03_gates/herg/herg_predictions.csv",
                "gate_table_n": int(len(gh)),
                "gate_table_sha256": sha256_file(gp),
                "pool_n": int(len(pool)),
                "in_scope": int(len(pool) - len(oos)),
                "out_of_scope": int(len(oos)),
                "out_of_scope_frozen_below_gate": int((oos.ADMET_hERG < thr).sum()),
                "out_of_scope_frozen_at_or_above_gate": int((oos.ADMET_hERG >= thr).sum()),
                "out_of_scope_below_gate_pct": r4(float((oos.ADMET_hERG < thr).mean()) * 100),
                "overlap_n": int(len(j)),
                "overlap_max_abs_diff": r6(float((j[hcol] - j.ADMET_hERG).abs().max())),
                "note": "落门表只覆盖四臂并集（外部运行的作用域）；并集外池内分子无预测值 → 按规则"
                        "保守判为不过门。本包冻结的 f3 全量表现已给出这些分子的 ADMET-AI hERG，"
                        "故可量化该保守口径的代价；**门与名单不因此改动**（改作用域属口径变更，"
                        "需重跑选名单与验收基线）。落门表与冻结表在同源交集内逐值一致（见 overlap 项）。",
            }
            log.info(f"门作用域缺口：池内并集外 {gap['out_of_scope']} 个，其中冻结表 hERG < "
                     f"{thr} 的有 {gap['out_of_scope_frozen_below_gate']} 个"
                     f"（{gap['out_of_scope_below_gate_pct']}%）；交集 {gap['overlap_n']} 个逐值 "
                     f"max|Δ|={gap['overlap_max_abs_diff']}")

        # ---------- 复现核对（与外部流程的数字对账）
        rc = spec["repro_check"]
        repro = {
            "soft_line_expected": rc["soft_line_expected"], "soft_line_observed": r6(soft_line),
            "soft_line_ok": bool(abs(float(soft_line) - float(rc["soft_line_expected"])) <= rc["tol"]),
            "library_head_tradeoff_expected": rc["library_head_tradeoff_expected"],
            "library_head_tradeoff_observed": r6(hard_line),
            "library_head_tradeoff_ok": bool(
                hard_line is not None
                and abs(float(hard_line) - float(rc["library_head_tradeoff_expected"])) <= rc["tol"]),
            "library_head_at_0.01_grid": a3[f"library_head@K={spec['primary_top_k']}"]["strictest_t_0.01_grid"],
            "note": rc["note"],
        }
        log.info(f"A2 软线 {r6(soft_line)}（期望 {rc['soft_line_expected']}，"
                 f"ok={repro['soft_line_ok']}）；A3 库级取舍线 {r6(hard_line)}"
                 f"（期望 {rc['library_head_tradeoff_expected']}，ok={repro['library_head_tradeoff_ok']}）；"
                 f"本包头部取舍线 {a3['pack_main_pool']['feasible_t']}")

        rep = {
            "role": spec["role"],
            "provenance": {"raw_table": al["raw_table"], "raw_sha256": got,
                           "library_rows": int(len(adt)), "producer": al["producer"]},
            "params": {"grid": spec["grid"], "round_step": step, "min_keep": mk,
                       "primary_top_k": spec["primary_top_k"], "alt_top_ks": spec["alt_top_ks"],
                       "ref_percentiles": spec["ref_percentiles"]},
            "A1_library_relative": a1,
            "A2_approved_drug_reference": a2,
            "A3_tradeoff": a3,
            "decision": {"frozen_gate_threshold": rule["step3_gates"]["herg"]["threshold"],
                         "soft_line": r6(soft_line), "hard_line_candidate": r6(hard_line),
                         "hard_line_candidate_pack_head": a3["pack_main_pool"]["feasible_t"],
                         "sanity_window": [r6(window[0]), r6(window[1])], "sanity_ok": sanity_ok,
                         "resolution_branch": branch,
                         "bands": bands, "band_counts": band_counts,
                         "note": spec["role"]},
            "retention_curve_file": "03_gates/admet/tables/herg_retention_curve.csv",
            "anchors_file": "03_gates/admet/tables/herg_threshold_anchors.csv",
            "repro_check": repro,
            "gate_scope_gap": gap,
            "non_interference": {
                "statement": "本步为分析步：不改变 hERG 门（主口径仍 <0.5）、不改池、不改名单；"
                             "硬线候选仅供将来改门时决策引用",
                "gate_threshold": rule["step3_gates"]["herg"]["threshold"]},
        }
        write_json(OUT / "reports/herg_threshold_analysis.json", rep)

        md = ["# hERG 阈值选择 SOP（包内版本 · 步骤 3e / 分析步 10c）", "",
              "**目的**：把「用 hERG 过滤分子」的阈值选择变成一次可复现、可审计、可迁移的流程。",
              "**本版本依据**：包内冻结的 ADMET-AI f3 库表（49,963 分子 × 105 个 ADMET 列，"
              f"sha256 `{got[:16]}…`）与其中的已上市药参照百分位列。",
              "**性质**：分析步——主口径 `hERG < "
              f"{rule['step3_gates']['herg']['threshold']}` 不变，本步只产出决策证据。", "",
              "## 1 三个锚定各自回答什么（分工，不是三选一）", "",
              "| 锚定 | 回答的问题 | 角色 | 本包取值 |", "|---|---|---|---|",
              f"| ① 库内相对 | 在这批分子里它排在哪 | 排序层（不做硬线） | 库 75/90 分位 = "
              f"{a1['library_f3']['0.75']} / {a1['library_f3']['0.9']}；池 75/90 分位 = "
              f"{a1['pool_4990']['0.75']} / {a1['pool_4990']['0.9']} |",
              f"| ② 已上市药参照 | 相对已上市药算不算差 | 解释层（软线 + 合理性窗口） | "
              f"50/75/90 分位 = {a2[str(spec['ref_percentiles'][0])]} / "
              f"{a2[str(spec['ref_percentiles'][1])]} / {a2[str(spec['ref_percentiles'][2])]} |",
              f"| ③ 取舍驱动 | 过滤会牺牲多少对接头部 | 决策层（硬线候选） | 库级头部 "
              f"{r6(hard_line)}；本包主池头部 {a3['pack_main_pool']['feasible_t']} |", "",
              "## 2 执行步骤（确定性）", "",
              f"1. **定硬线（③）**：在 {spec['grid'][2]} 网格上取"
              f"「保住对接前 K={spec['primary_top_k']} 名 ≥ {mk:.0%}」的最严阈值 → "
              f"{r6(hard_line)}（0.01 网格精确解 "
              f"{a3['library_head@K=%d' % spec['primary_top_k']]['strictest_t_0.01_grid']}），"
              f"再按 {step} 网格取整以免贴边。",
              f"2. **合理性窗口（②）**：硬线须落在 [A2@50, A2@90] = "
              f"[{r6(window[0])}, {r6(window[1])}] → "
              f"{'通过' if sanity_ok else '需按 §3 冲突消解处理'}。",
              f"3. **定软线（②）**：软线 = 参照集 75 分位 = {r6(soft_line)}"
              f"（表述：「不差于 3/4 的已上市药」）。",
              "4. **输出排序列（①）**：附库内百分位，供系列内挑类似物使用。", "",
              "## 3 冲突消解规则（其他库/其他项目出现分歧时）", "",
              "- 若 ③ 解 > ② 的 90 分位 → **以 90 分位封顶**（不允许为了保对接头部突破"
              "「不劣于 90% 已上市药」底线）。",
              "- 若 ③ 解 < ② 的 50 分位 → 说明对接头部整体偏「比一半已上市药还差」，",
              "改用参照集 75 分位作硬线，并回头检查打分与化学空间。",
              "- 若 ① 与 ② 相差 > 20 个百分位 → 库整体偏移严重，对外一律用 ② 表述。",
              "- 任何情况下 ① 不单独决定过滤。", "",
              f"本轮分支判定：**{branch}**。", "",
              "## 4 本包当前判定（不改门，仅标注）", "",
              "| 档 | 规则 | 池内分子数 | 交付 48 候选 |", "|---|---|---|---|",
              f"| 高 | hERG ≥ {hard_line} | {band_counts['pool_4990'].get('高', 0)} | "
              f"{band_counts['delivered_list_48'].get('高', 0)} |",
              f"| 中 | [{soft_line}, {hard_line}) | {band_counts['pool_4990'].get('中', 0)} | "
              f"{band_counts['delivered_list_48'].get('中', 0)} |",
              f"| 低 | < {soft_line} | {band_counts['pool_4990'].get('低', 0)} | "
              f"{band_counts['delivered_list_48'].get('低', 0)} |", "",
              "主口径的门仍是 `hERG < "
              f"{rule['step3_gates']['herg']['threshold']}`；上表「中/高」档在交付名单里只作"
              "MD 前置复核的优先级标注，**不删分子**。", "",
              "## 4b 门作用域缺口（非判据，量化保守口径的代价）", "",
              (f"落门表 `{gap['gate_table']}` 只覆盖四臂并集 {gap['gate_table_n']} 个分子；池内"
               f"{gap['out_of_scope']} 个分子在并集之外（无预测值）→ 按规则保守判为不过门。"
               f"本包冻结的 f3 全量表给出了这些分子的 ADMET-AI hERG：其中 "
               f"**{gap['out_of_scope_frozen_below_gate']} 个（"
               f"{gap['out_of_scope_below_gate_pct']}%）实际 < 门阈值 {thr}**，即被保守口径排除。"
               f"落门表与冻结表在同源交集（{gap['overlap_n']} 个）内逐值一致"
               f"（max|Δ| = {gap['overlap_max_abs_diff']}）。"
               if gap else
               "（落门表缺失，此项不可得——门为 pending，不影响门与名单。）")
              if gap else "（落门表缺失，此项不可得。）", "",
              "这些分子**不在任何臂/名单内**，本项只作透明度披露：若要改作用域，属口径变更"
              "（需重跑选名单与验收基线），本步不改门、不改名单。", "",
              "## 5 复现与产物", "",
              "```bash",
              "python3 run_all.py --from 10c_herg_threshold",
              "```", "",
              "- `03_gates/admet/tables/herg_retention_curve.csv`（0.01 步长保留曲线，覆盖库/库级头部/"
              "池/主池/四臂并集/交付名单）",
              "- `03_gates/admet/tables/herg_threshold_anchors.csv`（各头部 × min_keep 档的取舍线）",
              "- `03_gates/admet/reports/herg_threshold_analysis.json`（含输入 sha256 与复现核对）", "",
              "## 6 变更控制", "",
              "- 必须重跑本流程的情形：合并/更换化合物库；ADMET 冻结表重算（模型或版本变化）；"
              "对接排名更新；`grid`/`min_keep`/`primary_top_k`/`round_step` 改动。",
              "- 阈值只在 "
              f"{step} 网格取值（避免贴边解）；每次变更在 decision JSON 与本文件中留痕。",
              "- 冻结项：软线、硬线候选、判定字典（除非流程重跑，不得手改）。", "",
              "## 7 局限", "",
              "- hERG 预测未做概率校准；硬线候选是「本包库 + 本对接目标」下的操作点，"
              "**迁移到其他项目必须重跑本流程**。",
              "- 库内相对锚定（①）在本包以 f3 49,963 行为口径；外部流程报告用的是"
              "54,932 去重并集，数字会有小差异（不影响②③）。",
              "- 本流程只做优先级与标注；最终结论需离子通道实验（膜片钳等）确认。"]
        write_text(OUT / "HERG_THRESHOLD_SOP.md", "\n".join(md) + "\n")
        log.info("阈值锚定完成：HERG_THRESHOLD_SOP.md + 2 个 CSV + 1 个 JSON（门未改动）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
