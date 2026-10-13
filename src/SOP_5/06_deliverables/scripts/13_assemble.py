#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""交付层（上）：汇总 6 组名单矩阵 + 对照与决策卡 + 规范化 SDF 子集 + 验收对账 + MANIFEST。

产出（06_deliverables/）：
  SOP-05_MD_shortlist_main.csv              主名单（50：48 候选 + 2 阴性对照）
  SOP-05_MD_shortlist_sensitivity_*.csv     5 个敏感性名单
  SOP-05_MD_shortlist_scaffoldcap{3,5}.csv  骨架封顶敏感性档
  SOP-05_decision_cards_main.csv            主名单决策卡（逐分子证据）
  SOP-05_negative_controls.csv              阴性对照卡（含重算规则与实测值）
  SOP-05_poses_main_list.sdf                主名单起始姿态子集
  acceptance_report.json                    验收对账（逐条判据 + 实测）
  MANIFEST.sha256                           交付文件哈希清单
  INDEX.md                                  交付索引
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor import chem_core as cc
from vendor.io_utils import (StepLog, manifest_lines, read_json, r4, write_csv, write_json,
                             write_text)

DELIV = ROOT / "06_deliverables"


def _rob_note(row: dict, rows: list, threshold: float) -> str:
    """A 臂主判据的实测转述：与门槛比较后如实描述通过/未通过，不写死结论与数字。"""
    obs = float(row["mean_overlap"])
    verdict = ("未通过：实测低于门槛 " if obs < float(threshold) else "通过：实测 ≥ 门槛 ")
    others = "、".join(f"σ={r['sigma']} 为 {r['mean_overlap']}"
                      for r in rows if r["sigma"] != row["sigma"])
    return (f"{verdict}{threshold}（实测 {obs}；min {row['min_overlap']}、P10 "
            f"{row['p10_overlap']}）；{others}。门槛与判据均未调整，缺口如实保留")

COPY_MAP = {
    "SOP-05_MD_shortlist_main.csv": "04_selection/tables/list_main_QED0.50_anchor4.csv",
    "SOP-05_MD_shortlist_sens_QED0.60_anchor4.csv": "04_selection/tables/list_sens_QED0.60_anchor4.csv",
    "SOP-05_MD_shortlist_sens_QED0.40_anchor4.csv": "04_selection/tables/list_sens_QED0.40_anchor4.csv",
    "SOP-05_MD_shortlist_sens_QED0.50_anchor3.csv": "04_selection/tables/list_sens_QED0.50_anchor3.csv",
    "SOP-05_MD_shortlist_sens_QED0.60_anchor3.csv": "04_selection/tables/list_sens_QED0.60_anchor3.csv",
    "SOP-05_MD_shortlist_sens_QED0.40_anchor3.csv": "04_selection/tables/list_sens_QED0.40_anchor3.csv",
    "SOP-05_MD_shortlist_scaffoldcap3.csv": "04_selection/tables/list_sens_scaffoldcap3_QED0.50_anchor4.csv",
    "SOP-05_MD_shortlist_scaffoldcap5.csv": "04_selection/tables/list_sens_scaffoldcap5_QED0.50_anchor4.csv",
    "SOP-05_decision_cards_main.csv": "04_selection/tables/decision_cards_main_QED0.50_anchor4.csv",
    "SOP-05_constraint_trigger_rates.csv": "04_selection/tables/constraint_trigger_rates.csv",
    "SOP-05_gate_stack.csv": "03_gates/tables/gate_stack.csv",
    "SOP-05_robustness.json": "05_robustness/reports/robustness.json",
    "SOP-05_determinism_compare.json": "05_robustness/reports/determinism_compare.json",
    "SOP-05_herg_input_union.csv": "03_gates/herg/SOP-05_herg_input_union.csv",
    "SOP-05_herg_status.json": "03_gates/reports/herg_status.json",
    "SOP-05_herg_cross_model_consistency.csv": "03_gates/herg/cross_model/tables/herg_cross_model_agreement.csv",
    "SOP-05_herg_cross_model_rank_consistency.csv": "03_gates/herg/cross_model/tables/herg_cross_model_rank_consistency.csv",
    "SOP-05_herg_cross_model_list_flips.csv": "03_gates/herg/cross_model/tables/herg_cross_model_list_flips.csv",
    "SOP-05_admet_md_precheck.csv": "03_gates/admet/tables/admet_md_precheck.csv",
    "SOP-05_admet_descriptor_crosscheck.csv": "03_gates/admet/tables/admet_descriptor_crosscheck.csv",
    "SOP-05_herg_threshold_anchors.csv": "03_gates/admet/tables/herg_threshold_anchors.csv",
    "SOP-05_herg_retention_curve.csv": "03_gates/admet/tables/herg_retention_curve.csv",
    "SOP-05_MD_batch20.csv": "07_md_batch/tables/md_batch20.csv",
    "SOP-05_MD_batch_excluded.csv": "07_md_batch/tables/md_batch_excluded.csv",
    "SOP-05_MD_batch20.md": "07_md_batch/MD_BATCH.md",
    "SOP-05_md_batch.json": "07_md_batch/reports/md_batch.json",
    "SOP-05_admet_annotation.json": "03_gates/admet/reports/admet_annotation.json",
    "SOP-05_admet_descriptor_crosscheck.json": "03_gates/admet/reports/admet_descriptor_crosscheck.json",
    "SOP-05_herg_threshold_analysis.json": "03_gates/admet/reports/herg_threshold_analysis.json",
    "SOP-05_admet_layer.md": "03_gates/admet/ADMET_LAYER.md",
    "SOP-05_HERG_THRESHOLD_SOP.md": "03_gates/admet/HERG_THRESHOLD_SOP.md",
}

# 需要追加 admet_* 注释列的交付物（列追加；行序与成员不变）
ANNOTATE = [
    "SOP-05_MD_shortlist_main.csv",
    "SOP-05_MD_shortlist_sens_QED0.60_anchor4.csv",
    "SOP-05_MD_shortlist_sens_QED0.40_anchor4.csv",
    "SOP-05_MD_shortlist_sens_QED0.50_anchor3.csv",
    "SOP-05_MD_shortlist_sens_QED0.60_anchor3.csv",
    "SOP-05_MD_shortlist_sens_QED0.40_anchor3.csv",
    "SOP-05_MD_shortlist_scaffoldcap3.csv",
    "SOP-05_MD_shortlist_scaffoldcap5.csv",
    "SOP-05_decision_cards_main.csv",
    "SOP-05_negative_controls.csv",
]
FIG_MAP = {
    "SOP-05_fig_pose_gate.png": "01_pose_qc/figures/fig_pose_gate.png",
    "SOP-05_fig_druglikeness_funnel.png": "02_druglikeness/figures/fig_druglikeness_funnel.png",
    "SOP-05_fig_gate_stack.png": "03_gates/figures/fig_gate_stack.png",
    "SOP-05_fig_list_overlap.png": "04_selection/figures/fig_list_overlap.png",
    "SOP-05_fig_gate_impact.png": "04_selection/figures/fig_gate_impact.png",
    "SOP-05_fig_robA_sigma.png": "05_robustness/figures/fig_robA_sigma.png",
    "SOP-05_fig_robB_matrix.png": "05_robustness/figures/fig_robB_matrix.png",
    "SOP-05_fig_robE_band.png": "05_robustness/figures/fig_robE_band.png",
    "SOP-05_fig_robF_clusters.png": "05_robustness/figures/fig_robF_clusters.png",
    "SOP-05_fig_md_batch.png": "07_md_batch/figures/fig_md_batch.png",
    "SOP-05_fig_herg_cross_model.png": "03_gates/herg/cross_model/figures/fig_herg_cross_model.png",
}


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    with StepLog(ROOT, "13_assemble") as log:
        DELIV.mkdir(parents=True, exist_ok=True)
        (DELIV / "figures").mkdir(parents=True, exist_ok=True)
        for dst, src in COPY_MAP.items():
            s = ROOT / src
            if s.is_file():
                shutil.copy2(s, DELIV / dst)
            else:
                log.info(f"[缺失] {src}（跳过 {dst}）")
        for dst, src in FIG_MAP.items():
            s = ROOT / src
            if s.is_file():
                shutil.copy2(s, DELIV / "figures" / dst)
                svg = s.with_suffix(".svg")
                if svg.is_file():
                    shutil.copy2(svg, DELIV / "figures" / dst.replace(".png", ".svg"))
            else:
                log.info(f"[缺失] {src}（跳过 {dst}）")

        # ---- 阴性对照卡
        neg = pd.read_csv(ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv",
                          low_memory=False)
        neg = neg[neg.candidate_type == "negative_control"]
        write_csv(DELIV / "SOP-05_negative_controls.csv", neg)
        exp = rule["step4_selection"]["controls"]
        neg_match = {
            "NEG1": bool(exp["negative_1"]["expected_cid"] in set(neg.cid.astype(str))),
            "NEG2": bool(exp["negative_2"]["expected_cid"] in set(neg.cid.astype(str))),
        }
        log.info(f"阴性对照规则复现：{neg_match}（NEG1/NEG2 = "
                 f"{exp['negative_1']['expected_cid']} / {exp['negative_2']['expected_cid']}）")

        # ---- 主名单起始姿态子集
        want = set(pd.read_csv(DELIV / "SOP-05_MD_shortlist_main.csv", low_memory=False)
                   .loc[lambda d: d.candidate_type == "candidate", "cid"].astype(str))
        out = open(DELIV / "SOP-05_poses_main_list.sdf", "w", encoding="utf-8")
        written = 0
        for block, props, title in cc.iter_sdf_blocks(
                ROOT / "01_pose_qc/outputs/SOP-05_poses_normalized.sdf"):
            if props.get("cid", title) in want:
                core = cc.molblock_data_lines(block)      # 保证以 M  END 收尾
                out.write("\n".join(core) + "\n$$$$\n")
                written += 1
        out.close()
        log.info(f"主名单候选姿态子集：{written}/{len(want)} 条 → "
                 f"06_deliverables/SOP-05_poses_main_list.sdf")

        # ---- ADMET 注释列追加（非判据；列追加、行序与成员不变）
        ann_p = ROOT / "03_gates/admet/tables/admet_annotation_cards.csv"
        ann_n = 0
        if ann_p.is_file() and rule.get("step3_gates", {}).get("admet_layer"):
            ann_cols = rule["step3_gates"]["admet_layer"]["annotation_columns"]
            cards = pd.read_csv(ann_p, low_memory=False).set_index("cid")
            ann_member_ok = True
            for name in ANNOTATE:
                path = DELIV / name
                if not path.is_file():
                    continue
                df = pd.read_csv(path, dtype=str, keep_default_na=False)
                before = df.copy()          # 追加前快照：用于断言「只加列、不动行与值」
                for c in ann_cols:
                    vals = []
                    for cid in df.cid.astype(str):
                        if cid in cards.index:
                            v = cards.at[cid, c]
                            if c == "admet_coverage":
                                vals.append(str(v))
                            elif isinstance(v, float):
                                vals.append("" if pd.isna(v) else f"{v:.6g}")
                            else:
                                vals.append(str(v) if str(v) not in ("", "nan") else "")
                        else:
                            vals.append("not_covered" if c == "admet_coverage" else "")
                    df[c] = vals
                same = (df.cid.tolist() == before.cid.tolist()
                        and all(df[c].tolist() == before[c].tolist() for c in before.columns))
                if not same:
                    ann_member_ok = False
                    raise RuntimeError(f"注释列追加改变了 {name} 的成员/顺序/原列取值")
                df.to_csv(path, index=False)
                ann_n += 1
            # 交付名单注释全表（逐分子全列，冻结表原精度）
            pool_p = ROOT / "03_gates/admet/tables/admet_annotation_pool.csv"
            main_cids = pd.read_csv(DELIV / "SOP-05_MD_shortlist_main.csv",
                                    low_memory=False).cid.astype(str).tolist()
            if pool_p.is_file():
                pool = pd.read_csv(pool_p, low_memory=False).set_index("cid")
                wide = pool.reindex([c for c in main_cids if c in pool.index]).reset_index()
                miss = [c for c in main_cids if c not in pool.index]
                for cid in miss:                     # 阴性对照：显式标注未覆盖
                    row = {k: "" for k in wide.columns}
                    row["cid"] = cid
                    row["admet_coverage"] = "not_covered"
                    wide = pd.concat([wide, pd.DataFrame([row])], ignore_index=True)
                order = {c: i for i, c in enumerate(main_cids)}
                wide = wide.sort_values(by="cid", key=lambda s: s.map(order)).reset_index(drop=True)
                wide.to_csv(DELIV / "SOP-05_admet_annotation_main.csv", index=False)
            log.info(f"ADMET 注释列追加：{ann_n} 个交付物（仅加列；cid 顺序与原列逐值不变已断言"
                     f"={ann_member_ok}；门/池/名单未改动）；"
                     f"注释全表 → SOP-05_admet_annotation_main.csv")
        else:
            ann_member_ok = None
            log.info("ADMET 注释表缺失（未运行 10b）→ 跳过注释列追加（不影响门与名单）")

        # ---- 验收对账
        summ = read_json(ROOT / "04_selection/reports/selection_summary.json")
        main = summ["lists"]["main_QED0.50_anchor4"]
        ace = rule["acceptance"]
        herg = read_json(ROOT / "03_gates/reports/herg_status.json")
        _detp = ROOT / "05_robustness/reports/determinism_compare.json"
        det = (read_json(_detp) if _detp.is_file() else
               {"mode": "pending", "identical": None,
                "note": "本轮首跑：G 臂（12_determinism）在整链末尾执行，验收表此行标 pending；"
                        "第二次整链运行起回填实测值"})
        rob = read_json(ROOT / "05_robustness/reports/robustness.json")
        armsig = [r for r in rob["armA_score_perturbation"]["rows"] if r["sigma"] == 0.10][0]
        fun_dl = read_json(ROOT / "02_druglikeness/reports/druglikeness_funnel.json")
        gs = read_json(ROOT / "03_gates/reports/gate_stack.json")
        checks = [
            {"criterion": "主名单锚点 100% = 4/4", "class": "hard",
             "expected": ace["anchor_all_4_in_main_list"],
             "observed": main["anchor_all_4"], "pass": bool(main["anchor_all_4"])},
            {"criterion": f"成对 Tc 上限 ≤ {ace['pairwise_tc_max_le']}", "class": "hard",
             "expected": ace["pairwise_tc_max_le"], "observed": main["max_pairwise_tc"],
             "pass": bool(main["max_pairwise_tc"] <= ace["pairwise_tc_max_le"])},
            {"criterion": "唯一骨架数（报告值）", "class": "reported", "expected": "报告",
             "observed": main["n_distinct_scaffold"], "pass": True},
            {"criterion": "hERG 区间", "class": "reported", "expected": "<0.5（外部表可得时）",
             "observed": main["herg_range"], "pass": bool(herg["status"] in ("applied", "pending")),
             "note": (f"门状态 {herg['status']}、作用域 {herg.get('scope_mode')}，"
                      f"覆盖 {float(herg.get('scope_coverage', 0)) * 100:.2f}%"
                      + ("；未评估分子记 not_evaluated，不当作不合格"
                         if herg.get("n_not_evaluated") else ""))},
            {"criterion": "已知物精确重合数", "class": "hard",
             "expected": ace["known_compound_exact_hits"], "observed": main["known_exact_hits"],
             "pass": bool(main["known_exact_hits"] == ace["known_compound_exact_hits"])},
            {"criterion": f"σ=0.10 名单平均重叠 ≥ {ace['mean_list_overlap_min_under_sigma0.10']}（主判据）",
             "class": "hard", "expected": ace["mean_list_overlap_min_under_sigma0.10"],
             "observed": armsig["mean_overlap"],
             "pass": bool(armsig["mean_overlap"] >= ace["mean_list_overlap_min_under_sigma0.10"]),
             "note": (_rob_note(armsig, rob["armA_score_perturbation"]["rows"],
                                ace["mean_list_overlap_min_under_sigma0.10"]))},
            {"criterion": "σ=0.10 头部 20 重叠（次要指标，如实报告、不设门）",
             "class": "reported", "expected": "报告（原 0.8 的表述歧义见规则 1.1.0 澄清）",
             "observed": armsig["mean_head20_overlap"], "pass": None,
             "note": f"实测 {armsig['mean_head20_overlap']}，未达 0.80：头部 20 名对分数扰动更敏感；"
                     f"缺口如实保留"},
            {"criterion": "每一步有漏斗 JSON", "class": "hard",
             "expected": ["02_druglikeness_funnel.json", "03_gates_gate_stack.json",
                          "04_selection_summary.json", "05_robustness.json"],
             "observed": ["有" if (ROOT / q).is_file() else "缺" for q in
                          ["02_druglikeness/reports/druglikeness_funnel.json",
                           "03_gates/reports/gate_stack.json",
                           "04_selection/reports/selection_summary.json",
                           "05_robustness/reports/robustness.json"]],
             "pass": all((ROOT / q).is_file() for q in
                         ["02_druglikeness/reports/druglikeness_funnel.json",
                          "03_gates/reports/gate_stack.json",
                          "04_selection/reports/selection_summary.json",
                          "05_robustness/reports/robustness.json"])},
            {"criterion": "确定性重跑 sha256 全量一致（G 臂）", "class": "hard",
             "expected": True, "observed": None, "pass": None,
             "status": "pending（由链末 G 臂判定；本步不自行判定）",
             "note": (f"本步（13_assemble）在 G 臂之前执行，读到的 "
                      f"mode={det.get('mode')}、n_diff={det.get('n_diff')} 是**上一轮链末**结果，"
                      f"不得当作本轮确定性通过证据；本轮判定写入 "
                      f"05_robustness/reports/determinism_compare.json（链末）与 "
                      f"tests/verify_pack.py 第 [17] 组")},
            {"criterion": "确定性基线与比对运行分离（mode=compared 才算判据）",
             "class": "hard",
             "expected": "mode=compared 且 n_diff=0（链末写入，见 determinism_verdict）",
             "observed": None, "pass": None,
             "status": "pending（由链末 G 臂判定；本步不自行判定）",
             "note": ("判据要求整链不可自证：写基线的一轮（mode=baseline_written）"
                      "determinism_verified=false。本表由交付层生成，位于 G 臂之前，"
                      "故两行在**本轮**恒为 pending；最终结果见 05_robustness/reports/"
                      "determinism_compare.json")},
            {"criterion": "阴性对照规则复现（包内重算）", "class": "hard",
             "expected": [True, True], "observed": neg_match,
             "pass": all(neg_match.values())},
            {"criterion": "阴性对照 2 个（固定占位）", "class": "hard", "expected": 2,
             "observed": int(len(neg)), "pass": bool(len(neg) == 2)},
        ]
        hard = [c for c in checks if c.get("class") == "hard"]
        hard_meas = [c for c in hard if c["pass"] is not None]
        n_pending = len(hard) - len(hard_meas)
        acc = {"rule_set_id": rule["rule_set_id"],
               "hard_gate_pass": all(c["pass"] for c in hard_meas),
               "n_hard": len(hard), "n_hard_pass": sum(1 for c in hard_meas if c["pass"]),
               "n_hard_pending": n_pending,
               "all_pass": (all(c["pass"] for c in hard_meas) and n_pending == 0),
               "checks": checks,
               "determinism_verdict": {
                   "judged_by": "05_robustness/scripts/12_determinism.py（整链末步）",
                   "source": "05_robustness/reports/determinism_compare.json",
                   "last_recorded": {"mode": det.get("mode"), "n_diff": det.get("n_diff"),
                                     "determinism_verified": det.get("determinism_verified"),
                                     "baseline_sha256": det.get("baseline_sha256")},
                   "note": ("本步在 G 臂之前执行，此处只登记上一轮链末读数；"
                            "本轮判据以链末写入的同名文件为准，并由 "
                            "tests/verify_pack.py 第 [17] 组断言（mode=compared 且 n_diff=0）"),
               },
               "context": {
                   "funnel_druglikeness": fun_dl["funnel"],
                   "gate_chain": gs["gate_chain"],
                   "arm_pool_sizes": gs["arm_pools"],
                   "hERG": gs["hERG"],
                   "known_ligand_check": gs["known_ligand_check"],
                   "main_list": main,
                   "constraint_trigger": pd.read_csv(
                       ROOT / "04_selection/tables/constraint_trigger_rates.csv").to_dict("records"),
               }}
        _ad_p = ROOT / "03_gates/admet/reports/admet_annotation.json"
        _ht_p = ROOT / "03_gates/admet/reports/herg_threshold_analysis.json"
        if _ad_p.is_file() and _ht_p.is_file():
            _ad_j, _ht_j = read_json(_ad_p), read_json(_ht_p)
            _xc_j = read_json(ROOT / "03_gates/admet/reports/admet_descriptor_crosscheck.json") \
                if (ROOT / "03_gates/admet/reports/admet_descriptor_crosscheck.json").is_file() else {}
            _main_p = ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv"
            _deliv_p = DELIV / "SOP-05_MD_shortlist_main.csv"
            _membership_ok = bool(
                _deliv_p.is_file()
                and pd.read_csv(_deliv_p, dtype=str, keep_default_na=False).cid.tolist()
                == pd.read_csv(_main_p, dtype=str, keep_default_na=False).cid.tolist())
            acc["context"]["admet_annotation_layer"] = {
                "non_gating": True,
                "source": _ad_j["provenance"]["raw_table"],
                "source_sha256": _ad_j["provenance"]["raw_sha256"],
                "pool_coverage": _ad_j["coverage"]["covered"],
                "not_covered": _ad_j["coverage"]["not_covered"],
                "cutpoints_recalibrated_on_pool": _ad_j["cutpoints_recalibrated_on_pool"],
                "delivered_list": _ad_j["delivered_list"],
                "descriptor_crosscheck": "03_gates/admet/reports/admet_descriptor_crosscheck.json",
                "descriptor_crosscheck_tolerance": _xc_j.get("tolerance"),
                "descriptor_crosscheck_all_within_tol": _xc_j.get("all_numeric_within_tol"),
                "descriptor_crosscheck_rule_agreement_pct": {
                    k: v.get("agreement_pct") for k, v in (_xc_j.get("rule_agreements") or {}).items()},
                "prediction_layer_recomputed_in_pack": False,
                "list_membership_unchanged": _membership_ok,
                "annotation_merge_membership_asserted": ann_member_ok,
                "herg_soft_line": _ht_j["decision"]["soft_line"],
                "herg_hard_line_candidate": _ht_j["decision"]["hard_line_candidate"],
                "herg_threshold_repro_check": _ht_j["repro_check"],
                "herg_gate_unchanged": _ht_j["decision"]["frozen_gate_threshold"],
                "note": "非判据：不改门（hERG<0.5）/不改池/不改名单；注释列为交付物追加列",
            }
        write_json(DELIV / "acceptance_report.json", acc)
        log.info(f"验收：硬门 {acc['n_hard_pass']}/{len(hard)} 通过"
                 f"（pending {n_pending}）；次要指标 "
                 f"{sum(1 for c in checks if c.get('class') == 'reported')} 项如实报告；"
                 f"all_pass={acc['all_pass']}")

        # ---- MANIFEST + INDEX
        files = [f for f in DELIV.rglob("*") if f.is_file() and f.name != "MANIFEST.sha256"]
        write_text(DELIV / "MANIFEST.sha256",
                   "\n".join(manifest_lines(ROOT, files, "sha256")) + "\n\n"
                   + "\n".join(f"# md5  {l.split('  ', 1)[0][:32]}  {l.split('  ', 1)[1]}"
                               for l in manifest_lines(ROOT, files, "md5")))
        idx = ["# SOP_5 交付索引（06_deliverables）", "",
               f"- 规则集：`{rule['rule_set_id']}` v{rule['rule_version']}"
               f"（rules_sha256 `{rule.get('rules_sha256', '')[:16]}…`）",
               f"- 验收：**硬门 {acc['n_hard_pass']}/{acc['n_hard']} 通过"
               f"{'' if acc['n_hard_pending'] == 0 else f'（pending {acc['n_hard_pending']}）'}**；"
               f"次要指标 "
               f"{sum(1 for c in checks if c.get('class') == 'reported')} 项如实报告"
               f"（详见 `acceptance_report.json`）", "",
               "| 文件 | 说明 |", "|---|---|"]
        notes = {
            "SOP-05_MD_shortlist_main.csv": "主名单 50（48 候选 + 2 阴性对照），含逐分子门与注释",
            "SOP-05_decision_cards_main.csv": "主名单决策卡（含最近的已知活性、pKa/cLogD、PAINS/Brenk 标注）",
            "SOP-05_negative_controls.csv": "阴性对照卡（包内重算 NEG1/NEG2）",
            "SOP-05_poses_main_list.sdf": ("主名单候选起始姿态子集（数据行与上传姿态逐字节一致，"
                                           "块以 M  END 收尾）"),
            "acceptance_report.json": "验收对账：逐条判据 + 实测 + 上下文（硬门 / 次要指标分类）",
            "SOP-05_report.md": "最终报告（14_report.py 生成）",
            "SOP-05_downstream_MD_plan.md": "下游 MD 接口与判据（含机时预算单列）",
            "SOP-05_gate_stack.csv": "逐分子门叠加结果（4,990 条）",
            "SOP-05_herg_input_union.csv": "hERG 补算清单（四臂并集，外部接口）",
            "SOP-05_herg_status.json": "hERG 门状态（缺表时 pending）",
            "SOP-05_robustness.json": "稳健性 A–G 全量结果",
            "SOP-05_determinism_compare.json": "确定性比对报告（G 臂；最近一次结果）",
            "SOP-05_constraint_trigger_rates.csv": "选序约束触发率（Tc/骨架/化学型/风险档）",
            "SOP-05_admet_annotation_main.csv": "交付名单 ADMET 注释全表（冻结表原精度，50 行）",
            "SOP-05_admet_md_precheck.csv": "MD 前置复核清单（交付候选中带 ≥1 项警告者）",
            "SOP-05_admet_descriptor_crosscheck.csv": "描述符独立复算交叉核验（逐分子偏差与判定一致性）",
            "SOP-05_admet_descriptor_crosscheck.json": "交叉核验汇总（容差与已知差异说明）",
            "SOP-05_admet_annotation.json": "注释层证据（覆盖/重标定切点/触发率/分诊层/复核清单）",
            "SOP-05_admet_layer.md": "注释层说明（步骤 10b）",
            "SOP-05_herg_threshold_anchors.csv": "hERG 阈值三锚定表（各头部 × min_keep 档取舍线）",
            "SOP-05_herg_retention_curve.csv": "hERG 阈值—保留曲线（0.01 步长）",
            "SOP-05_herg_threshold_analysis.json": "阈值锚定证据（含输入 sha256 与复现核对）",
            "SOP-05_HERG_THRESHOLD_SOP.md": "hERG 阈值选择 SOP（包内版本，含变更控制）",
            "SOP-05_MD_shortlist_scaffoldcap3.csv": "敏感性名单：骨架封顶 ≤3",
            "SOP-05_MD_shortlist_scaffoldcap5.csv": "敏感性名单：骨架封顶 ≤5",
            "SOP-05_MD_shortlist_sens_QED0.60_anchor4.csv": "敏感性名单：QED≥0.60 ∧ anchor=4",
            "SOP-05_MD_shortlist_sens_QED0.40_anchor4.csv": "敏感性名单：QED≥0.40 ∧ anchor=4",
            "SOP-05_MD_shortlist_sens_QED0.50_anchor3.csv": "敏感性名单：QED≥0.50 ∧ anchor≥3",
            "SOP-05_MD_shortlist_sens_QED0.60_anchor3.csv": "敏感性名单：QED≥0.60 ∧ anchor≥3",
            "SOP-05_MD_shortlist_sens_QED0.40_anchor3.csv": "敏感性名单：QED≥0.40 ∧ anchor≥3",
            "INDEX.md": "本索引",
            "MANIFEST.sha256": "交付文件 sha256（附 md5 对照）",
        }
        for nm in files:
            fn = nm.name
            if fn.startswith("SOP-05_fig_"):
                notes.setdefault(fn, "图（PNG+SVG，SVG 已去元数据以保证可复现）")
            if nm.suffix == ".py":
                notes.setdefault(fn, "本包脚本副本（便于交付目录内追溯生成逻辑）")
        rows_py, rows_data = [], []
        for f in sorted(files, key=lambda x: str(x.relative_to(DELIV))):
            if f.name == "acceptance_report.json":
                continue
            row = f"| `{f.relative_to(DELIV)}` | {notes.get(f.name, '见对应步骤文档')} |"
            (rows_py if f.suffix == ".py" else rows_data).append(row)
        idx += rows_data + ["", "脚本与辅助文件：", "", "| 文件 | 说明 |", "|---|---|"] + rows_py
        for name, note in notes.items():
            if not (DELIV / name).is_file():
                idx.append(f"| `{name}` | {note}（由后续步骤生成） |")
        write_text(DELIV / "INDEX.md", "\n".join(idx))
        log.info(f"交付层汇总完成：{len(files)} 个文件 → 06_deliverables/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
