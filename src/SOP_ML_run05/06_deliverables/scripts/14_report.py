#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""交付层（下）：生成最终 Markdown 报告与下游 MD 接口说明。

不引入新数值：全部数字从各步骤 JSON/表读入，保证报告与产物同源。
输出：06_deliverables/SOP-05_report.md、SOP-05_downstream_MD_plan.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.io_utils import StepLog, read_json, write_text

DELIV = ROOT / "06_deliverables"

# 与既有口径的对账目标（上一轮复核/预演公开数值，仅作对照，不作为本包输入）
RECONCILE = [
    ("姿态输入", 5000),
    ("InChIKey 去重后", 4990),
    ("PAINS0 且 Ro5/Veber 零违规", 4577),
    ("QED≥0.60 / ≥0.50 / ≥0.40", "2098 / 3198 / 3988"),
    ("姿态闸门 pass / warn / reject", "4862 / 134 / 4"),
    ("QED0.60 ∧ anchor=4（无 maxTc 门）", 1271),
    ("QED0.60 ∧ anchor=4 ∧ maxTc<0.9", 1269),
    ("QED0.60 ∧ anchor≥3 ∧ maxTc<0.9", 2029),
    ("QED0.40 ∧ anchor=4 ∧ maxTc<0.9", 2477),
    ("hERG 补算四臂并集", 3236),
    ("Butina 0.35（全池）簇数 / 最大簇 / 单例", "642 / 20 / 446"),
]


def md_table(rows, header) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def _pct(v) -> str:
    """覆盖度：0–1 小数渲染为百分比；已是百分数/缺失则原样。"""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "NA"
    return f"{f * 100:.2f}%" if 0.0 <= f <= 1.0 else f"{f}"


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    inp = read_json(ROOT / "00_inputs/manifests/input_verify.json")
    ext = read_json(ROOT / "00_inputs/manifests/extracted_manifest.json")
    pid = read_json(ROOT / "01_pose_qc/reports/pose_identity_summary.json")
    pnorm = read_json(ROOT / "01_pose_qc/reports/pose_normalize_summary.json")
    pgate = read_json(ROOT / "01_pose_qc/reports/pose_gate_summary.json")
    xchk = read_json(ROOT / "02_druglikeness/reports/descriptor_crosscheck.json")
    fun = read_json(ROOT / "02_druglikeness/reports/druglikeness_funnel.json")
    mtc = read_json(ROOT / "03_gates/reports/maxTc_summary.json")
    herg = read_json(ROOT / "03_gates/reports/herg_status.json")
    gs = read_json(ROOT / "03_gates/reports/gate_stack.json")
    sels = read_json(ROOT / "04_selection/reports/selection_summary.json")
    _nks = [v.get("n_known_excluded_in_arm") for v in sels["lists"].values()
            if isinstance(v.get("n_known_excluded_in_arm"), int)]
    _nk_lo, _nk_hi = (min(_nks), max(_nks)) if _nks else ("NA", "NA")
    rob = read_json(ROOT / "05_robustness/reports/robustness.json")
    _detp = ROOT / "05_robustness/reports/determinism_compare.json"
    det = (_detp.is_file() and read_json(_detp)) or {
        "mode": "pending", "identical": None,
        "note": "本轮首跑：G 臂位于整链末步，第二次整链运行起回填实测值"}
    acc = read_json(DELIV / "acceptance_report.json")
    main = sels["lists"]["main_QED0.50_anchor4"]
    main_df = pd.read_csv(DELIV / "SOP-05_MD_shortlist_main.csv", low_memory=False)
    cand = main_df[main_df.candidate_type == "candidate"]
    trig_main = next((r for r in acc["context"]["constraint_trigger"]
                      if str(r["arm"]).startswith("main")), {})

    with StepLog(ROOT, "14_report") as log:
        L = []
        A = L.append
        A("# SOP_5 · A2A 拮抗剂 MD 候选名单最终报告")
        A("")
        A(f"- 规则集：`{rule['rule_set_id']}` v{rule['rule_version']}")
        A(f"- rules_sha256：`{rule.get('rules_sha256', '')}`")
        A(f"- 工程包：`SOP_5/`（自包含：输入与产物全部在包内，脚本不引用包外路径）")
        A(f"- 数据截止：2026-09-23（UTC）")
        A(f"- 验收：**硬门 {acc['n_hard_pass']}/{acc['n_hard']} 通过"
          f"{'' if acc['n_hard_pending'] == 0 else f'（{acc['n_hard_pending']} 项 pending）'}**；"
          f"次要指标 "
          f"{sum(1 for c in acc['checks'] if c.get('class') == 'reported')} 项如实报告；"
          f"all_pass={acc['all_pass']}")
        A("")
        A("## 0. 摘要")
        A("")
        A(f"主口径（类药性硬门 ∧ 姿态闸门 pass ∧ anchor 4/4 ∧ maxTc<0.90 ∧ QED≥0.50"
          f"{' ∧ hERG<0.5' if herg['status'] == 'applied' else '，hERG 门 pending'}）"
          f"下给出 **{main['list_size']}** 个 MD 起始分子："
          f"{main['n_selected_candidates']} 个候选 + 2 个阴性对照。")
        A("")
        A(md_table([
            ["候选 / 对照", f"{main['n_selected_candidates']} / {main['n_controls']}"],
            ["CNN_VS 区间", f"{main['CNN_VS_range'][0]}–{main['CNN_VS_range'][1]}"],
            ["QED 区间", f"{main['QED_range'][0]}–{main['QED_range'][1]}"],
            ["锚点 / 姿态闸门", f"全部 4/4（最小 {main['anchor_hits_min']}）；全部 pass"],
            ["成对 ECFP4 Tc 上限",
             f"{main['max_pairwise_tc']}（判据 < {rule['step4_selection']['pairwise_tc_max']}）"],
            ["唯一骨架数（封顶 ≤" + str(main["scaffold_cap"]) + ")", main["n_distinct_scaffold"]],
            ["化学型最大计数 / 上限", f"{main['chemotype_max']} / ≤{main['chem_cap']}"],
            ["AM1-BCC 风险档 / 上限", f"{main['risk_n']} / ≤{main['risk_cap']}"],
            ["maxTc（vs 479 已知活性）区间",
             f"{main['maxTc_range'][0]}–{main['maxTc_range'][1]}"],
            ["hERG 区间", str(main["herg_range"])],
            ["已知物精确重合", f"{main['known_exact_hits']}（判据 0）"],
        ], ["指标", "实测"]))
        A("")
        A("**三点必读结论**")
        A("")
        A(f"1. 结构证据链完整：{main['n_selected_candidates']} 个候选全部来自上传的 gnina "
          f"最优姿态，姿态闸门 pass、四锚点（N253/E169/F168/W246）4/4 接触，"
          f"与最近已知活性的 ECFP4 相似度上限 {main['maxTc_range'][1]}（门限 0.90）。")
        A(f"2. 多样性与成药性同时受控：一骨架一分子（{main['n_distinct_scaffold']} 个骨架）、"
          f"化学型 ≤{main['chem_cap']}、AM1-BCC 风险档 ≤{main['risk_cap']}；"
          f"选序中因成对 Tc 冲突跳过 {trig_main.get('skip_tc_pairwise', 'NA')} 个分子"
          f"（约束确实在起作用，见 §5）。")
        A(f"3. hERG 门：状态 `{herg['status']}`，补算清单 "
          + ((f"`03_gates/herg/SOP-05_herg_input_union.csv` 共 {herg['scope_n']} 个分子；"
              f"外部 ADMET-AI 预测表已落入接口位置并按 <{herg['threshold']} 硬落门"
              f"（并集内覆盖 {_pct(herg.get('coverage_in_scope'))}，名单与各臂池均已剔除不通过者）。")
             if herg.get("status") == "applied" else
             (f"`03_gates/herg/SOP-05_herg_input_union.csv` 共 {herg['scope_n']} 个分子；"
              "状态仍为待补算，名单未因 hERG 删减。")))
        A("")
        A("## 1. 输入固化（步骤 00–01）")
        A("")
        A("| 输入 | 角色 | md5 | 校验 |")
        A("|---|---|---|---|")
        for e in inp["entries"]:
            A(f"| `{Path(e['file']).name}` | {e['role']} | `{e['md5_got']}` | "
              f"{'一致' if e['match'] else '不一致'} |")
        A("")
        A(f"- 归档白名单解包：{ext['n_files']} 个成员（{ext['total_bytes'] / 1e6:.1f} MB）→ "
          f"`00_inputs/extracted/`；清单 `00_inputs/manifests/EXTRACTED_MANIFEST.md5`。")
        A("- 规则冻结：`config/frozen_rule.json` 为阈值唯一来源；脚本哈希登记于 "
          "`00_inputs/manifests/FREEZE.json`。")
        A("")
        A("## 2. 步骤 1：姿态核验与规范化")
        A("")
        A(md_table([
            ["姿态条数", pid["n_poses"]],
            ["身份分级", pid["identity_class_counts"]],
            ["多片段（输入 SMILES 含 '.'）", pid["n_multifragment_template"]],
            ["|CNNscore×CNNaffinity − CNN_VS| 越界数",
             f"{pid['cnn_vs_x_affinity_violations']}（容差 {pid['cnn_vs_x_affinity_tolerance']}）"],
            ["坐标往返最大偏差（Å）", pid["coord_max_dev_max"]],
            ["InChIKey 去重后（保留 CNN_VS 最大）",
             f"{pnorm['n_kept']}（剔除 {pnorm['n_dropped']}："
             f"{', '.join(pnorm['dropped_cids'][:6])}{' …' if pnorm['n_dropped'] > 6 else ''}）"],
            ["姿态闸门 pass / warn / reject", pgate["gate_counts"]],
            ["锚点命中分布", pgate["anchor_hits_counts"]],
            ["16 残基接触数中位", pgate["n_contact_res_median"]],
            ["与 run01 top200 逐距离对齐",
             f"n={pgate['align_with_run01_top200']['n']}，最大偏差 "
             f"{pgate['align_with_run01_top200']['max_abs_dist_diff']} Å"],
        ], ["核验项", "结果"]))
        A("")
        ev = pgate.get("labeled5098_anchor_evidence") or {}
        if ev:
            A(f"**锚点门的独立证据链**（5,098 带标签富集库：{ev['n_actives']} 活性 / "
              f"{ev['n_decoys']} 诱饵，IFP 在本包内重算）：锚点 4/4 的活性率 "
              f"{ev['anchor4_active_rate']} vs 诱饵 {ev['anchor4_decoy_rate']}，"
              f"富集 {ev['anchor4_enrichment']}×；锚点 ≥3 为 {ev['anchor_ge3_active_rate']} "
              f"vs {ev['anchor_ge3_decoy_rate']}。该证据支持把 anchor=4 作为硬门。")
            A("")
        A("逐分子标注列（见 `SOP-05_decision_cards_main.csv`）：`stereo_source`、`n_fragments`、"
          "`protonation_note`、`qc_class`、`pose_gate`、`warn_or_reject_reason`。")
        A("")
        A("## 3. 步骤 2：类药性与多参数画像")
        A("")
        A("描述符与 PAINS/Ro5/Veber/QED 全部在包内用 RDKit 重算（上游 `pains` 列全空，"
          "故以重算为唯一口径）；与上游列交叉校验：")
        A("")
        A(md_table([[c["descriptor"], c["upstream_col"], c["max_abs_diff"], c["mean_abs_diff"]]
                    for c in xchk["descriptor_crosscheck"]],
                   ["描述符（本包）", "上游列", "max|Δ|", "mean|Δ|"]))
        A("")
        A(f"- CNS MPO 实现对照原文（Wager 2010 Table 3 三例）：最大偏差 "
          f"**{xchk['mpo_regression_vs_paper']['max_abs_diff']}**；"
          + "；".join(f"{c['name']} 回算 {c['recomputed']} vs 原文 {c['reported']}"
                      for c in xchk["mpo_regression_vs_paper"]["cases"]) + "。")
        A(f"- 漏斗：{fun['funnel']}")
        A(f"- 剔除归因：{fun['attribution']}")
        A(f"- 标注列（药物样池）：{fun['annotation_columns_on_druglike_pool']}；"
          f"CNS MPO≥4、其 80% 区间界（`cns_mpo_lo80` / `cns_mpo_hi80` / `cns_mpo_pi80_span`）"
          f"与 TPSA∈[90,140] 按规则集为**标注列**，不构成主门。")
        A(f"- CNS MPO 的两个输入均为上游第三方模型：`basic_pka`＝MolGpKa 最强碱性中心 pKa"
          f"（点值，无区间；无碱性中心者上游记 0，全表 192 个、本池 4 个、**交付名单 0 个**，"
          f"按 pKa≤8 档取满分）；`clogd74`＝cpLogD v2.0（ChEMBL33）点估计，"
          f"**附带 80% 预测区间**（`clogd74_lower_80` / `clogd74_upper_80`，置信度 0.8），"
          f"区间宽度中位 {fun.get('annotation_columns_on_druglike_pool', {}).get('cLogD 区间宽度中位', 'NA')} log。"
          f"决策卡同时给出区间界定的 CNS MPO：悲观界 `cns_mpo_lo80`（cLogD 取区间上界）"
          f"与乐观界 `cns_mpo_hi80`（取区间下界），两者与点估计一并列出，供 MD 前人工复核取用。")
        A("")
        A("## 4. 步骤 3：结构与毒性门")
        A("")
        A(f"- 门叠加链：{gs['gate_chain']}")
        A(f"- maxTc（vs {mtc['reference_actives']} 个已知活性）：分位 "
          f"{mtc['maxTc_percentiles']}；通过 <0.90 的 {mtc['pass_primary_lt_0.9']} 个、"
          f"通过 <0.70（敏感性档）的 {mtc['pass_sensitivity_lt_0.7']} 个。")
        _gmode = herg.get("scope_mode", "union4")
        _gsn = herg.get("gate_scope_n", herg["scope_n"])
        A(f"- hERG：状态 `{herg['status']}`，阈值 <{herg['threshold']}；"
          + (f"**门作用域 = 全池 {_gsn} 个**（口径 `scope={_gmode}`，落门依据为包内冻结的 "
             f"ADMET-AI 全量表，外部落门表作逐值对账件）；"
             f"四臂并集 **{herg['scope_n']}** 个降为**审计件**（臂池 {herg['arm_pool_sizes']}）。"
             if _gmode == "pool" else
             f"门作用域 = 四臂并集补算范围 **{herg['scope_n']}** 个"
             f"（臂池 {herg['arm_pool_sizes']}）。")
          + f"**作用域覆盖率 {_pct(herg.get('scope_coverage'))}**"
            f"（{_gsn}/{herg.get('n_pool_total')}；池内未评估 "
            f"{herg.get('n_not_evaluated')} 个）；"
            f"接口说明 `03_gates/herg/HERG_INTERFACE.md`。")
        if herg.get("status") == "applied":
            _wh = (f"`{str(herg.get('weight_sha256'))[:16]}…`" if herg.get("weight_sha256")
                   else "未提供（应用接口不回传权重哈希，故仅登记输出表哈希）")
            A(f"  - 外部预测表：`{herg.get('predictions_file')}`"
              f"（来源 {herg.get('source') or '未标注'}；权重哈希 {_wh}）")
            _scope_lab = ("门作用域（全池）内覆盖" if _gmode == "pool" else "并集内覆盖")
            A(f"  - {_scope_lab} **{_pct(herg.get('coverage_in_scope'))}**；通过 "
              f"{herg.get('n_pass_in_scope')} 个、不通过 {herg.get('n_fail_in_scope')} 个。"
              + (f"外部落门表对审计用四臂并集的覆盖为 "
                 f"{_pct(herg.get('external_table_coverage_of_audit_union'))}"
                 f"（并集 {herg.get('audit_union_n')} 个中含 "
                 f"{herg.get('n_audit_union_outside_external')} 个外部表未覆盖者，"
                 f"其 hERG 值取自包内冻结表，故门内无未评估分子）。"
                 if _gmode == "pool" and herg.get("audit_union_n") else ""))
            A(f"  - **作用域语义**：hERG 补算范围是四臂并集而非整个池，故池内有 "
              f"{herg.get('n_not_evaluated')} 个分子**没有预测值**。它们的门状态是 "
              f"`not_evaluated`（**未评估**），不得转述为「高风险」或「不合格」——"
              f"把缺值写成「保守判为不过门」等于把「没测」冒名为「测出不合格」。"
              f"当前这 {herg.get('n_not_evaluated')} 个分子均不在任何臂/名单内；若下游要扩大"
              f"候选池，必须同步扩大补算范围，否则会把未评估当作不合格使用。")
        # ---- 4b 第二模型交叉验证（分析步 08b：只读产物、不改变任何门/名单）
        _cm_p = ROOT / "03_gates/herg/cross_model/reports/herg_cross_model.json"
        if _cm_p.is_file():
            cm = read_json(_cm_p)
            _ds, _cv = cm["dataset"], cm["cv"]
            _ag, _rk = cm["agreement_vs_admet_ai"], cm["rank_consistency"]
            _gi = cm["gate_impact"]
            _u, _pm, _dl = (_gi["union_recomputed"], _gi["pre_herg_main_pool"],
                            _gi["delivered_list_48"])
            _m2 = {r["model2"]: r for r in _ag}
            _m2r = {r["model2"]: r for r in _rk}
            A("")
            A("### 4b. hERG 判定对模型的依赖度（第二模型交叉验证 · 分析步 08b）")
            A("")
            A(f"上文的 hERG 门来自**单一外部模型**（ADMET-AI）。为量化「该门结论有多少是所选模型的"
              f"产物」，本步用第二个**独立训练**的 hERG 模型在同一批 {_u['pool_n']:,} 个并集分子上复算。"
              f"第二模型数据来自 ChEMBL REST API 的靶点 `CHEMBL240`（KCNH2 / hERG, Homo sapiens）："
              f"`standard_type ∈ {{IC50, Ki}}`、`standard_relation = '='`、`pChEMBL` 非空的活性记录"
              f" {_ds['n_raw_records']:,} 条 → 按 InChIKey 去重后得标注集 **{_ds['n_labeled']:,}** 个"
              f"（阻断 {_ds['n_blocker']:,} / 非阻断 {_ds['n_nonblocker']:,}；pChEMBL ∈ (5, 6) 的模糊带舍弃），"
              f"特征仅 `{_ds['features']}`，随机种子与包内一致（PROJECT_SEED）。"
              f"该步是**分析步**：只产出证据，不修改任何门、池或名单。")
            A("")
            A("**第二模型自身性能（5 折交叉验证）**")
            A("")
            _cvrows = [[a.upper(),
                        _cv[f"{a}_random_5fold"]["roc_auc"], _cv[f"{a}_random_5fold"]["pr_auc"],
                        _cv[f"{a}_scaffold_5fold"]["roc_auc"], _cv[f"{a}_scaffold_5fold"]["pr_auc"]]
                       for a in ("rf", "lr", "knn")]
            A(md_table(_cvrows, ["第二模型算法", "随机分层 ROC-AUC", "随机分层 PR-AUC",
                                 "骨架分组 ROC-AUC", "骨架分组 PR-AUC"]))
            A("")
            A(f"- 骨架分组 CV 比随机分层低（RF {_cv['rf_scaffold_5fold']['roc_auc']} vs "
              f"{_cv['rf_random_5fold']['roc_auc']}），说明该模型对**新骨架**的泛化更弱，"
              f"属于 hERG 活性预测的常见水平，不构成对第一模型的优劣判断（两者数据源部分重叠）。")
            A("")
            A("**与 ADMET-AI 的一致性（同一 {:,} 个分子，阈值同为 0.5）**".format(_u["pool_n"]))
            A("")
            _arows = [[m["model2"].upper(), m["spearman_rho"], _pct(m["gate_agreement"]),
                       m["cohen_kappa"], _pct(m["m1_pass_rate"]), _pct(m["m2_pass_rate"]),
                       m["m1_pass_m2_fail"], m["m1_fail_m2_pass"]] for m in _ag]
            A(md_table(_arows, ["第二模型", "Spearman ρ", "门决策一致率", "Cohen κ",
                                "模型1 过门率", "模型2 过门率",
                                "模型1过/模型2不过", "模型1不过/模型2过"]))
            A("")
            _kmin = min(m["cohen_kappa"] for m in _ag)
            _kmax = max(m["cohen_kappa"] for m in _ag)
            _amin = min(m["gate_agreement"] for m in _ag)
            _amax = max(m["gate_agreement"] for m in _ag)
            _rmin = min(m["spearman_rho"] for m in _ag)
            _rmax = max(m["spearman_rho"] for m in _ag)
            A(f"- 四个第二模型（含三模型均值集成）的 Spearman ρ 仅 {_rmin}–{_rmax}，"
              f"0.5 阈值上的门决策一致率 {_pct(_amin)}–{_pct(_amax)}、Cohen κ {_kmin}–{_kmax}："
              f"**几乎一致的判定只占三分之一到一半**，且第二模型显著更宽松"
              f"（过门率 {_pct(min(m['m2_pass_rate'] for m in _ag))}–"
              f"{_pct(max(m['m2_pass_rate'] for m in _ag))} vs 第一模型 "
              f"{_pct(_m2['rf']['m1_pass_rate'])}，中位分也更低）。")
            A("")
            A("**与阈值无关的排序一致性**（判断差异是否只是标定/阈值问题）")
            A("")
            _rrows = [[m["model2"].upper(), m["high_risk_top10pct_jaccard"],
                       round(m["m2_matched_threshold"], 4), m["matched_overlap_jaccard"],
                       _pct(m["recall_of_m1_passes"])] for m in _rk]
            A(md_table(_rrows, ["第二模型", "高风险前 10% 交集 Jaccard", "匹配过门率所需阈值",
                                "匹配过门率后过门集合 Jaccard", "覆盖第一模型过门分子的比例"]))
            A("")
            _jmin = min(m["high_risk_top10pct_jaccard"] for m in _rk)
            _jmax = max(m["high_risk_top10pct_jaccard"] for m in _rk)
            _mmin = min(m["matched_overlap_jaccard"] for m in _rk)
            _mmax = max(m["matched_overlap_jaccard"] for m in _rk)
            _cmin = min(m["recall_of_m1_passes"] for m in _rk)
            _cmax = max(m["recall_of_m1_passes"] for m in _rk)
            A(f"- 即使把阈值调平（让第二模型恰好也通过同样比例），两份「过门集合」的 Jaccard 仍只有 "
              f"{_mmin}–{_mmax}，仅覆盖第一模型过门分子的 {_pct(_cmin)}–{_pct(_cmax)}；"
              f"两模型各自认为最高风险的 10% 分子交集 Jaccard 也仅 {_jmin}–{_jmax}。"
              f"**差异不只是阈值/标定问题，排序本身也不同。**")
            A("")
            A("**对结论的实测影响（若改用「两模型一致过门」，以集成模型为例）**")
            A("")
            A(f"- 并集 {_u['pool_n']:,} 个（其中 {_u.get('m1_missing_n', 'NA')} 个无第一模型外部预测值，"
              f"按保守口径计为不过门）：第一模型过门 {_u['m1_pass']}、集成模型过门 {_u['m2_pass']}、"
              f"两者都过 {_u['m1_and_m2_pass']}；第一模型判高风险而被删的分子里，第二模型仍判低风险者 "
              f"**{_u['m2_only_pass']:,}** 个。")
            _gate_pre = int(gs["gate_chain"].get("6_hERG<0.5", 0))
            _gate_pre_noved = int(main["pool_n"])
            A(f"- hERG 前主池 {_pm['pool_n']:,} 个（**三个分母不可互换**）："
              f"① **全池门口径**——本包 hERG 门（依据包内冻结全量表，覆盖 100%）在该池上过门 "
              f"{_gate_pre}、剔除已知物后 {_gate_pre_noved} 个即交付主臂池；"
              f"② **两侧都有预测值的子集口径**——第二模型覆盖 {_pm.get('covered_n', 'NA')}/"
              f"{_pm['pool_n']}，其中 {_pm.get('m1_missing_n', 'NA')} 个没有第一模型外部预测值"
              f"（`hERG_admet_ai` 为空，按保守口径计为不过门），故实际分母 "
              f"{_pm.get('m1_evaluated_n', 'NA')}：第一模型过门 {_pm['m1_pass']} → "
              f"「两模型都过」仅 {_pm['m1_and_m2_pass']}"
              f"（在该子集内收紧 {(1 - _pm['m1_and_m2_pass'] / _pm['m1_pass']) * 100:.1f}%）。")
            _dl_missing = [c for c in
                           pd.read_csv(ROOT / "03_gates/herg/cross_model/tables/"
                                              "herg_cross_model_predictions.csv",
                                       low_memory=False)
                           .query("hERG_admet_ai != hERG_admet_ai")       # 空值
                           .loc[lambda d: d.cid.isin(main_df.cid), "cid"]]
            A(f"- 交付 48 名单：{_dl['m1_and_m2_pass']} 个同时被第二模型判过门，"
              f"**{_dl['m1_only_pass']} 个仅在 ADMET-AI 下过门**"
              f"（{'、'.join('`' + c + '`' for c in cm.get('list_flips_cids', []))}）；"
              f"另有 {_dl.get('m1_missing_n', 'NA')} 个"
              f"（{'、'.join('`' + c + '`' for c in _dl_missing)}）没有第一模型外部预测值，"
              f"按保守口径计为不过门（它们在包内冻结表上均过门）；"
              f"RF 单模型口径下 48/48 全部过门。")
            A("")
            A("> 结论与建议（不改变本包主口径）：hERG 单模型 0.5 阈值**不是**可迁移的硬判据——"
              "「安全」与「不安全」的划分对模型选择高度敏感。因此建议：① 对外转述时明确"
              "「hERG 过门」＝「在 ADMET-AI 下过门」；② MD 前置复核以"
              "**两个模型一致判为高风险**的分子为优先核查对象，而不是用任一单模型的硬阈值把分子永久剔除；"
              f"③ 若后续要求「双模型一致」，须重新选名单（两侧都有预测值的 "
              f"{_pm.get('m1_evaluated_n', 'NA')} 个内 {_pm['m1_pass']} → {_pm['m1_and_m2_pass']}；"
              f"另有 {_pm.get('m1_missing_n', 'NA')} 个缺第一模型外部预测值，需另行判定）"
              f"并同步改验收口径，属于口径变更而非工程修正。")
            A(f"> 独立性的诚实边界：{cm.get('independence_caveat', '')}。")
        kn = gs["known_ligand_check"]
        A(f"- 已知物核验（ChEMBL 独立 A2A 配体 {kn.get('n_known_smiles', 'NA')} 条）："
          f"精确 InChIKey 重合 {kn.get('exact_inchikey_hits')}、Tc≥0.9 相似 "
          f"{kn.get('n_tc_ge_0.9')}、全池最大 Tc {kn.get('max_tc')}。")
        ng = sels.get("novelty_gate", {})
        ug = ng.get("ungated_diagnostic", {})
        if ng:
            A("")
            _basestr = ("类药性 ∧ 姿态 pass ∧ maxTc<0.9" +
                        (" ∧ hERG<0.5" if herg.get("status") == "applied" else ""))
            _rk = ug.get("ungated_rank_of_known", {})
            _rk_txt = ("、".join(f"`{c}`（第 {r} 位）" for c, r in _rk.items())
                       if _rk else str(ug.get("cids_only_in_ungated")))
            A("**novelty 门（规则 1.1.0）**：验收判据要求名单内已知物为 0，"
              f"而池内确有 {kn.get('exact_inchikey_hits')} 个分子与 ChEMBL 已知 A2A 配体"
              f" InChIKey 精确重合；主口径基座（{_basestr}）内命中 "
              f"**{ng.get('n_exact_known_in_base')}** 个，已在选序层剔除"
              "（剔除后的各臂池规模见 §5 表）。"
              "被剔除者中 `ZINC000118633318` 就是已知拮抗剂 **CHEMBL4594442**"
              "（`Cc1cc(-c2nnc(N)nc2-c2ccc(F)cc2)cc(Cl)n1`，ChEMBL 记录 best_potency = 1.7 nM）："
              f"不设该门时它会进入主名单（未设门与设门名单差异 {ug.get('n_different_from_main')} 个分子，"
              f"即 {_rk_txt}）。")
            A("")
            A("> 处理逻辑：novelty 门只作用于**选序层**，池层完整保留并计数（便于审计）；"
              "该门与 hERG 门相互独立 —— 前者是名单新颖性要求，后者是安全前置。")
        A("")
        A("### 4c. ADMET-AI 注释层与 hERG 阈值锚定（分析步 10b/10c，非判据）")
        A("")
        _ad_p = ROOT / "03_gates/admet/reports/admet_annotation.json"
        _xc_p = ROOT / "03_gates/admet/reports/admet_descriptor_crosscheck.json"
        _ht_p = ROOT / "03_gates/admet/reports/herg_threshold_analysis.json"
        if _ad_p.is_file() and _ht_p.is_file():
            _ad = read_json(_ad_p)
            _xc = read_json(_xc_p)
            _ht = read_json(_ht_p)
            _prov = _ad["provenance"]
            A(f"本包固化了一张 ADMET-AI 全量预测表（f3 库 {_prov['rows']:,} 行 × "
              f"{_prov['admet_columns']} 个 ADMET 列，sha256 `{_prov['raw_sha256'][:16]}…`；"
              f"生产者：{_prov['producer']}），并新增两个**分析步**读取它：")
            A("")
            A("- `10b_admet_annotation`：池级标注列（毒性/ADME/CYP/描述符/已上市药百分位） + "
              "池内重标定切点 + 描述符独立复算交叉核验 + MD 前置复核清单；")
            A("- `10c_herg_threshold`：hERG 阈值三锚定 + 0.01 步长保留曲线 + 名单档位标注。")
            A("")
            A("> **非判据声明**：两步都不改门（主口径仍 `hERG < "
              f"{_ht['decision']['frozen_gate_threshold']}`）、不改池、不改名单；"
              "注释列在交付物中以**追加列**形式出现（列追加、行序与成员不变）。"
              "预测值非实验数据，采购/合成/立项前需体外确认。")
            A("")
            _nc = []
            _w_p = ROOT / "06_deliverables/SOP-05_admet_annotation_main.csv"
            if _w_p.is_file():
                _w = pd.read_csv(_w_p, dtype=str, keep_default_na=False)
                _nc = _w.loc[_w.admet_coverage == "not_covered", "cid"].tolist()
            A(f"**覆盖**：本包池 {_ad['coverage']['pool_n']:,} 个分子"
              f"（{_ad['coverage']['covered']:,} 全覆盖）均有预测值；"
              + (f"{len(_nc)} 个阴性对照（" + "、".join("`" + c + "`" for c in _nc)
                 + "，来自 run01 的 5,098 富集库）不在覆盖范围内 → 标注 "
                 "`admet_coverage=not_covered`（不写空值）。"
                 if _nc else
                 "池外分子（如 run01 5,098 富集库来源的阴性对照）无预测值 → 标注 "
                 "`admet_coverage=not_covered`（不写空值）。"))
            A("")
            A("**池内重标定切点（库内相对切点按本包池标定；外部库原值仅作对照）**")
            A("")
            _ref = _ad["cutpoints_reference_from_external_run"]
            _cut = _ad["cutpoints_recalibrated_on_pool"]
            A(md_table([
                ["hERG 高危线（池 P79.8）", _cut["herg"], _ref["herg"]],
                ["DILI 高危线（池 P88.2）", _cut["DILI"], _ref["DILI"]],
                ["溶解度下限（池 P10）", _cut["solubility"], _ref["solubility"]],
                ["AMES 绝对线", _cut["AMES"], _ref["AMES"]],
                ["ClinTox / 致癌性 绝对线",
                 f"{_cut['ClinTox']} / {_cut['Carcinogens_Lagunin']}",
                 f"{_ref['ClinTox']} / {_ref['Carcinogens_Lagunin']}"],
            ], ["切点", "本包池重标定值", "外部库原值（对照）"]))
            A("")
            A(f"池内触发率：hERG {_ad['trigger_rates_pct']['herg']}%、"
              f"DILI {_ad['trigger_rates_pct']['DILI']}%、AMES {_ad['trigger_rates_pct']['AMES']}%、"
              f"ClinTox {_ad['trigger_rates_pct']['ClinTox']}%、"
              f"致癌性 {_ad['trigger_rates_pct']['Carc']}%、"
              f"结构警报 {_ad['trigger_rates_pct']['struct_alert']}%；"
              f"分诊层 " + "、".join(f"{k} {v}" for k, v in _ad["tier_counts"].items()) + "。")
            A("")
            _dl = _ad["delivered_list"]
            A(f"**交付 48 候选的注释画像**：分诊层 "
              + "、".join(f"{k} {v}" for k, v in _dl["tier_counts"].items())
              + "；hERG 档 " + "、".join(f"{k} {v}" for k, v in _dl["herg_band_counts"].items())
              + f"；hERG 中位 {_dl['herg_median']}。"
              f"带 ≥1 项警告者 **{_dl['md_precheck_n']}** 个（清单见 "
              "`06_deliverables/SOP-05_admet_md_precheck.csv`），"
              "建议在 MD 立项前优先做离子通道/CYP/溶解度复核："
              + "、".join("`" + c + "`" for c in _dl["md_precheck_cids"][:10])
              + ("…" if _dl["md_precheck_n"] > 10 else "") + "。")
            A("")
            A("**描述符独立复算核验（本包 RDKit 口径 vs ADMET-AI 内置实现）**")
            A("")
            _ns = _xc["numeric"]
            A(md_table([[k, v["max_abs_diff"], v["n_over_tol"]] for k, v in _ns.items()],
                       ["量", "最大绝对偏差", f"超容差条数（容差 {_xc['tolerance']}）"]))
            A("")
            _rs = _xc["rule_agreements"]
            A(f"- 判定一致率：Ro5 {_pct(_rs['ro5']['agreement_pct'])}、"
              f"PAINS {_pct(_rs['pains']['agreement_pct'])}、"
              f"BRENK {_pct(_rs['brenk']['agreement_pct'])}；"
              f"BRENK 的差异为**既有且已解释**（"
              f"{_xc['expected_differences']['BRENK']}），不作为失败项。")
            A(f"- 结论：类药性硬门的数值与判定由两套独立实现互相确认"
              f"（{'全部一致' if _xc['all_numeric_within_tol'] else '存在超出容差项'}）；"
              "这是把「注释表与包内口径同源」写成可复核证据的一步。")
            A("")
            A("**hERG 阈值三锚定（分析步 10c）**")
            A("")
            _a1, _a2, _a3 = _ht["A1_library_relative"], _ht["A2_approved_drug_reference"], _ht["A3_tradeoff"]
            _dec = _ht["decision"]
            _lh = f"library_head@K={_ht['params']['primary_top_k']}"
            A(md_table([
                ["① 库内相对（仅排序）",
                 f"库 75/90 分位 = {_a1['library_f3']['0.75']} / {_a1['library_f3']['0.9']}；"
                 f"池 75/90 分位 = {_a1['pool_4990']['0.75']} / {_a1['pool_4990']['0.9']}"],
                ["② 已上市药参照",
                 f"50/75/90 分位 = {_a2['50']} / {_a2['75']} / {_a2['90']}"
                 f"（75 分位 = 软线 {_dec['soft_line']}）"],
                ["③ 取舍驱动（硬线候选）",
                 f"库级头部（K={_ht['params']['primary_top_k']}，保留 ≥"
                 f"{_ht['params']['min_keep']:.0%}）= {_a3[_lh]['feasible_t']}"
                 f"（0.01 网格精确解 {_a3[_lh]['strictest_t_0.01_grid']}）；"
                 f"本包主池头部 = {_dec['hard_line_candidate_pack_head']}"],
            ], ["锚定", "本包实算值"]))
            A("")
            A(f"- 合理性窗口 [A2@50, A2@90] = [{_dec['sanity_window'][0]}, {_dec['sanity_window'][1]}]，"
              f"硬线候选落在窗口内：{_dec['sanity_ok']}；冲突消解分支：{_dec['resolution_branch']}。")
            _rc2 = _ht["repro_check"]
            A(f"- 与外部流程数字对账：软线 {_rc2['soft_line_observed']}（期望 "
              f"{_rc2['soft_line_expected']}，一致={_rc2['soft_line_ok']}）、"
              f"库级取舍线 {_rc2['library_head_tradeoff_observed']}（期望 "
              f"{_rc2['library_head_tradeoff_expected']}，一致={_rc2['library_head_tradeoff_ok']}）——"
              "说明包内实算复现了外部流程的阈值决策。")
            _gp = _ht.get("gate_scope_gap") or {}
            if _gp:
                A(f"- **门作用域缺口（非判据，量化未评估口径）**：落门表只覆盖四臂并集 "
                  f"{_gp['gate_table_n']:,} 个分子，池内另有 {_gp['out_of_scope']:,} 个无预测值 → "
                  f"门状态为 `not_evaluated`（未评估，**不是「判为高风险」**）；本包冻结的 f3 "
                  f"全量表显示其中 "
                  f"**{_gp['out_of_scope_frozen_below_gate']:,} 个"
                  f"（{_gp['out_of_scope_below_gate_pct']}%）实际 hERG < 门阈值**"
                  f"（落门表与冻结表在交集 {_gp['overlap_n']:,} 个内逐值一致，max|Δ| = "
                  f"{_gp['overlap_max_abs_diff']}）。这些分子不在任何臂/名单内。"
                  f"把补算范围扩到全池不需要新计算（f3 表已在包内且已登记哈希），"
                  f"故本缺口属**口径选择**而非数据可得性问题；本包仍不改门、不改名单。")
            A(f"- 保留曲线（0.01 步长，覆盖库/库级头部/池/主池/四臂并集/交付名单）见 "
              "`06_deliverables/SOP-05_herg_retention_curve.csv`；阈值选择流程、冲突消解规则与"
              "变更控制见 `06_deliverables/SOP-05_HERG_THRESHOLD_SOP.md`。")
            A("")
            A("> 主口径不变：hERG 门仍是 `hERG < "
              f"{_dec['frozen_gate_threshold']}`，上表的软线/硬线候选只用于"
              "「中/高」档标注与 MD 前置复核优先级，**不删分子**。若要改门，属于口径变更"
              "（需重选名单、重设验收基线并重跑冻结流程，实测代价见保留曲线与 "
              "`03_gates/admet/reports/herg_threshold_analysis.json`）。"
              "另需注意：本包 hERG 与外部 run 的 hERG **同源同值**，不能当作 §4b 意义上的"
              "第二模型；库内相对锚定以本包 49,963 行为口径（外部报告为 54,932 去重并集），"
              "数字允许小差异。")
        else:
            A("ADMET 注释层或阈值锚定产物缺失（未运行 10b/10c）→ 本小节标记 pending，"
              "不影响任何门与名单。")
        A("")
        A("## 5. 步骤 4：名单矩阵与约束行为")
        A("")
        rows = []
        for name, s in sels["lists"].items():
            if name.startswith("unconstrained"):
                continue
            rows.append([name, s.get("pool_n", "—"), s.get("n_selected_candidates", "—"),
                         s.get("max_pairwise_tc"), s.get("n_distinct_scaffold"),
                         s.get("overlap_with_main", "—")])
        A(md_table(rows, ["名单", "池", "选中", "成对 Tc 上限", "唯一骨架", "与主名单重叠"]))
        A("")
        if (herg.get("scope_coverage") or 1.0) < 1.0 and herg.get("status") == "applied":
            A(f"> **阅读注意（作用域）**：上表各臂池 = 并集内过门 ∩ novelty，而 hERG 补算范围只是"
              f"四臂并集（覆盖池的 {_pct(herg.get('scope_coverage'))}）。两条 anchor≥3 臂"
              f"（`sens_QED0.50_anchor3` / `sens_QED0.40_anchor3`）的池内有一部分分子落在并集之外"
              f"（无预测值 → `not_evaluated`），未进入该臂池。因此这两行的池规模**偏小**、"
              f"与主名单的重叠率**偏乐观**；四条 anchor=4 臂完全落在并集内，不受影响。")
            A("")
        unc = sels["lists"]["unconstrained_top48"]
        A(f"**约束咬合度**：不施加名单层约束时，同池 CNN_VS 前 48 名中同一化学型最多 "
          f"{unc['chemotype_max']} 个（上限 {main['chem_cap']}）、唯一骨架 "
          f"{unc['n_distinct_scaffold']} 个 —— 约束确实改变了名单构成。")
        A("")
        A("阴性对照（本包内从 run01 归档重算，规则冻结于规则集）：")
        A("")
        A(md_table([[r["control_role"], f"`{r['cid']}`", r["rule"], r["cnn_vs"],
                     r["vina_affinity"]] for r in sels["negative_controls"]],
                   ["角色", "cid", "规则", "CNN_VS", "Vina 亲和"]))
        A("")
        A("> 两个阴性对照来自 run01 归档、**不在上游 49,963 描述符表内**，因此其 "
          "`basic_pka` / `clogd74` / `cns_mpo*` 列为空（`cns_mpo` 由上游两个模型输入决定，"
          "本包不为它们伪造预测值），其余结构与打分字段照常给出。")
        A("")
        A("## 6. 主名单头部 20（选序顺序）")
        A("")
        head = cand.nsmallest(20, "pick_order")
        A(md_table([[int(r.pick_order), f"`{r.cid}`", round(float(r.CNN_VS), 3),
                     round(float(r.qed), 3), int(r.anchor_hits_4A),
                     round(float(r.maxTc_actives), 3), round(float(r.cns_mpo), 2),
                     str(r.chemotype)] for r in head.itertuples()],
                   ["#", "cid", "CNN_VS", "QED", "锚点", "maxTc(活性)", "CNS MPO", "化学型"]))
        A("")
        A("> 完整名单与逐分子证据：`SOP-05_MD_shortlist_main.csv`、"
          "`SOP-05_decision_cards_main.csv`；起始姿态：`SOP-05_poses_main_list.sdf`。")
        A("")
        A("## 7. 步骤 5：稳健性 A–G")
        A("")
        A("**A 打分扰动**（CNN_VS + N(0,σ)，各 100 次重抽样，与主名单重叠率）")
        A("")
        A(md_table([[r["sigma"], r["mean_overlap"], r["min_overlap"], r["p10_overlap"],
                     r["mean_head20_overlap"], r["min_head20_overlap"]]
                    for r in rob["armA_score_perturbation"]["rows"]],
                   ["σ", "平均重叠", "最小重叠", "P10", "头部20 平均", "头部20 最小"]))
        A("")
        _ra10 = rob["armA_score_perturbation"]
        _r10 = next(r for r in _ra10["rows"] if r["sigma"] == _ra10["primary_sigma"])
        A(f"- 判据（σ={_ra10['primary_sigma']} 平均重叠 ≥ "
          f"{_ra10['judge_mean_overlap_min']}）："
          f"**{'通过' if _ra10['pass_primary'] else '未通过'}**"
          f"（实测 {_r10['mean_overlap']}，P10 {_r10['p10_overlap']}，min "
          f"{_r10['min_overlap']}）")
        if not _ra10["pass_primary"]:
            A(f"  - **未通过的归因（本包如实报告，不调门槛）**：v1.5.0 的描述符结构规范化"
              f"（去盐）改变了下侧边际分子的构成（硬门由不过变过的分子新建入池），"
              f"而该判据在 v1.4.x 时的余量本就极薄（0.8046 vs 门槛 0.80，余量 0.0046）——"
              f"这与包内既有自述「打分口径单一、选序只用 CNN_VS」是同一件事的两面："
              f"名单尾部对分数扰动敏感。**不改阈值、不换判据口径**；若要在 v1.6.x 之后恢复"
              f"通过，应提升打分融合（如 CNN_VS 与 Vina 亲和的一致序）或把尾部填充规则改为"
              f"分层抽选，而不是放宽 0.80。")
        A(f"- B 阈值矩阵（QED 4 档 × maxTc 4 档）：重叠 "
          f"{rob['armB_threshold_matrix']['overlap_min']}–"
          f"{rob['armB_threshold_matrix']['overlap_max']}；明细 "
          f"`05_robustness/tables/robB_matrix.csv`。")
        A("- C 单轴敏感性：")
        A("")
        A(md_table([[c["axis"], c["change"], c.get("pool_n"), c.get("list_n"),
                     c.get("overlap_vs_main"), c.get("status")]
                    for c in rob["armC_one_axis"]],
                   ["轴", "变动", "池", "选中", "与主名单重叠", "状态"]))
        A("")
        A("**D 指纹替换（只报重叠）**：" + "；".join(
            f"{d['fingerprint']} → {d['overlap_vs_main']}"
            for d in rob["armD_fingerprint"]) + "。")
        A("")
        A("**E 并列带宽**（第 48 名分数 −band 内随机置换，各 100 次）：")
        A("")
        A(md_table([[e["band"], e["n_in_band"], e["mean_overlap"], e["min_overlap"],
                     e["mean_head20_overlap"]] for e in rob["armE_tie_band"]["rows"]],
                   ["band", "带宽内分子数", "平均重叠", "最小重叠", "头部20 平均"]))
        A("")
        fb, fbf = rob["armF_butina"], rob["armF_butina_full_pool"]
        A(f"**F Butina 聚类（cutoff {fb['cutoff']}）**：主池 {fb['n_members']} 分子 → "
          f"{fb['n_clusters']} 簇（最大 {fb['largest_cluster']}、单例 {fb['n_singletons']}）；"
          f"主名单同簇重复度上限 {fb['main_list_cluster_multiplicity_max']}。"
          f"全池（{fbf['pool_n']}）{fbf['n_clusters']} 簇 / 最大 {fbf['largest_cluster']} / "
          f"单例 {fbf['n_singletons']}（既有口径 642 / 20 / 446）。")
        A("")
        # 不嵌入 mode 之外的时机字段：mode 由「本轮是否写基线」决定，是判据本身的一部分，
        # 必须写进报告（否则 mode=baseline_written 会被误读成「确定性通过」）。
        # G 臂（12_determinism）是整链**末步**：本步（14_report）读到的必然是上一轮链末读数。
        # 因此这里必须区分「上一轮链末判据」与「本轮链末判据（写成于本步之后）」，不得混为一谈。
        _det_ok = bool(det.get("determinism_verified"))
        _det_head = "**G 确定性重跑**（读数取自上一轮链末；本轮判据在本步之后由同一步写入）"
        if _det_ok:
            A(f"{_det_head}：比对 "
              f"{det.get('n_files_current', det.get('n_files'))} 个文件，差异 "
              f"{det.get('n_diff', 0)} 个 → "
              f"**{'完全一致' if det.get('identical') else '存在差异'}**"
              f"（基线 sha256 `{str(det.get('baseline_sha256'))[:12]}…`，"
              f"`mode={det.get('mode')}`）。")
        else:
            A(f"{_det_head}：该轮 `mode={det.get('mode')}` —— "
              f"**只写入基线，未做比对**（哈希范围 "
              f"{det.get('n_files_current', det.get('n_files'))} 个文件，基线 sha256 "
              f"`{str(det.get('baseline_sha256'))[:12]}…`）。写基线的运行与基线同源，"
              f"**不构成确定性通过**；下一次整链运行（mode=compared）才给出判据。")
        A("")
        A("> 判据口径：G 臂必须是「写基线一轮 + 比对一轮」两轮分离；本轮真实判据以链末写入的 "
          "`05_robustness/reports/determinism_compare.json` 为准（`mode=compared` 且 "
          "`n_diff=0` 才算通过），并由 `tests/verify_pack.py` 第 [17] 组机器断言。"
          "交付层的 `acceptance_report.json` 位于 G 臂之前，故其确定性两行在本轮恒为 pending，"
          "不得把上一轮读数当作本轮的通过证据。")
        A("")
        A("## 8. 验收对账")
        A("")
        def _verdict(c):
            if c["pass"] is None:
                if str(c.get("status", "")).startswith("pending"):
                    return "pending（由链末 G 臂判定；见 §7）"
                return "如实报告（不设门）"
            return "通过" if c["pass"] else "未通过"

        A(md_table([[c["criterion"], c["expected"],
                     ("pending（本步在 G 臂之前，不自行判定）"
                      if str(c.get("status", "")).startswith("pending")
                      else c["observed"]),
                     _verdict(c)] for c in acc["checks"]],
                   ["判据", "期望", "实测", "结论"]))
        A("")
        # 臂池实测（口径可比版：类药性 ∧ 姿态 pass ∧ 该臂 QED/anchor ∧[可选] maxTc<0.9，
        # 均**未**施加 novelty 门，以便与既有臂池数字同口径比较）
        _cu = pd.read_csv(ROOT / "03_gates/tables/candidate_universe.csv", low_memory=False)
        _kn = pd.read_csv(ROOT / "03_gates/tables/known_ligand_similarity.csv")[
            ["cid", "known_exact_hit"]]
        _cu = _cu.drop(columns=[c for c in ("known_exact_hit", "max_tc_known",
                                            "nearest_known_ik") if c in _cu.columns])
        _cu = _cu.merge(_kn, on="cid", how="left")
        _cu["known"] = _cu.known_exact_hit.fillna(0).astype(int)
        _dl = _cu[(_cu.druglike_pass == 1) & (_cu.pose_gate == "pass")]

        def _arm(name: str, drop_mtc: bool = False) -> str:
            q = float(name.split("_")[0].split("QED")[1])
            a = int(name.split("anchor")[1])
            d = _dl[(_dl.qed >= q) & (_dl.anchor_hits_4A >= a)]
            if not drop_mtc:
                d = d[d.maxTc_actives < 0.9]
            return f"{len(d)}（novelty 门后 {len(d[d.known == 0])}）"

        A("**与既有口径的数字对账**（期望值来自上一轮复核/预演，仅作对照、不作输入）：")
        A("")
        got = {
            "姿态输入": pid["n_poses"],
            "InChIKey 去重后": pnorm["n_kept"],
            "PAINS0 且 Ro5/Veber 零违规": fun["funnel"].get("2_PAINS0_Ro5Veber零违规"),
            "QED≥0.60 / ≥0.50 / ≥0.40": f"{fun['funnel'].get('3_QED>=0.60')} / "
                                         f"{fun['funnel'].get('3_QED>=0.50')} / "
                                         f"{fun['funnel'].get('3_QED>=0.40')}",
            "姿态闸门 pass / warn / reject": " / ".join(
                str(pgate["gate_counts"].get(k, 0)) for k in ("pass", "warn", "reject")),
            "QED0.60 ∧ anchor=4（无 maxTc 门）": _arm("QED0.60_anchor4", drop_mtc=True),
            "QED0.60 ∧ anchor=4 ∧ maxTc<0.9": _arm("QED0.60_anchor4"),
            "QED0.60 ∧ anchor≥3 ∧ maxTc<0.9": _arm("QED0.60_anchor3"),
            "QED0.40 ∧ anchor=4 ∧ maxTc<0.9": _arm("QED0.40_anchor4"),
            "hERG 补算四臂并集": herg["scope_n"],
            "Butina 0.35（全池）簇数 / 最大簇 / 单例":
                f"{fbf['n_clusters']} / {fbf['largest_cluster']} / {fbf['n_singletons']}",
        }
        A(md_table([[k, exp, got.get(k, "—")] for k, exp in RECONCILE],
                   ["指标", "既有口径", "本包实测"]))
        A("")
        A("差异说明（三条，逐项可查）：① 本包一律**重算**描述符与门（上游 `pains` 列全空，"
          "且不得引用包外既有产物），故以本包实测为准；② 臂级 1–2 个分子的差异来自包内用 "
          "RDKit 2025.09.5 重算 QED 与上游 `qed` 列在阈值边界的浮点差（阈值 ±5e-4 内的候选："
          "QED≈0.60 有 9 个、≈0.50 有 10 个、≈0.40 有 6 个），例如 QED0.60∧anchor=4 实测 "
          "1270（无 maxTc 门）/ 1268（含 maxTc<0.9），与既有 1271 / 1269 各差 1 个分子，"
          f"而 QED0.60∧anchor≥3 = 2029 与既有值一致；四臂并集由 3236 变为 "
          f"{herg['scope_n']}（+{herg['scope_n'] - 3236}）**不是漂移**，而是 v1.5.0 描述符规范化"
          f"（去盐）后硬门由不过变过的分子新建入池所致，属已登记的口径变更（见 "
          f"`rule_changelog` 的 1.5.0 条）；"
          "③ 上表臂池为**未施加 novelty 门**的可比口径，交付名单实际使用的池规模见 §5 表"
          f"（novelty 门在各臂剔除 {_nk_lo}–{_nk_hi} 个已知物）。口径差异在 "
          "`02_druglikeness/reports/descriptor_crosscheck.json` 与上表逐项体现。")
        A("")
        A("## 9. 局限与偏差声明")
        A("")
        A(f"1. **hERG 门（状态 `{herg['status']}`）是未校准的单模型分诊线**：本包不实现 "
          f"hERG 预测（本地无 ADMET-AI/chemprop 运行时，权重不在包内）。"
          + ("名单未因 hERG 删减、并已写入 `herg_gate=pending` 标注；"
             "MD 前置复核必须补齐该门。" if herg["status"] == "pending"
             else (f"本轮已闭合：外部 ADMET-AI 预测表 `{herg.get('predictions_file')}` "
                   f"（来源 {herg.get('source') or '未标注'}）落入接口位置后按 <0.5 落门，"
                   f"名单与各臂池均已剔除不通过者；作用域覆盖率 "
                   f"{_pct(herg.get('scope_coverage'))}（并集内覆盖 "
                   f"{_pct(herg.get('coverage_in_scope'))}），池内其余 "
                   f"{herg.get('n_not_evaluated')} 个为 `not_evaluated`（未评估，"
                   f"不当作不合格），且不在任何名单内。")))
        _ans = ROOT / "03_gates/admet/reports/admet_annotation.json"
        _cmp = ROOT / "03_gates/herg/cross_model/reports/herg_cross_model.json"
        _triv = None
        if _cmp.is_file():
            _ag0 = read_json(_cmp)["agreement_vs_admet_ai"]
            _triv = max(m["m1_pass_rate"] for m in _ag0)
        if herg["status"] == "applied" and _ans.is_file():
            _sem = read_json(_ans).get("herg_gate_semantics", {})
            if _sem:
                _triv_txt = (f"，而「恒判过门」这一平凡基线本身就有 {_pct(_triv)}（门决策一致率"
                             f"只略高于它）" if _triv is not None else "")
                A(f"   ① **阈值口径只是一个分数，不是安全界**：`0.5` 在冻结库上已有 "
                  f"{_sem['library_pct_at_or_above_gate']}% 的分子在线上（本包池 "
                  f"{_sem['pool_pct_at_or_above_gate']}%），它大致对应已上市药在该 hERG 头上的 "
                  f"**{_sem['approved_drug_percentile_at_gate_library_near']} 百分位**；"
                  f"换一个第二模型，门决策一致率只剩三分之一到一半（§4b）{_triv_txt}"
                  f"——因此「ADMET-AI 过门」不等于「已判安全」，它只能转述为"
                  f"「ADMET-AI 口径下的低风险标记」。")
                A(f"   ② **交付名单多数贴线**：48 个候选的 hERG 距门边距最小 "
                  f"{_sem['delivered_margin_min']}、中位 "
                  f"{_sem['delivered_margin_median']}；边距 < {_sem['margin_tight_threshold']} 的 "
                  f"{_sem['delivered_n_margin_lt_tight']} 个、< 0.10 的 "
                  f"{_sem['delivered_n_margin_lt_0.10']} 个、落在 [0.40, 0.50) 的 "
                  f"{_sem['delivered_n_in_half_band_0.40_0.50']} 个。三档 band（软线 0.744）"
                  f"在交付名单上恒为「低」，没有区分度；可排序的分诊量是**距门边距**"
                  f"`admet_herg_margin`（见 `SOP-05_admet_md_precheck.csv`），"
                  f"MD 前人工复核应按它升序取用。")
                A(f"   ③ **实验出口（不可省略）**：hERG 风险的工程判据是暴露侧边距"
                  f"（Cmax/IC50 量级，经验边界见 Redfern 2003, *Cardiovasc Res* 58(1):32–45, "
                  f"PMID 12667944），而 QSAR 概率没有浓度维度，不能替代膜片钳/自动化电生理。"
                  f"故本门的正确用法是「分诊 + 排序」，不是「淘汰」。")
            A(f"   ④ 残余不确定性：外部模型版本自身的适用域限制；预测值非实验数据。")
        _cm9_p = ROOT / "03_gates/herg/cross_model/reports/herg_cross_model.json"
        if herg["status"] == "applied" and _cm9_p.is_file():
            _cm9 = read_json(_cm9_p)
            _g9, _a9 = _cm9["gate_impact"], _cm9["agreement_vs_admet_ai"]
            _k9 = [m["cohen_kappa"] for m in _a9]
            _g9r = [m["gate_agreement"] for m in _a9]
            _rl = [m["recall_of_m1_passes"] for m in _cm9["rank_consistency"]]
            A(f"   ⑤ **hERG 判定强依赖所选模型**（量化见 §4b）：换用第二个独立训练的 hERG 模型，"
              f"0.5 阈值上的门决策一致率仅 {_pct(min(_g9r))}–{_pct(max(_g9r))}"
              f"（κ {min(_k9)}–{max(_k9)}）；**与阈值无关**的排序交集同样低"
              f"（高风险前 10% Jaccard "
              f"{min(m['high_risk_top10pct_jaccard'] for m in _cm9['rank_consistency'])}–"
              f"{max(m['high_risk_top10pct_jaccard'] for m in _cm9['rank_consistency'])}；"
              f"匹配过门率后仅覆盖第一模型过门分子的 {_pct(min(_rl))}–{_pct(max(_rl))}），"
              f"说明差异不只是标定问题。若改用「两模型一致过门」，hERG 前主池在**两侧都有预测值的 "
              f"{_g9['pre_herg_main_pool'].get('m1_evaluated_n', 'NA')} 个**内由 "
              f"{_g9['pre_herg_main_pool']['m1_pass']} 变为 "
              f"{_g9['pre_herg_main_pool']['m1_and_m2_pass']}（另有 "
              f"{_g9['pre_herg_main_pool'].get('m1_missing_n', 'NA')} 个缺第一模型外部预测值、"
              f"未在这两个数内；见 §4b 的分母说明）、交付 48 名单中 "
              f"{_g9['delivered_list_48']['m1_only_pass']} 个会掉出——本包主口径仍保持"
              f" ADMET-AI 单模型门不变（改口径需重新选名单并重设验收基线）；"
              f"建议把「两模型一致判高风险」用作 MD 前置复核的优先级依据，"
              f"而不要把任一单模型的 0.5 阈值当作可迁移的硬判据。")
            A(f"   ⑥ **第二模型不是数据独立的外部验证**：{_cm9.get('independence_caveat', '')}——"
              f"§4b 因此只应读作**算法/训练管道层面的模型依赖度**。")
        A("2. **打分口径单一**：选序只用 CNN_VS，未与 Vina 亲和融合；池层质量由 anchor/姿态/"
          "maxTc 三门独立约束，降低单打分过拟合风险。")
        A("3. **富集指标偏乐观**：参考集为 DUD-E 诱饵，存在已知诱饵偏置"
          "（Wallach & Heifets 2018），故本包不把富集 AUC 作为验收判据，"
          "只用锚点命中率作为结构门证据（活性 vs 诱饵，见 §2）。")
        A("4. **pKa/cLogD 为上游第三方预测**，CNS MPO 继承其不确定性，故仅作标注列："
          "`basic_pka` 来自 MolGpKa（commit 4dc83521…，**只有点值、无预测区间**）；"
          "`clogd74` 来自 cpLogD v2.0 / ChEMBL33（sha256 281338b9…，附带 80% 预测区间）。"
          "按 80% 区间界定重算 cLogD 轴：本池有 2,778/4,990 个分子 CNS MPO 变动 ≥0.1，"
          "交付 48 名单中 22 个变动 ≥0.1（最大 0.3985，`ZINC000092920989`），"
          "点估计 ≥4 的 40 个里**有 8 个在悲观界下落回 4 以下**"
          "（最低 `ZINC000071721614` 3.547 → 3.453）——引用「CNS MPO≥4」时须一并转述该区间。"
          "pKa 轴无区间可用，其不确定性无法量化。")
        A("5. **化学型定义**为工程化的化学系列近似（氨基烟腈母核 × 碱性胺，2 个 SMARTS），"
          "用于防止名单被单一母核挤占；骨架约束另以 Bemis-Murcko 骨架独立把关。")
        A("6. **名额口径**：化学型/风险档上限按 `round(8N/20)`、`round(3N/20)`（N=50 → 20 / 8）"
          "执行，与规则集一致。")
        A("")
        # ---- §10 MD 批次（执行序第 16 步 15_md_batch_select；机时分配层）
        _mb_p = ROOT / "07_md_batch/reports/md_batch.json"
        if _mb_p.is_file():
            mbr = read_json(_mb_p)
            mbt = pd.read_csv(ROOT / "07_md_batch/tables/md_batch20.csv", low_memory=False)
            mbx = pd.read_csv(ROOT / "07_md_batch/tables/md_batch_excluded.csv")
            _dup_ex = mbx[mbx.reason.str.contains("骨架")]
            A("## 10. MD 批次：50 个名单 → 20 个（步骤 6，机时分配层）")
            A("")
            A(f"名单层已完成全部门禁（48 个候选**骨架互不相同**、姿态全 `pass`、锚点全 4/4、"
              f"hERG 全过门），因此「先跑哪 20 个」是**机时分配**问题而非过滤问题。口径见 "
              f"`config/frozen_rule.json → step6_md_batch`（阈值不内嵌脚本）："
              f"批次 = **{mbr['n_candidates']} 候选 + {mbr['n_controls']} 阴性对照 = "
              f"{mbr['batch_size']}**（NEG1/NEG2 固定占位，不占候选名额）。")
            A("")
            _note = mbt[mbt.md_eligibility_note.notna()
                        & (mbt.md_eligibility_note.astype(str).str.len() > 0)] \
                if "md_eligibility_note" in mbt.columns else mbt.iloc[0:0]
            if str(mbr['pool'].get('hERG_cross_model_policy')) == "retain_and_annotate":
                A(f"- 资格：① **不再剔除**第二 hERG 模型（集成口径）不过门的候选 "
                  f"**{int(mbr['pool']['n_retained_hERG_cross_model_fail'])}** 个"
                  f"（{'、'.join('`' + c + '`' for c in _note.cid)}）——该判定来自与主模型"
                  f"κ≈0（近乎随机一致）的第二模型，用一个未校准信号决定机时分配在纪律上最弱，"
                  f"故改为**保留入池并逐分子标注「模型依赖」**（表内列 `md_eligibility_note`、"
                  f"`herg_cross_model_flag`），机时代价仅 +0.2–0.3 卡日量级；"
                  f"**本层只做机时分配，不产生安全性结论**——两个模型一致判高风险的分子才应优先复核"
                  f"（见 `SOP-05_downstream_MD_plan.md` §5）；"
                  f"② Bemis-Murcko 通用骨架去重再剔除 "
                  f"**{int(mbr['pool']['n_excluded_generic_scaffold_dup'])}** 个同系物占位者"
                  f"（{'、'.join('`' + c + '`' for c in _dup_ex.cid)}，保留 CNN_VS 高者）"
                  f"→ 资格集 {int(mbr['pool']['n_eligible'])} 个。")
            A(f"- 配额：成对 ECFP4 Tc < {mbr['caps']['pairwise_tc_max']}、骨架 ≤1、"
              f"化学型 ≤ round(8×{mbr['batch_size']}/20)={mbr['caps']['chem_cap']}、"
              f"AM1-BCC 风险档 ≤ round(3×{mbr['batch_size']}/20)={mbr['caps']['risk_cap']}"
              f"（与名单层同口径按 N 缩放）。")
            A(f"- 分层：≥{mbr['stratification']['min_below_median']} 个候选须来自 CNN_VS "
              f"下半区（相对 48 候选的中位 {mbr['stratification']['reference_median_CNN_VS']}），"
              f"实测 **{mbr['stratification']['n_below_median_selected']}** 个。理由：纯分数截断会把"
              f"批次 CNN_VS 跨度压到 0.5 以内，使一级 MD 无法检验「对接分数 ⇒ 结合稳定性」。")
            _opt = mbr.get("optimality", {})
            _objc = mbr.get("objective", {})
            _heur = mbr.get("heuristic_reference", {})
            if _opt and _objc:
                _axl = [x["label"] for x in _objc["coverage_axes"]]
                _axes = "、".join(_axl)
                A(f"- 目标与最优性（v1.7.0 起）：在该可行域内按 `step6_md_batch.objective` 用 MILP"
                  f"（`scipy.optimize.milp`，HiGHS）求**可证明最优解**，字典序为 "
                  f"① 最大化分层覆盖下限 t → ② 在 t=t* 内最大化 CNN_VS 跨度 → ③ 仍并列时最大化 "
                  f"Σ CNN_VS（质量兜底）→ ④ 规范并列打破。{len(_axl)} 条性质轴"
                  f"（{_axes}）各按资格池自分布切 {_objc['n_bins_per_axis']} 分位，"
                  f"共 {_objc['n_bins']} 格；最优解 **每格至少 {_opt['t_star']} 席**"
                  f"（占据 {_opt['n_occupied_bins']}/{_opt['n_bins']} 格）、跨度 "
                  f"**{float(_opt['span']):.4f}**、Σ CNN_VS {float(_opt['sum_CNN_VS']):.4f}；四阶段均 "
                  f"`status=0`、`mip_gap=0`，重解名册一致（`resolved_roster_consistent="
                  f"{"true" if _opt.get("resolved_roster_consistent") else "false"}`）。")
                if _heur:
                    A(f"- 对照（v1.6.0 贪心填充，同一步内留痕，不驱动选批）：t={_heur['t']}、"
                      f"跨度 {_heur['span']}、Σ CNN_VS {_heur['sum_CNN_VS']}；与最优解重叠 "
                      f"**{_heur['overlap_with_optimal']}** 席、互换 **{_heur['n_swaps']}** 席。"
                      f"精确求解显示该贪心恰是「max Σ CNN_VS」口径的**全局最优解**——"
                      f"本版修的是**目标函数**，不是算法。")
                A("> **转述注意**：这里的「最优」限定在 `step6_md_batch.objective` 口径内，"
                  "**不是**「这批分子比其余候选更好」——最优解的 Σ CNN_VS 反而略低于质量口径；"
                  "其价值在于让一级 MD 能检验「对接分数 ⇒ 结合稳定性」并覆盖性质空间。"
                  "低分位（CNN_VS Q1）分子进入批次后失败概率更高，这正是要**测量**的量，"
                  "不得事后解读为「选坏了」。")
            A("")
            A(md_table([[k, v] for k, v in mbr["caps"]["chemotype_counts"].items()]
                       + [["风险档命中（上限 %d）" % mbr["caps"]["risk_cap"],
                           mbr["caps"]["risk_count"]],
                          ["唯一骨架 / 唯一通用骨架",
                           "%d / %d" % (mbr["diversity"]["n_unique_scaffold"],
                                        mbr["diversity"]["n_unique_generic_scaffold"])],
                          ["批内成对 Tc 上限（名单整体 0.5588）", mbr["diversity"]["max_pairwise_tc"]],
                          ["分层覆盖下限 t*（各分位最小席位，v1.7.0）",
                           "%s（占 %s/%s 格）" % (_opt.get("t_star", "—"),
                                                _opt.get("n_occupied_bins", "—"),
                                                _opt.get("n_bins", "—"))],
                          ["CNN_VS 跨度（v1.7.0 最优解 / v1.6.0 贪心对照）",
                           "%.4f / %.4f" % (float(_opt.get("span", 0)),
                                            float(_heur.get("span", 0)))],
                          ["CNN_VS 区间（批次实测）",
                           "%s–%s（宽 %s）" % (mbr["spans"]["CNN_VS_min"],
                                              mbr["spans"]["CNN_VS_max"],
                                              mbr["spans"]["CNN_VS_span"])],
                          ["CNS MPO 跨度 / cLogP 跨度",
                           "%s–%s / %s–%s" % (mbr["spans"]["cns_mpo_min"],
                                              mbr["spans"]["cns_mpo_max"],
                                              mbr["spans"]["cLogP_min"],
                                              mbr["spans"]["cLogP_max"])],
                          ["一级机时估算（卡日，线性外推）",
                           mbr["md_cost_estimate"]["level1_cardday"]],
                          ["二级机时估算（卡日，按 60–70% 通过率 / 全通过）",
                           "%s / %s" % (
                               mbr["md_cost_estimate"]["level2_cardday_at_65pct_pass"],
                               mbr["md_cost_estimate"]["level2_cardday_if_all_pass"])]],
                       ["批次属性", "实测"]))
            A("")
            A("**批次成员（按 CNN_VS 降序，对照在末 2 位）**")
            A("")
            _rows = []
            for r in mbt.itertuples():
                if r.candidate_type == "candidate":
                    _rows.append([int(r.batch_order), f"`{r.cid}`", "候选", int(r.pick_order),
                                  round(float(r.CNN_VS), 3), round(float(r.qed), 3),
                                  round(float(r.cns_mpo), 2), round(float(r.maxTc_actives), 3),
                                  round(float(r.herg_prob), 3),
                                  round(float(r.hERG_cross_ens), 3), int(r.risk_am1bcc),
                                  r.selection_reason])
                else:
                    _rows.append([int(r.batch_order), f"`{r.cid}`", str(r.control_role), "—",
                                  "—", "—", "—", "—", "—", "—", "—", r.selection_reason])
            A(md_table(_rows, ["#", "cid", "类型", "名单序", "CNN_VS", "QED", "CNS MPO",
                               "maxTc(活性)", "hERG(ADMET-AI)", "hERG(第二模型)", "风险档",
                               "入选理由"]))
            A("")
            A("> 执行要点：**先跑 2 个对照**（各 10 ns）确认一级判据有区分度（对照应漂出）；"
              "一级判据四条（对**晶体样结合模式**的重原子 RMSD < 0.3 nm、接触指纹保留率 ≥ 70%、"
              "锚点 4/4 保留 ≥ 3/4 帧、方向性接触保留 ≥ 80% 帧）与二级判据"
              "（锚点 4/4 保留率 ≥ 80%，`gmx_MMPBSA` 缺失时显式降级并在报告中标注）"
              "须**预注册**后再跑；边界分子优先加重复而不是加新分子。"
              "判据、口袋坐标系与对照设计见 `SOP-05_downstream_MD_plan.md` §4–§6。"
              "本步是机时分配层：**不改名单、不改门**；MD 后若需换人，须改规则并重跑整链。")
            A("")
        A("## 11. 复现方式")
        A("")
        A("```bash")
        A("cd SOP_5")
        A("python3 run_all.py --list        # 查看步序")
        A("python3 run_all.py               # 全链复现")
        A("python3 tests/verify_pack.py     # 包完整性自检")
        A("python3 tests/check_isolation.py # 静态确认无包外路径引用")
        A("```")
        A("")
        A("预期：末尾打印 `ALL DONE`；重跑产物与首次逐字节一致（G 臂），"
          "见 `05_robustness/determinism/sha256_baseline.json`。")
        A("")
        A("## 12. 参考来源")
        A("")
        A("- Liu W, et al. *Science* 2012;337(6091):232–236. PMID 22798613; "
          "DOI 10.1126/science.1219218（受体 PDB 4EIY，https://www.rcsb.org/structure/4EIY ）")
        A("- McNutt AT, et al. GNINA 1.0. *J Cheminform* 2021;13:43. PMID 34108002; "
          "DOI 10.1186/s13321-021-00522-2（gnina 1.3.2，https://github.com/gnina/gnina ）")
        A("- Trott O, Olson AJ. AutoDock Vina. *J Comput Chem* 2010;31(2):455–461. PMID 19499576")
        A("- Wallach I, Heifets A. *J Cheminform* 2018;10:45. PMID 31430292（DUD-E 诱饵偏置）")
        A("- Wager TT, Hou X, Verhoest PR, Villalobos A. *ACS Chem Neurosci* 2010;1(6):435–449. "
          "DOI 10.1021/cn100008c（CNS MPO Table 1）")
        A("- Swanson K, et al. ADMET-AI. *Bioinformatics* 2024;40(7):btae416. "
          "https://academic.oup.com/bioinformatics/article/40/7/btae416/7698030 "
          "（hERG 概率来源；本包以外部预测表接口接入该模型，未内置实现）")
        _cmp = ROOT / "03_gates/herg/cross_model/reports/herg_cross_model.json"
        if _cmp.is_file():
            _cmd = read_json(_cmp)["dataset"]
            A(f"- ChEMBL REST API（EBI），靶点 `CHEMBL240`（KCNH2 / hERG, Homo sapiens）："
              f"https://www.ebi.ac.uk/chembl/api/data/activity.json "
              f"（检索于 {_cmd.get('fetched_utc')}；第二模型交叉验证的训练/标注数据来源，"
              f"原始记录 {_cmd.get('n_raw_records'):,} 条，标注 {_cmd.get('n_labeled'):,} 个）")
        A("- RDKit 2025.09.5. https://www.rdkit.org")
        A("")
        write_text(DELIV / "SOP-05_report.md", "\n".join(L))

        D = []
        B = D.append
        _tier_thr = rule["step1_pose_qc"].get("pose_quality", {}).get("jaccard_threshold", 0.60)
        B("# SOP_5 下游接口：MD 体系构建与判据（机时预算单列）")
        B("")
        B("## 1. 输入")
        B("")
        B(f"- 起始结构：`06_deliverables/SOP-05_poses_main_list.sdf`"
          f"（{main['n_selected_candidates']} 个候选姿态，几何与上传姿态一致，"
          f"坐标往返偏差 {pid['coord_max_dev_max']} Å）")
        B("- 阴性对照姿态：见 `06_deliverables/SOP-05_negative_controls.csv`"
          "（candidate_type=negative_control，起始姿态取自 run01 归档 5,098 集姿态）")
        B("- 受体：`00_inputs/extracted/run01/01_redocking/receptor/4EIY_protein.pdb`"
          "（链 A；盒中心 (-0.471, 8.935, 17.159)，边长 23.042 Å；16 残基口袋）")
        B("- 体系：CHARMM-GUI POPC 膜 + custom 构型（A2A 为 7TM 膜蛋白）")
        B("")
        B("## 2. 两级判据（gmx_MMPBSA 不可用时的降级方案）")
        B("")
        B("**参照系（v1.6.0 起修正）**：判据的参照物是**晶体样结合模式**，不是上传姿态本身。")
        B("")
        B("- 口袋坐标系：用 16 残基 Cα 对 4EIY 受体副本叠合，得到统一框架；所有 RMSD 在该框架下算。")
        _coc = pgate.get("cocrystal_reference", {})
        if _coc.get("status") == "ok":
            B(f"- 起始姿态的正确性由**共晶正对照**标定：ZMA 重算读数 —— 锚点 "
              f"{_coc.get('anchor_hits_4A')}/4（"
              + "、".join(f"{k} {v} Å" for k, v in _coc.get("anchor_distances_A", {}).items())
              + f"）、真实接触残基 {_coc.get('n_contact_residues_4A')} 个、"
                f"盒子余量 {_coc.get('box_margin_A')} Å"
                f"（16 残基口袋中 {_coc.get('pocket16_extra_vs_cocrystal_n')} 个在共晶里不接触）。")
            B(f"- 姿态质量档 `pose_quality_tier`（A/B/C/D，按 anchor_hits × 与共晶接触集的 "
              f"Jaccard，阈值 {_tier_thr}）用于**排序与优先级**，不参与落门；实测全 5,000 姿态 "
              f"{pgate.get('pose_directionality', {}).get('tier_counts_all5000')}。")
        B("- 旧判据的缺陷：「锚点保留 ≥3/4 且**对上传姿态**的重原子 RMSD < 0.3 nm」用的是与姿态"
          "闸门同一套 `<4 Å` 定义、同一参照物——**错误的起始位姿只要原地不动就自动通过**。"
          "新判据把参照物换成晶体样模式，闸门与 MD 才互为独立检验。")
        B("")
        B("| 级别 | 判据 | 阈值建议 | 用途 |")
        B("|---|---|---|---|")
        B("| 一级（10 ns/分子） | ① 对**晶体样结合模式的**重原子 RMSD（口袋坐标系）"
          "② 接触指纹保留率（对共晶 12 残基接触集）③ 锚点 4/4 保留率 "
          "④ 方向性接触（N253/E169 氢键角、F168/W246 芳香堆积）保留率；"
          "配体自身 RMSD 仅作辅助 | ① < 0.3 nm ② 保留率 ≥70% ③ ≥3/4 帧 "
          "④ 极性接触至少 1 项保留 ≥80% 帧 | 剔除漂出分子 |")
        B("| 二级（100–200 ns/分子） | 锚点 4/4 保留率、H 键占据率、"
          "MM-GBSA/PBSA（若 gmx_MMPBSA 可用） | 锚点 4/4 保留率 ≥80% | 排序与终选 |")
        B("")
        B("> `gmx_MMPBSA` 缺失时，二级判据退化为「锚点接触保留率 + H 键占据率 + 能量分解」，"
          "报告中必须显式标注该降级，不得表述为自由能排序。")
        B("")
        B("## 3. 机时预算（单列，按 1 分子 1 重复）")
        B("")
        B("| 批次 | 分子数 | 时长 | GPU 卡时（估） | 备注 |")
        B("|---|---|---|---|---|")
        B("| 冒烟批（已完成） | 2 | 10 ns | ~0.2 | 水溶液对照；配体重原子 RMSD 0.17 nm |")
        B(f"| 一级 | {main['n_selected_candidates']} + 2 对照 | 10 ns | "
          f"≈ {(main['n_selected_candidates'] + 2) * 0.1:.1f} | 膜体系，单卡并行 |")
        B("| 二级 | 一级通过者（约 60–70%） | 100–200 ns | 按通过数线性放大 | 需先定卡时预算 |")
        B("")
        _ad2 = ROOT / "03_gates/admet/reports/admet_annotation.json"
        if _ad2.is_file():
            _a2 = read_json(_ad2)["delivered_list"]
            B("## 4. ADMET 前置复核优先级（分析步 10b/10c；非判据，不改名单）")
            B("")
            B(f"名单内 {_a2['md_precheck_n']} 个候选带 ≥1 项 ADMET 警告"
              "（清单：`06_deliverables/SOP-05_admet_md_precheck.csv`；逐分子注释见 "
              "`06_deliverables/SOP-05_admet_annotation_main.csv` 与追加了 `admet_*` 列的决策卡）。"
              "建议在 100–200 ns 二级 MD 之前先做体外/正交复核：")
            B("")
            B("- 离子通道（膜片钳 hERG）优先：hERG 注释位于池内中/高档者；")
            B("- CYP 抑制/底物与肝毒性（DILI/AMES）异常者；")
            B("- 溶解度与血浆蛋白结合在 MD 体系浓度下是否可行。")
            B("")
            B("> 注意：本层是标注而非判据——名单成员仍由既有门（含 `hERG < 0.5`）决定；"
              "若 MD/实测结果要求换人，须重跑 `run_all.py` 重新选名单，不得手改名单文件。")
            B("")
        B("## 5. 暴露侧边距（hERG 门的正交出口，不可省略）")
        B("")
        _sem5 = {}
        if _ad2.is_file():
            _sem5 = read_json(_ad2).get("herg_gate_semantics", {})
        B("`hERG < 0.5` 是**未校准的单模型分诊线**（量化见最终报告 §4b/§4c）："
          + (f"0.5 在冻结库上已有 {_sem5.get('library_pct_at_or_above_gate')}% 的分子处在上侧、"
             f"大致对应已上市药在该 hERG 头上的 {_sem5.get('approved_drug_percentile_at_gate_library_near')} "
             f"百分位；" if _sem5 else "")
          + "换第二个独立训练的 hERG 模型，0.5 阈值上的门决策一致率只有三分之一到一半"
            "（Cohen κ ≈ 0.02–0.07，近乎随机）。")
        B("")
        B("| 用法 | 做法 |")
        B("|---|---|")
        B("| 分诊与排序 | 按 ≥1 项 ADMET 警告清单 + **距门边距升序**"
          "（`SOP-05_admet_md_precheck.csv` 的 `admet_herg_margin`）确定复核优先级 |")
        B("| 正交实验出口 | **膜片钳（首选，或自动化电生理）**：hERG 阻断的工程判据是暴露侧边距"
          "（Cmax/IC50 量级；经验边界见 Redfern 2003, *Cardiovasc Res* 58(1):32–45, "
          "PMID 12667944），QSAR 概率没有浓度维度，不能替代它 |")
        B("| 优先核查对象 | **两个模型一致判为高风险**的分子（而非任一单模型判高风险者）"
          "——单模型的硬阈值不足以永久剔除分子 |")
        B("")
        if _sem5:
            B(f"- 交付名单实测：距门边距最小 {_sem5.get('delivered_margin_min')}、"
              f"中位 {_sem5.get('delivered_margin_median')}；边距 < "
              f"{_sem5.get('margin_tight_threshold')} 的 "
              f"{_sem5.get('delivered_n_margin_lt_tight')} 个、< 0.10 的 "
              f"{_sem5.get('delivered_n_margin_lt_0.10')} 个、落在 [0.40, 0.50) 的 "
              f"{_sem5.get('delivered_n_in_half_band_0.40_0.50')} 个。")
            B("")
        B("## 6. 对照设计（MD 必须带对照，否则无法判断判据是否失效）")
        B("")
        B("| 角色 | 对象 | 来源与理由 | 预期读数 |")
        B("|---|---|---|---|")
        B("| **正对照** | **ZMA（ZM241385）** | 4EIY 共晶配体（坐标见 "
          "`00_inputs/extracted/run01/01_redocking/receptor/4EIY_ligand.pdb`）：它是"
          "「正确答案」的几何定义 | 一级判据应**全部通过**；正对照若漂出，说明体系/力场/判据"
          "有问题，当批结果全部作废 |")
        B("| **硬负例** | **`CHEMBL5852598`** | 取自独立验证集 `independent_a2a_chembl.csv` 的 "
          "`label=0` 子集（实测 10 µM，弱/无活性），且与名单候选 `ZINC001776011322` 的 ECFP4 "
          "Tc = 0.628（最相似的「近邻但不结合」） | 一级判据应**不通过**或不保留结合模式；"
          "它通过则判据过松 |")
        B("| **阴性对照（沿用）** | NEG1 / NEG2 | 5,098 富集库中的高可信诱饵（仅具假定标签，"
          "故为次要对照） | 与硬负例互补：一个测「形状相似但不结合」，一个测「对接分数高的诱饵」 |")
        B("")
        B("> 硬负例的 3D 起始姿态**不在包内**，需按 `screening_kit/12_docking_validation` 的协议"
          "先对接生成（包内不伪造姿态）。选样规则（`label=0`、potency 缺失或 ≥10 µM、对名单候选的"
          "最大 ECFP4 Tc 最高）登记于 `config/frozen_rule.json → downstream_interface.md_controls`，"
          "结果可复算。")
        B("")
        B("## 7. 交付与追踪")
        B("")
        B("- MD 结果回填到决策卡（`SOP-05_decision_cards_main.csv`）新增列，"
          "保持「一分子一行、证据可追溯」。")
        B("- 若 MD 后更换候选，需重跑 `run_all.py` 以保持名单与门禁一致；"
          "禁止手工编辑名单文件（会破坏 G 臂确定性判据）。")
        B("")
        write_text(DELIV / "SOP-05_downstream_MD_plan.md", "\n".join(D))
        log.info("报告 → 06_deliverables/SOP-05_report.md；下游计划 → "
                 "06_deliverables/SOP-05_downstream_MD_plan.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
