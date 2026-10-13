#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""包完整性自检（交付前机器核验）。

检查项：
  ① 原始输入 md5 与 frozen_rule 一致（5 项）
  ② 关键派生输入哈希（受体 PDB、标签库、ChEMBL 独立集）
  ③ 规则一致性：frozen_rule 规范哈希 == rules_sha256；脚本哈希与 FREEZE.json 一致
  ④ 步骤产物齐备且规模合理（文件存在 + 行数/条数下限）
  ⑤ 名单质量硬判据（锚点 100%、成对 Tc 上限、名单规模、阴性对照 2 个）
  ⑥ CNS MPO 实现回归（原文三例，|Δ| ≤ 0.06）
  ⑦ 每个步骤都有漏斗 JSON
  ⑧ 结论文件（报告、验收对账、MANIFEST）存在且非空
  ⑭ ADMET-AI 注释层与 hERG 阈值锚定（分析步 10b/10c，非判据）：冻结表哈希登记、
     池覆盖、注释列齐备、交叉核验容差、阈值锚定复现、名单成员与顺序未变
  ⑲ 步序拓扑：读取者必须晚于写出者（按 run_order 读写路径静态核验，例外须登记）
  ⑳ 重跑鲁棒性：整目录删除后可由 run_all 预建目录重建（含功能测试）
退出码：0 = 全部通过；7 = 存在未通过项。
"""
from __future__ import annotations

import json
import sys

import numpy as np
from pathlib import Path

import pandas as pd

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.cns_mpo import regression_check
from vendor.io_utils import canonical_hash, md5_file, read_json, sha256_file

FAIL = []


def chk(name: str, cond: bool, detail: str = "") -> bool:
    print(f"  [{'OK  ' if cond else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not cond:
        FAIL.append(name)
    return cond


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    print("== SOP_5 自检 ==")

    print("\n[1] 原始输入哈希")
    for rel, spec in rule["inputs"].items():
        p = ROOT / "00_inputs/raw" / ("archives/" + rel if rel.endswith(".zip")
                                      else "inputs/" + rel)
        chk(f"input {rel}", p.is_file() and md5_file(p) == spec["md5"],
            f"md5={md5_file(p)[:12] if p.is_file() else 'NA'}…")

    print("\n[2] 关键派生输入")
    derive = {
        "00_inputs/extracted/run01/01_redocking/receptor/4EIY_protein.pdb":
            rule["receptor"]["md5_expected"],
        "00_inputs/extracted/run01/02_enrichment_5098/library/quick_5100_library.csv": None,
        "00_inputs/extracted/screening_kit/07_independent_validation/independent_a2a_chembl.csv": None,
    }
    for rel, exp in derive.items():
        p = ROOT / rel
        ok = p.is_file() and (exp is None or md5_file(p) == exp)
        chk(f"derived {Path(rel).name}", ok, f"{p.stat().st_size / 1e3:.0f} kB" if p.is_file() else "缺失")

    print("\n[3] 规则与脚本一致性")
    body = {k: v for k, v in rule.items() if k != "rules_sha256"}
    chk("frozen_rule.rules_sha256 匹配", canonical_hash(body) == rule.get("rules_sha256"),
        str(rule.get("rules_sha256"))[:16] + "…")
    fz = ROOT / "00_inputs/manifests/FREEZE.json"
    if fz.is_file():
        fr = read_json(fz)
        bad = [k for k, v in fr["scripts"].items()
               if not (ROOT / k).is_file() or sha256_file(ROOT / k) != v]
        chk("脚本哈希与 FREEZE.json 一致", len(bad) == 0,
            (f"变化/缺失 {len(bad)}：{bad[:3]}（FREEZE.json 由 01_verify_inputs 每轮重写；"
             f"若刚改过脚本又用 --from/--steps 局部跑，请先跑一次含 01_verify_inputs 的链以刷新登记）")
            if bad else f"{len(fr['scripts'])} 个脚本")
    else:
        chk("FREEZE.json 存在", False, "需先运行 01_verify_inputs")

    print("\n[4] 步骤产物齐备")
    products = {
        "01_pose_qc/tables/pose_pool.csv": 4900,
        "01_pose_qc/outputs/SOP-05_poses_normalized.sdf": 1,
        "01_pose_qc/tables/pose_gate.csv": 4900,
        "01_pose_qc/tables/labeled5098_contact_fp.csv": 1000,
        "02_druglikeness/tables/descriptors.csv": 4900,
        "02_druglikeness/tables/fingerprints_ecfp4.npz": None,
        "03_gates/tables/maxTc_actives.csv": 4900,
        "03_gates/tables/candidate_universe.csv": 4900,
        "03_gates/herg/SOP-05_herg_input_union.csv": 100,
        "04_selection/tables/list_main_QED0.50_anchor4.csv": 50,
        "05_robustness/tables/robA_sigma_sweep.csv": 3,
        # 2026-09-24（换数据 run05）：第二模型预测表按新表的四臂并集重建，规模随数据变；
        # 这里只守住“非空可读”下限，覆盖一致性由下方 _np == _n_union 断言精确保证。
        "03_gates/herg/cross_model/tables/herg_cross_model_predictions.csv": 100,
        "07_md_batch/tables/md_batch20.csv": 20,
    }
    for rel, min_rows in products.items():
        p = ROOT / rel
        if not p.is_file():
            chk(f"product {rel}", False, "缺失")
            continue
        if min_rows is None or p.suffix == ".npz":
            chk(f"product {rel}", p.stat().st_size > 0, f"{p.stat().st_size / 1e3:.0f} kB")
        else:
            n = sum(1 for _ in open(p, encoding="utf-8", errors="replace")) - 1
            chk(f"product {rel}", n >= min_rows, f"{n} 行（下限 {min_rows}）")

    print("\n[5] 名单质量硬判据")
    lp = ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv"
    if lp.is_file():
        df = pd.read_csv(lp, low_memory=False)
        cand = df[df.candidate_type == "candidate"]
        ctrl = df[df.candidate_type == "negative_control"]
        chk("主名单 50 行", len(df) == rule["step4_selection"]["list_size"], f"{len(df)} 行")
        chk("候选 48 + 对照 2", len(cand) == 48 and len(ctrl) == 2,
            f"{len(cand)} + {len(ctrl)}")
        chk("候选锚点全部 4/4", bool((cand.anchor_hits_4A == 4).all()),
            f"min={int(cand.anchor_hits_4A.min())}")
        chk("候选姿态闸门全部 pass", bool((cand.pose_gate == "pass").all()))
        chk("候选 QED 全部 ≥ 0.50", bool((cand.qed >= 0.50).all()),
            f"min={cand.qed.min():.3f}")
        chk("候选 maxTc 全部 < 0.90", bool((cand.maxTc_actives < 0.90).all()),
            f"max={cand.maxTc_actives.max():.3f}")
        chk("对照 cid 与规则期望一致",
            set(ctrl.cid.astype(str)) == {rule["step4_selection"]["controls"]["negative_1"]["expected_cid"],
                                          rule["step4_selection"]["controls"]["negative_2"]["expected_cid"]},
            str(sorted(ctrl.cid.astype(str))))
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from vendor import chem_core as cc
        fps = {}
        des = pd.read_csv(ROOT / "02_druglikeness/tables/descriptors.csv", low_memory=False).set_index("cid")
        for c in cand.cid:
            fps[c] = cc.fp_ecfp4(Chem.MolFromSmiles(des.loc[c, "smiles"]))
        from vendor.picker import max_pairwise_tc
        tcm = max_pairwise_tc(cand, fps)
        chk("成对 Tc 上限 ≤ 0.65", tcm <= 0.65, f"max={tcm}")
        chk("一骨架一分子（主口径）", int(cand.scaffold.nunique()) == len(cand),
            f"{cand.scaffold.nunique()} 骨架 / {len(cand)} 分子")

    print("\n[6] CNS MPO 回归（原文 Table 3）")
    rc = regression_check(verbose=False)
    chk("|Δ| ≤ 0.06", rc["max_abs_diff"] <= 0.06, f"max|Δ|={rc['max_abs_diff']}")

    print("\n[7] 每步漏斗 JSON")
    funnels = ["02_druglikeness/reports/druglikeness_funnel.json",
               "03_gates/reports/gate_stack.json",
               "04_selection/reports/selection_summary.json",
               "05_robustness/reports/robustness.json",
               "01_pose_qc/reports/pose_identity_summary.json",
               "01_pose_qc/reports/pose_gate_summary.json"]
    for f in funnels:
        p = ROOT / f
        chk(f"funnel {Path(f).name}", p.is_file() and p.stat().st_size > 200)

    print("\n[8] 结论文件")
    for f in ["06_deliverables/SOP-05_report.md", "06_deliverables/INDEX.md",
              "06_deliverables/acceptance_report.json", "06_deliverables/MANIFEST.sha256",
              "06_deliverables/SOP-05_MD_shortlist_main.csv",
              "06_deliverables/SOP-05_poses_main_list.sdf",
              "06_deliverables/SOP-05_downstream_MD_plan.md"]:
        p = ROOT / f
        chk(f"deliverable {Path(f).name}", p.is_file() and p.stat().st_size > 100,
            f"{p.stat().st_size / 1e3:.0f} kB" if p.is_file() else "缺失")

    print("\n[9] SDF 结构合法性与几何一致性")
    from vendor import chem_core as cc
    up = {}
    _sdf_in = ROOT / "00_inputs/raw/inputs/top5000_best_pose.sdf"
    if not _sdf_in.is_file():
        # 精简骨架未含大体积输入（见 RESTORE_LIST.md）：登记缺失并继续，
        # 后续小件检查仍能给出读数，不以 traceback 中止
        chk("输入姿态 SDF 存在", False,
            f"缺失 {_sdf_in.relative_to(ROOT)}（骨架未含大体积输入，见 RESTORE_LIST.md）")
    else:
        for block, props, title in cc.iter_sdf_blocks(_sdf_in):
            up[props.get("cid", title)] = cc.molblock_data_lines(block)
    for rel, n_exp in [("01_pose_qc/outputs/SOP-05_poses_normalized.sdf", None),
                       ("06_deliverables/SOP-05_poses_main_list.sdf", None)]:
        p = ROOT / rel
        if not p.is_file():
            chk(f"SDF {Path(rel).name}", False, "缺失")
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        n_dollar = txt.count("$$$$")
        n_mend = sum(1 for l in txt.splitlines() if l.strip() == "M  END")
        chk(f"SDF {Path(rel).name} 每块有 M  END", n_dollar == n_mend,
            f"$$$$={n_dollar} / M_END={n_mend}")
        # 几何一致性：抽样比对数据行（标题/counts/原子/键/M END）
        blocks = [b for b in txt.split("$$$$") if b.strip()]
        step = max(len(blocks) // 50, 1)
        bad = 0
        for b in blocks[::step]:
            core = cc.molblock_data_lines(b.strip("\n"))
            cid = core[0].strip()
            if cid in up and core[1:] != up[cid][1:]:
                bad += 1
        chk(f"SDF {Path(rel).name} 几何与上传一致（抽样 {len(blocks[::step])} 条）", bad == 0,
            f"不一致 {bad} 条")

    print("\n[10] 外部接口输入（hERG 预测表）")
    hspec = rule.get("step3_gates", {}).get("herg", {})
    hp = ROOT / hspec.get("external_table", "03_gates/herg/herg_predictions.csv")
    if hp.is_file():
        import pandas as _pd
        _h = _pd.read_csv(hp, low_memory=False)
        cols = [c for c in ("hERG", "herg_prob") if c in _h.columns] or \
               [c for c in _h.columns if "herg" in c.lower()]
        chk("hERG 表含 cid + hERG 列", ("cid" in _h.columns) and bool(cols),
            f"列={list(_h.columns)[:6]}")
        if cols:
            v = _h[cols[0]]
            chk("hERG 取值在 [0,1]", bool(v.dropna().between(0, 1).all()),
                f"min={v.min():.4f} max={v.max():.4f} n={v.notna().sum()}")
        reg = hspec.get("external_sha256")
        if reg:
            chk("hERG 表哈希与规则登记一致", sha256_file(hp) == reg, f"{sha256_file(hp)[:12]}…")
        else:
            chk("hERG 表已在规则中登记哈希", False, "external_sha256 为空")
    else:
        chk("hERG 表存在（缺失 → 该门 pending，不算失败）", True, "pending 降级模式")

    print("\n[11] CNS MPO 输入与 80% 区间标注列")
    _card = ROOT / "06_deliverables/SOP-05_MD_shortlist_main.csv"
    if _card.is_file():
        import pandas as _pd2
        _c = _pd2.read_csv(_card, low_memory=False)
        _need = ["basic_pka", "clogd74", "clogd74_lower_80", "clogd74_upper_80",
                 "cns_mpo_lo80", "cns_mpo_hi80", "cns_mpo_pi80_span"]
        _miss = [c for c in _need if c not in _c.columns]
        chk("决策卡含 pKa/cLogD 输入与区间标注列", not _miss, f"缺={_miss}")
        if not _miss:
            _s = _c.dropna(subset=["cns_mpo_lo80", "cns_mpo_hi80"])
            chk("区间界包含点估计（lo80 ≤ mpo ≤ hi80，容差 1e-6）",
                bool(((_s.cns_mpo_lo80 <= _s.cns_mpo + 1e-6).all())
                     and ((_s.cns_mpo <= _s.cns_mpo_hi80 + 1e-6).all())),
                f"n={len(_s)}")
            chk("span == hi80 − lo80",
                bool(((_s.cns_mpo_pi80_span - (_s.cns_mpo_hi80 - _s.cns_mpo_lo80))
                      .abs() <= 1e-6).all()))
            _q = _c.dropna(subset=["clogd74", "clogd74_lower_80", "clogd74_upper_80"])
            chk("cLogD 区间包含点值（lower ≤ 点值 ≤ upper）",
                bool(((_q.clogd74_lower_80 <= _q.clogd74 + 1e-9)
                      & (_q.clogd74 <= _q.clogd74_upper_80 + 1e-9)).all()),
                f"n={len(_q)}（阴性对照不在上游描述符表内，无此三列 → 已排除）")
    # master 取数不得出现重复列后缀（否则下游「列存在才写入」会静默丢列）
    try:
        import sys as _sys
        _sys.path.insert(0, str(ROOT))
        from vendor.master import build_master as _bm
        _mm = _bm(ROOT)
        _dup = [c for c in _mm.columns if c.endswith("_x") or c.endswith("_y")]
        chk("主表无重复列（pKa/cLogD 以单一来源合并）", not _dup, f"重复列={_dup[:5]}")
    except Exception as _e:                                   # noqa: BLE001
        chk("主表可构建", False, str(_e)[:80])

    print("\n[12] 包级清单自洽性")
    _pm = ROOT / "PACK_MANIFEST.sha256"
    if _pm.is_file():
        _lines = [l for l in _pm.read_text(encoding="utf-8").splitlines() if l.strip()]
        chk("PACK_MANIFEST 不含自身条目（排除自指）",
            not any(l.split(None, 1)[-1].strip() == "PACK_MANIFEST.sha256" for l in _lines),
            f"条目数={len(_lines)}")
        _idx = ROOT / "INDEX_PACK.md"
        if _idx.is_file():
            _n = [int(x) for x in __import__("re").findall(r"清单文件数[^0-9]*(\d+)",
                                                            _idx.read_text(encoding="utf-8"))]
            chk("INDEX_PACK 记录的文件数与清单条目一致",
                bool(_n) and _n[0] == len(_lines), f"index={_n[:1]} manifest={len(_lines)}")
    else:
        chk("PACK_MANIFEST 存在", False, "缺失")

    print("\n[13] hERG 第二模型交叉验证（分析步 08b）")
    _cm_dir = ROOT / "03_gates/herg/cross_model"
    _cm_js = _cm_dir / "reports/herg_cross_model.json"
    if not _cm_js.is_file():
        chk("cross-model 结论 JSON 存在", False, "缺失")
    else:
        _cm = json.loads(_cm_js.read_text(encoding="utf-8"))
        _need = ("dataset", "cv", "agreement_vs_admet_ai", "rank_consistency",
                 "gate_impact", "independence_caveat")
        _miss = [k for k in _need if k not in _cm]
        chk("cross-model JSON 字段齐备", not _miss, f"缺 {_miss}" if _miss else "6/6")
        _n_union = _cm["gate_impact"]["union_recomputed"]["pool_n"]

        # 原始数据哈希三方对账：规则 ↔ fetch_meta ↔ 实际文件
        _rule_cm = (read_json(ROOT / "config/frozen_rule.json")
                    .get("step3_gates", {}).get("herg", {}).get("cross_model", {}))
        _raw = ROOT / str(_rule_cm.get("raw_table", ""))
        _meta_p = _cm_dir / "raw/fetch_meta.json"
        chk("cross-model 原始数据表存在", _raw.is_file(), _rule_cm.get("raw_table", ""))
        if _raw.is_file():
            _h_md5 = md5_file(_raw)
            _h_sha = sha256_file(_raw)
            chk("原始数据 md5 与规则一致", _h_md5 == _rule_cm.get("raw_md5"),
                f"{_h_md5}")
            chk("原始数据 sha256 与规则一致", _h_sha == _rule_cm.get("raw_sha256"),
                f"{_h_sha[:16]}…")
            if _meta_p.is_file():
                _meta = json.loads(_meta_p.read_text(encoding="utf-8"))
                chk("fetch_meta 哈希与规则一致",
                    _meta.get("md5") == _rule_cm.get("raw_md5")
                    and _meta.get("sha256") == _rule_cm.get("raw_sha256"),
                    f"n_rows={_meta.get('n_rows')}")
                chk("fetch_meta 记录数与 JSON 一致",
                    _meta.get("n_rows") == _cm["dataset"].get("n_raw_records"),
                    f"meta={_meta.get('n_rows')} json={_cm['dataset'].get('n_raw_records')}")

        # 预测表规模与并集一致
        _pred = _cm_dir / "tables/herg_cross_model_predictions.csv"
        if _pred.is_file():
            _np = sum(1 for _ in open(_pred, encoding="utf-8", errors="replace")) - 1
            chk("第二模型预测表覆盖整个并集", _np == _n_union,
                f"{_np} 行 vs 并集 {_n_union}")
        else:
            chk("第二模型预测表存在", False, "缺失")

        # 计数自洽：both_pass + m1_only = m1_pass；both_pass + m2_only = m2_pass；四格和 = pool_n
        _gi = _cm["gate_impact"]
        _bad = []
        _fields = {"pool_n", "m1_pass", "m2_pass", "m1_and_m2_pass",
                   "m1_only_pass", "m2_only_pass"}
        _arms = [k for k, v in _gi.items()
                 if isinstance(v, dict) and _fields <= set(v)]
        chk("gate_impact 含 6 个池（键名改动不影响本断言）", len(_arms) == 6,
            "、".join(_arms))
        for _arm in _arms:
            _g = _gi.get(_arm)
            if not _g:
                _bad.append(f"{_arm} 缺失")
                continue
            if _g["m1_and_m2_pass"] + _g["m1_only_pass"] != _g["m1_pass"]:
                _bad.append(f"{_arm} m1 计数")
            if _g["m1_and_m2_pass"] + _g["m2_only_pass"] != _g["m2_pass"]:
                _bad.append(f"{_arm} m2 计数")
            if (_g["m1_pass"] + _g["m2_only_pass"]) > _g["pool_n"]:
                _bad.append(f"{_arm} 超池")
        chk("gate_impact 计数自洽（6 个池）", not _bad, "；".join(_bad) if _bad else "6/6")

        # 一致性表口径：m1_pass_m2_fail / m1_fail_m2_pass 与 both_* 互补，且行数 = 4 个第二模型
        _ag = _cm["agreement_vs_admet_ai"]
        chk("一致性表含 4 个第二模型",
            len(_ag) == 4 and {r["model2"] for r in _ag} == {"rf", "lr", "knn", "ens"},
            "/".join(r["model2"] for r in _ag))
        _sums = [r["both_pass"] + r["both_fail"] + r["m1_pass_m2_fail"] + r["m1_fail_m2_pass"]
                 for r in _ag]
        chk("一致性混淆四格和 = 并集规模",
            all(s == _n_union for s in _sums), f"{_sums}")
        _pmj = _cm["gate_impact"].get("pre_herg_main_pool", {})
        chk("前主池登记外部表覆盖基数（covered_n / uncovered_n）",
            isinstance(_pmj.get("covered_n"), int) and isinstance(_pmj.get("uncovered_n"), int)
            and _pmj["covered_n"] + _pmj["uncovered_n"] == _pmj.get("pool_n"),
            f"覆盖 {_pmj.get('covered_n')} + 未覆盖 {_pmj.get('uncovered_n')} "
            f"= {_pmj.get('pool_n')}")
        # 第一模型外部预测覆盖：与其缺值约定必须显式登记，且 m1 计数不得跨基数
        _m1bad = []
        for _k in _arms:
            _g2 = _cm["gate_impact"][_k]
            if _g2.get("m1_evaluated_n", 0) + _g2.get("m1_missing_n", 0) != _g2.get("pool_n"):
                _m1bad.append(f"{_k} covered_n 拆分")
            if _g2.get("m1_pass", 0) > _g2.get("m1_evaluated_n", 0):
                _m1bad.append(f"{_k} m1_pass>{_k}.m1_evaluated_n")
        chk("m1 覆盖拆分自洽（evaluated + missing = pool_n）且 m1_pass 不跨基数",
            not _m1bad, "；".join(_m1bad) if _m1bad else
            f"前主池 evaluated={_pmj.get('m1_evaluated_n')} + missing={_pmj.get('m1_missing_n')}")
        chk("登记第一模型缺值约定（不得静默当不过门）",
            "m1_missing_convention" in _cm["gate_impact"],
            _cm["gate_impact"].get("m1_missing_convention", "")[:48])
        chk("gate_impact 含键名口径说明（不再内嵌规模数字）",
            "key_note" in _cm["gate_impact"], _cm["gate_impact"].get("key_note", "")[:40])
        chk("模糊带/标注规模未变（阻断+非阻断 = 标注）",
            _cm["dataset"]["n_blocker"] + _cm["dataset"]["n_nonblocker"]
            == _cm["dataset"]["n_labeled"],
            f"{_cm['dataset']['n_labeled']}")

        # 交付层副本
        _ri = read_json(ROOT / "REPRO_INFO.json")
        chk("REPRO_INFO 额外输入计数与清单一致",
            _ri.get("counts", {}).get("declared_extra_inputs") == len(_ri.get("extra_inputs", [])),
            f"{_ri.get('counts', {}).get('declared_extra_inputs')} vs {len(_ri.get('extra_inputs', []))}")
        _ro = read_json(ROOT / "config/run_order.json")
        chk("REPRO_INFO 步数与步序文件一致",
            _ri.get("steps") == len(_ro.get("steps", [])),
            f"{_ri.get('steps')} vs {len(_ro.get('steps', []))}")

        for _f in ("SOP-05_herg_cross_model_consistency.csv",
                   "SOP-05_herg_cross_model_rank_consistency.csv",
                   "SOP-05_herg_cross_model_list_flips.csv",
                   "figures/SOP-05_fig_herg_cross_model.png"):
            chk(f"交付副本 {_f}", (ROOT / "06_deliverables" / _f).is_file())

    print("\n[14] ADMET-AI 注释层与 hERG 阈值锚定（分析步 10b/10c，非判据）")
    _al = rule.get("step3_gates", {}).get("admet_layer", {})
    _ht_rule = rule.get("step3_gates", {}).get("herg_threshold_analysis", {})
    _raw_p = ROOT / str(_al.get("raw_table", ""))
    if not _raw_p.is_file():
        chk("ADMET 冻结表存在", False, "缺失 → 10b 无法运行")
    else:
        _sha = sha256_file(_raw_p)
        chk("冻结表 sha256 与规则登记一致", _sha == _al.get("raw_sha256"), f"{_sha[:16]}…")
        chk("冻结表 md5 与规则登记一致", md5_file(_raw_p) == _al.get("raw_md5"))
        with open(_raw_p, encoding="utf-8", errors="replace") as _f:
            _nrow = sum(1 for _ in _f) - 1
        chk("冻结表行数与规则登记一致", _nrow == _al.get("n_rows_expected"),
            f"{_nrow} vs {_al.get('n_rows_expected')}")
        _cols = pd.read_csv(_raw_p, nrows=0).columns
        _na = len([c for c in _cols if c.startswith("ADMET_")])
        chk("冻结表 ADMET 列数与规则登记一致",
            _na == _al.get("n_admet_columns_expected"),
            f"{_na} vs {_al.get('n_admet_columns_expected')}")

    _cards_p = ROOT / "03_gates/admet/tables/admet_annotation_cards.csv"
    _pool_p = ROOT / "03_gates/admet/tables/admet_annotation_pool.csv"
    _ad_js = ROOT / "03_gates/admet/reports/admet_annotation.json"
    _need_ann = list(_al.get("annotation_columns", []))
    if not (_cards_p.is_file() and _pool_p.is_file() and _ad_js.is_file()):
        chk("注释层产物齐备（cards/pool/json）", False, "缺失 → 需运行 10b")
    else:
        _cards = pd.read_csv(_cards_p, low_memory=False)
        _pool = pd.read_csv(_pool_p, low_memory=False)
        _miss = [c for c in _need_ann if c not in _cards.columns]
        chk("注释列齐备（Cards 表）", not _miss, f"缺 {_miss}" if _miss else f"{len(_need_ann)} 列")
        chk("池覆盖数与规则期望一致",
            len(_pool) == _al.get("expected_pool_coverage"),
            f"{len(_pool)} vs {_al.get('expected_pool_coverage')}")
        _ad = json.loads(_ad_js.read_text(encoding="utf-8"))
        chk("注释层 JSON 覆盖计数自洽",
            _ad["coverage"]["covered"] == len(_pool)
            and _ad["coverage"]["covered"] + len(_ad["coverage"]["not_covered"])
            == _ad["coverage"]["pool_n"],
            f"covered={_ad['coverage']['covered']} pool={_ad['coverage']['pool_n']}")
        chk("注释层声明不干预门/池/名单",
            bool(_ad.get("non_interference", {}).get("statement"))
            and _ad["non_interference"]["gate_threshold"]
            == rule["step3_gates"]["herg"]["threshold"],
            f"gate={_ad.get('non_interference', {}).get('gate_threshold')}")

        # 描述符独立复算交叉核验：数值超容差 = 0；判定一致率（BRENK 差异已登记）
        _xc_p = ROOT / "03_gates/admet/reports/admet_descriptor_crosscheck.json"
        if not _xc_p.is_file():
            chk("交叉核验 JSON 存在", False, "缺失")
        else:
            _xc = json.loads(_xc_p.read_text(encoding="utf-8"))
            _over = {k: v["n_over_tol"] for k, v in _xc["numeric"].items()
                     if v.get("n_over_tol")}
            chk("描述符数值交叉核验无超容差项", not _over and _xc["all_numeric_within_tol"],
                f"超容差 {_over}" if _over else f"{len(_xc['numeric'])} 个量，容差 {_xc['tolerance']}")
            _ra = _xc["rule_agreements"]
            # 口径：本包已按最大片段去盐、ADMET-AI 未去盐 → 严格一致率只在**单片段**子集上判定；
            # 多片段子集单列（numeric_multi_fragment / n_multi_fragment），不判失败。
            _r5s, _pns = _ra["ro5"].get("agreement_pct_single_fragment"), \
                _ra["pains"].get("agreement_pct_single_fragment")
            chk("Ro5/PAINS 判定在单片段子集上 100% 一致",
                _r5s == 100.0 and _pns == 100.0,
                f"单片段：Ro5={_r5s} PAINS={_pns}；全池（含未去盐的多片段子集）"
                f"Ro5={_ra['ro5']['agreement_pct']} PAINS={_ra['pains']['agreement_pct']}")
            # 2026-09-24（换数据 run05）：本批 4,979 个 ML 分子的上游 SMILES 全为单片段
            # → n_multi_fragment = 0，“多片段差异组”不适用。判据改为：口径必须始终登记；
            # 差异组存在时才要求它单列且带 numeric_multi_fragment。
            _nmf = int(_xc.get("n_multi_fragment", 0))
            _mf_reg = "multi_fragment" in _xc.get("expected_differences", {})
            chk("多片段子集单列且登记为预期差异",
                _mf_reg and (_nmf == 0 or bool(_xc.get("numeric_multi_fragment"))),
                (f"多片段 0 个 → 本批上游 SMILES 全为单片段，差异组不适用；口径已登记={_mf_reg}"
                 if _nmf == 0 else
                 f"多片段 {_nmf} 个，已登记={_mf_reg}、单列块={'numeric_multi_fragment' in _xc}"))
            # BRENK：允许两种合法状态——① 有差异且已登记（原判据）；② 两边命中集完全相同
            # （一致率 100%，无差异可登记）；其余（有差异未登记）一律不判通过。
            _bk = _ra["brenk"]
            _bk_no_diff = int(_bk.get("pack_only", 0)) == 0 and int(_bk.get("admet_only", 0)) == 0
            chk("BRENK 差异已登记为已知不判失败",
                "BRENK" in _xc.get("expected_differences", {})
                and (_bk["agreement_pct"] < 100.0 or _bk_no_diff),
                f"BRENK={_bk['agreement_pct']}%"
                + ("（本批两边命中集相同，无差异）" if _bk_no_diff else ""))
            _mr = _xc.get("main_fragment_redundancy", {}) if "main_fragment_redundancy" in _xc else \
                {}
            if _mr:
                chk("去盐后身份冗余已登记（盐/游离碱同物种对）",
                    _mr.get("n_groups", 0) > 0 and _mr.get("n_groups_involving_multi_fragment", 0)
                    == _mr.get("n_groups", 0),
                    f"{_mr.get('n_groups')} 组 / {_mr.get('n_molecules')} 个分子")

    # hERG 阈值锚定：与外部流程数字复现（软线 / 库级取舍线）
    _ht_js = ROOT / "03_gates/admet/reports/herg_threshold_analysis.json"
    if not _ht_js.is_file():
        chk("阈值锚定 JSON 存在", False, "缺失 → 需运行 10c")
    else:
        _ht = json.loads(_ht_js.read_text(encoding="utf-8"))
        _tol = float(_ht_rule.get("repro_check", {}).get("tol", 0.005))
        _rc = _ht["repro_check"]
        chk("软线复现（0.744 ± %s）" % _tol,
            abs(_ht["decision"]["soft_line"] - _ht_rule.get("soft_line")) <= _tol
            and _rc["soft_line_ok"],
            f"observed={_rc['soft_line_observed']} expected={_rc['soft_line_expected']}")
        chk("库级取舍线复现（0.80 ± %s）" % _tol,
            abs(_ht["A3_tradeoff"][f"library_head@K={_ht['params']['primary_top_k']}"]
                ["feasible_t"] - _ht_rule.get("library_head_tradeoff_expected")) <= _tol
            and _rc["library_head_tradeoff_ok"],
            f"observed={_rc['library_head_tradeoff_observed']} "
            f"expected={_rc['library_head_tradeoff_expected']}")
        chk("hERG 门阈值未变（帧冻值 0.5）",
            _ht["decision"]["frozen_gate_threshold"] == rule["step3_gates"]["herg"]["threshold"],
            f"{_ht['decision']['frozen_gate_threshold']}")
        _gap = _ht.get("gate_scope_gap") or {}
        chk("门作用域缺口证据齐备且计数自洽",
            bool(_gap) and _gap["in_scope"] + _gap["out_of_scope"] == _gap["pool_n"]
            and (_gap["out_of_scope_frozen_below_gate"]
                 + _gap["out_of_scope_frozen_at_or_above_gate"] == _gap["out_of_scope"])
            and _gap["overlap_max_abs_diff"] <= 1e-5,
            f"并集外 {_gap.get('out_of_scope')}（其中 {_gap.get('out_of_scope_frozen_below_gate')} "
            f"个 hERG 低于门）、交集 max|Δ|={_gap.get('overlap_max_abs_diff')}")
        _rcs = ROOT / "03_gates/admet/tables/herg_retention_curve.csv"
        if _rcs.is_file():
            _rc_df = pd.read_csv(_rcs)
            _heads = [c for c in _rc_df.columns if c.startswith("keep_")]
            _mono = all(bool((_rc_df[c].diff().dropna() >= -1e-9).all()) for c in _heads)
            _grid = _rc_df.t.round(2).diff().dropna().round(6)
            chk("保留曲线覆盖 9 个头部、0.01 步长、保留率随阈值单调不减",
                len(_heads) >= 8 and _mono and len(_rc_df) > 50 and bool((_grid == 0.01).all()),
                f"{len(_heads)} 头部 / {len(_rc_df)} 档（t {_rc_df.t.min()}–{_rc_df.t.max()}）")
            # 独立复算：取曲线最大档，用注释表重算库/池/交付名单三个头部的保留率
            _t = float(_rc_df.t.max())
            _obs = {}
            _pl = pd.read_csv(_pool_p, low_memory=False)
            _obs["keep_pool_4990_pct"] = float((_pl.ADMET_hERG < _t).mean()) * 100
            _wild = pd.read_csv(ROOT / "06_deliverables/SOP-05_admet_annotation_main.csv",
                                low_memory=False)
            _cov = _wild[_wild.admet_coverage == "covered"]
            _obs["keep_delivered_list_48_pct"] = float((_cov.ADMET_hERG < _t).mean()) * 100
            _row = _rc_df.loc[_rc_df.t == _rc_df.t.max()].iloc[0]
            _bad = {k: (round(float(_row[k]), 3), round(v, 3))
                    for k, v in _obs.items() if abs(float(_row[k]) - v) > 0.01}
            chk(f"保留曲线最大档（t={_t}）独立复算一致", not _bad,
                f"不一致 {_bad}" if _bad else "池 / 交付名单两个头部一致（±0.01 pp）")

    # 交付层：注释列已追加，且成员/顺序与原名单一致（非判据不删分子）
    _src_map = {
        "SOP-05_MD_shortlist_main.csv": "04_selection/tables/list_main_QED0.50_anchor4.csv",
        "SOP-05_decision_cards_main.csv": "04_selection/tables/decision_cards_main_QED0.50_anchor4.csv",
        "SOP-05_negative_controls.csv": None,
    }
    for _dst, _src in _src_map.items():
        _dp = ROOT / "06_deliverables" / _dst
        if not _dp.is_file():
            chk(f"交付件 {_dst}", False, "缺失")
            continue
        if _src is None:
            continue
        _d = pd.read_csv(_dp, dtype=str, keep_default_na=False)
        _s = pd.read_csv(ROOT / _src, dtype=str, keep_default_na=False)
        chk(f"{_dst} 成员与顺序与原名单一致",
            _d.cid.tolist() == _s.cid.tolist(), f"{len(_d)} 行")
        _newc = [c for c in _d.columns if c not in _s.columns]
        chk(f"{_dst} 追加列均为注释列（{len(_newc)} 列）",
            set(_newc) <= set(_need_ann) and len(_newc) == len(_need_ann),
            f"非注释新增列={sorted(set(_newc) - set(_need_ann))[:3]}")
    _wide_p = ROOT / "06_deliverables/SOP-05_admet_annotation_main.csv"
    if _wide_p.is_file():
        _wd = pd.read_csv(_wide_p, dtype=str, keep_default_na=False)
        _ml = pd.read_csv(ROOT / "06_deliverables/SOP-05_MD_shortlist_main.csv",
                          dtype=str, keep_default_na=False)
        chk("注释全表行序与主名单一致",
            _wd.cid.tolist() == _ml.cid.tolist(), f"{len(_wd)} 行")
        _nc = _wd.loc[_wd.admet_coverage == "not_covered", "cid"].tolist()
        chk("阴性对照标注 not_covered（无预测值不写空值）",
            len(_nc) == 2 and all(str(c) in set(_wd.cid) for c in _nc), str(_nc))
    else:
        chk("注释全表存在", False, "缺失")

    # 报告小节与交付件互指有效
    _rep = ROOT / "06_deliverables/SOP-05_report.md"
    if _rep.is_file():
        _rt = _rep.read_text(encoding="utf-8")
        chk("报告含 4c 非判据小节", "### 4c." in _rt and "非判据" in _rt)
        _cited = ["06_deliverables/SOP-05_admet_md_precheck.csv",
                  "06_deliverables/SOP-05_herg_retention_curve.csv",
                  "06_deliverables/SOP-05_HERG_THRESHOLD_SOP.md"]
        _bad = [c for c in _cited if c in _rt and not (ROOT / c).is_file()]
        chk("报告引用的注释层交付件均存在", not _bad, f"缺 {_bad}" if _bad else f"{len(_cited)} 项")
    _acc_p = ROOT / "06_deliverables/acceptance_report.json"
    if _acc_p.is_file():
        _acc = json.loads(_acc_p.read_text(encoding="utf-8"))
        _acx = _acc.get("context", {}).get("admet_annotation_layer", {})
        chk("验收对账含注释层非判据上下文",
            bool(_acx) and _acx.get("non_gating") is True
            and _acx.get("herg_gate_unchanged") == rule["step3_gates"]["herg"]["threshold"],
            f"coverage={_acx.get('pool_coverage')} gate={_acx.get('herg_gate_unchanged')}")
        chk("验收对账记录「名单成员未变 + 交叉核验容差 + 阈值复现」",
            _acx.get("list_membership_unchanged") is True
            and _acx.get("annotation_merge_membership_asserted") is True
            and _acx.get("descriptor_crosscheck_all_within_tol") is True
            and _acx.get("prediction_layer_recomputed_in_pack") is False
            and bool(_acx.get("herg_threshold_repro_check", {}).get("soft_line_ok")),
            f"membership={_acx.get('list_membership_unchanged')} "
            f"tol={_acx.get('descriptor_crosscheck_tolerance')} "
            f"soft_ok={_acx.get('herg_threshold_repro_check', {}).get('soft_line_ok')}")
    else:
        chk("验收对账存在", False, "缺失")

    print("\n[15] MD 批次选择（执行序第 16 步 15_md_batch_select）")
    _mb_rule = rule.get("step6_md_batch", {})
    _mb_dir = ROOT / "07_md_batch"
    _mb_csv = _mb_dir / "tables/md_batch20.csv"
    _mb_js = _mb_dir / "reports/md_batch.json"
    if not _mb_rule:
        chk("规则登记 step6_md_batch", False, "缺失")
    else:
        chk("规则登记 step6_md_batch", True,
            f"batch_size={_mb_rule.get('batch_size')} "
            f"chem_cap={_mb_rule.get('caps', {}).get('chem_cap_value')}")
    if not (_mb_csv.is_file() and _mb_js.is_file()):
        chk("MD 批次产物存在", False, "缺失 → 需运行 15_md_batch_select")
    else:
        _mb = pd.read_csv(_mb_csv, low_memory=False)
        _mbj = json.loads(_mb_js.read_text(encoding="utf-8"))
        _b_cand = _mb[_mb.candidate_type == "candidate"]
        _b_ctrl = _mb[_mb.candidate_type == "negative_control"]
        chk("批次规模与规则一致",
            len(_mb) == _mb_rule.get("batch_size")
            and len(_b_cand) == _mb_rule.get("n_candidates")
            and len(_b_ctrl) == _mb_rule.get("list", {}).get("n_negative_controls"),
            f"{len(_mb)} = {len(_b_cand)} 候选 + {len(_b_ctrl)} 对照")
        _exp_ctrl = {v["expected_cid"] for k, v in
                     rule["step4_selection"]["controls"].items() if k.startswith("negative")}
        chk("阴性对照沿用主名单 NEG1/NEG2",
            set(_b_ctrl.cid) == _exp_ctrl, f"{sorted(set(_b_ctrl.cid))}")
        # 成员合法性：候选必须来自主名单，且不得是第二模型 hERG 不过门者
        _main = pd.read_csv(ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv",
                            low_memory=False)
        _main_ids = set(_main.loc[_main.candidate_type == "candidate", "cid"])
        chk("批次候选均来自主口径名单",
            set(_b_cand.cid).issubset(_main_ids),
            f"越界 {sorted(set(_b_cand.cid) - _main_ids)[:3]}")
        _flips_p = ROOT / "03_gates/herg/cross_model/tables/herg_cross_model_list_flips.csv"
        if _flips_p.is_file():
            _flips = set(pd.read_csv(_flips_p).cid)
            _excl = rule["step6_md_batch"]["eligibility"].get(
                "exclude_hERG_cross_model_fail", True)
            if _excl:
                chk("第二 hERG 模型不过门者未入批",
                    set(_b_cand.cid).isdisjoint(_flips),
                    f"混入 {sorted(set(_b_cand.cid) & _flips)}")
            else:
                chk("第二 hERG 模型不过门者保留且逐分子带标注（v1.6.0 策略）",
                    all(bool(str(_mb.loc[_mb.cid == c, "md_eligibility_note"].iloc[0]).strip())
                        for c in (set(_b_cand.cid) & _flips))
                    if (set(_b_cand.cid) & _flips) else True,
                    f"策略=retain_and_annotate；清单 {sorted(_flips)}；"
                    f"入批者 {sorted(set(_b_cand.cid) & _flips)}")
        else:
            chk("第二模型偏离清单存在", False, "缺失")
        # 配额
        _chem_cap = int(_mb_rule["caps"]["chem_cap_value"])
        _cnt = _b_cand.chemotype.value_counts().to_dict()
        chk(f"化学型 ≤{_chem_cap}", max(_cnt.values()) <= _chem_cap, f"{_cnt}")
        _risk_cap = int(_mb_rule["caps"]["risk_cap_value"])
        _riskn = int(_b_cand.risk_am1bcc.sum())
        chk(f"AM1-BCC 风险档 ≤{_risk_cap}", _riskn <= _risk_cap, f"{_riskn}")
        chk("批次内 Bemis-Murcko / 通用骨架均 ≤1",
            _b_cand.scaffold.nunique() == len(_b_cand)
            and _b_cand.scaffold_generic.nunique() == len(_b_cand),
            f"{len(_b_cand)} 个分子 / {_b_cand.scaffold.nunique()} 骨架")
        # 分层
        _med = float(_main.loc[_main.candidate_type == "candidate", "CNN_VS"].median())
        _nlow = int((_b_cand.CNN_VS < _med).sum())
        chk("CNN_VS 下半区候选数 ≥ 规则下限",
            _nlow >= int(_mb_rule["stratification"]["min_below_median_candidates"]),
            f"{_nlow} vs {_mb_rule['stratification']['min_below_median_candidates']}"
            f"（中位 {round(_med, 4)}）")
        # 成对 Tc（从包内指纹重算，独立于脚本内部计算）
        _z = np.load(ROOT / "02_druglikeness/tables/fingerprints_ecfp4.npz", allow_pickle=True)
        _fps = {}
        from rdkit import DataStructs
        from rdkit.DataStructs.cDataStructs import ExplicitBitVect as _EBV
        for _c, _bits in zip(_z["cid"], _z["bits"]):
            _bv = _EBV(2048)
            _bv.SetBitsFromList(np.flatnonzero(_bits).tolist())
            _fps[str(_c)] = _bv
        _ids = _b_cand.cid.tolist()
        _tc = max(DataStructs.TanimotoSimilarity(_fps[_ids[i]], _fps[_ids[j]])
                  for i in range(len(_ids)) for j in range(i + 1, len(_ids)))
        chk(f"批内成对 ECFP4 Tc < {_mb_rule['caps']['pairwise_tc_max']}",
            _tc < float(_mb_rule["caps"]["pairwise_tc_max"]), f"max={round(_tc, 4)}")
        # JSON 与 CSV 同源
        chk("批次 JSON 与 CSV 计数一致",
            _mbj["n_candidates"] == len(_b_cand) and _mbj["batch_size"] == len(_mb),
            f"{_mbj['n_candidates']}/{_mbj['batch_size']}")
        _pol = _mbj["pool"].get("hERG_cross_model_policy")
        _exp_flip = _mb_rule["eligibility"]["expected_n_hERG_cross_model_fail"]
        chk("批次 JSON 记录 hERG 交叉模型判定与留痕",
            _mbj["pool"]["n_excluded_hERG_cross_model_fail"]
            + _mbj["pool"].get("n_retained_hERG_cross_model_fail", 0) == _exp_flip
            and (_mb_dir / "tables/md_batch_excluded.csv").is_file()
            and (_mb_dir / "MD_BATCH.md").is_file(),
            f"策略={_pol}；剔除 {_mbj['pool']['n_excluded_hERG_cross_model_fail']}、"
            f"保留并标注 {_mbj['pool'].get('n_retained_hERG_cross_model_fail')}、"
            f"骨架去重 {_mbj['pool']['n_excluded_generic_scaffold_dup']}")
        chk("hERG 交叉模型策略与规则一致",
            (_pol == "retain_and_annotate")
            == (_mb_rule["eligibility"].get("exclude_hERG_cross_model_fail") is False),
            f"规则 exclude={_mb_rule['eligibility'].get('exclude_hERG_cross_model_fail')}，"
            f"实测策略={_pol}")
        if _pol == "retain_and_annotate":
            _flips_p2 = ROOT / "03_gates/herg/cross_model/tables/herg_cross_model_list_flips.csv"
            if _flips_p2.is_file():
                _fl = set(pd.read_csv(_flips_p2).cid)
                _in_batch = set(_b_cand.cid) & _fl
                chk("跨模型翻转分子若入批则带标注",
                    all("标注" in str(_mb.loc[_mb.cid == c, "md_eligibility_note"].iloc[0])
                        for c in _in_batch) if _in_batch else True,
                    f"入批翻转分子 {sorted(_in_batch)}（未入批的分子仍保留在资格池）")
        chk("批次声明不干预名单/门",
            "不改名单" in _mbj.get("note", "") or "不改门" in _mbj.get("note", ""),
            _mbj.get("note", "")[:40])
        # 步序：新步必须在 13_assemble 之前，且 12_determinism 仍是末步
        _ord = json.loads((ROOT / "config/run_order.json").read_text(encoding="utf-8"))["steps"]
        _ids_order = [x["id"] for x in _ord]
        chk("步序：15_md_batch_select 在 13_assemble 之前且 12_determinism 末位",
            "15_md_batch_select" in _ids_order
            and _ids_order.index("15_md_batch_select") < _ids_order.index("13_assemble")
            and _ids_order[-1] == "12_determinism",
            f"共 {len(_ids_order)} 步，位置 "
            f"{_ids_order.index('15_md_batch_select') + 1}")
    for _f in ("SOP-05_MD_batch20.csv", "SOP-05_MD_batch_excluded.csv",
               "SOP-05_MD_batch20.md", "SOP-05_md_batch.json",
               "figures/SOP-05_fig_md_batch.png"):
        chk(f"交付副本 {_f}", (ROOT / "06_deliverables" / _f).is_file())

    print("\n[16] hERG 作用域披露与门语义（评审 P0-1 / P0-2）")
    _hs = ROOT / "03_gates/reports/herg_status.json"
    if _hs.is_file():
        _h = json.loads(_hs.read_text(encoding="utf-8"))
        _need = ("scope_coverage", "n_not_evaluated", "n_pool_total", "scope_n")
        _miss = [k for k in _need if k not in _h]
        chk("herg_status 含作用域披露字段", not _miss, f"缺 {_miss}" if _miss else "4/4")
        if not _miss:
            _gsn = _h.get("gate_scope_n", _h["scope_n"])
            chk("scope_coverage == gate_scope_n / n_pool_total",
                abs(_h["scope_coverage"] - _gsn / _h["n_pool_total"]) <= 1e-4,
                f"{_h['scope_coverage']} vs {_gsn}/{_h['n_pool_total']}"
                f"（口径 {_h.get('scope_mode')}）")
            chk("n_not_evaluated == 池 − 作用域",
                _h["n_not_evaluated"] == _h["n_pool_total"] - _gsn,
                f"{_h['n_not_evaluated']}")
            if _h.get("scope_mode") == "pool":
                chk("全池口径下作用域覆盖率 = 100%（无未评估分子）",
                    _h["scope_coverage"] == 1.0 and _h["n_not_evaluated"] == 0,
                    f"{_h['scope_coverage']} / {_h['n_not_evaluated']}")
                chk("全池口径下落门依据与外部表逐值对账 ≤1e-6",
                    _h.get("external_vs_gate_max_abs_diff", 1) <= 1e-6,
                    f"交集 {_h.get('external_vs_gate_overlap_n')}，"
                    f"max|Δ|={_h.get('external_vs_gate_max_abs_diff')}")
        # 报告必须与 JSON 同步（防止文本与数据漂移）
        _rp = ROOT / "06_deliverables/SOP-05_report.md"
        if _rp.is_file():
            _txt = _rp.read_text(encoding="utf-8", errors="replace")
            _pct = f"{_h.get('scope_coverage', 0) * 100:.2f}%"
            chk("报告的作用域覆盖率与 herg_status 一致", _pct in _txt,
                f"期望出现 {_pct}")
            chk("报告不再把并集外分子断言为「保守不过门」",
                "按规则保守" not in _txt and "not_evaluated" in _txt,
                "已改为 not_evaluated（未评估）；引用旧说法仅作为被否定的表述保留")
    else:
        chk("herg_status.json 存在", False, "缺失")

    _wsem = ROOT / "03_gates/admet/reports/admet_annotation.json"
    if _wsem.is_file():
        _s = json.loads(_wsem.read_text(encoding="utf-8"))
        chk("注释层含 hERG 门语义块", "herg_gate_semantics" in _s,
            "键=" + ",".join(sorted(_s.get("herg_gate_semantics", {}))[:4]) + "…")
        _ann = ROOT / "06_deliverables/SOP-05_admet_annotation_main.csv"
        if _ann.is_file():
            _a = pd.read_csv(_ann, dtype=str, keep_default_na=False)
            chk("交付注释表含距门边距列",
                "admet_herg_margin" in _a.columns, f"列数={len(_a.columns)}")
            if "admet_herg_margin" in _a.columns:
                _hc = ("admet_hERG" if "admet_hERG" in _a.columns else "ADMET_hERG")
                _num = pd.to_numeric(_a.admet_herg_margin, errors="coerce")
                _hv = pd.to_numeric(_a[_hc], errors="coerce")
                _ok = ((_num - (0.5 - _hv)).abs() <= 1e-6) | _hv.isna()
                chk("距门边距 == 0.5 − admet_hERG", bool(_ok.all()),
                    f"n={int(_num.notna().sum())}")
    else:
        chk("admet_annotation.json 存在", False, "缺失")

    _pg = ROOT / "01_pose_qc/reports/pose_gate_summary.json"
    if _pg.is_file():
        _p = json.loads(_pg.read_text(encoding="utf-8")).get("cocrystal_reference", {})
        chk("姿态步含共晶正对照参照读数", _p.get("status") == "ok",
            f"ligand={_p.get('ligand')}")
        if _p.get("status") == "ok":
            chk("正对照锚点命中 4/4（口径自证）", _p.get("anchor_hits_4A") == 4,
                f"{_p.get('anchor_hits_4A')}/4；距离 "
                f"{_p.get('anchor_distances_A')}")
            chk("正对照给出盒子余量", isinstance(_p.get("box_margin_A"), (int, float)),
                f"{_p.get('box_margin_A')} Å（半边长 {_p.get('box_half_edge_A')}）")
        _gt = ROOT / "01_pose_qc/tables/pose_gate.csv"
        if _gt.is_file():
            _g = pd.read_csv(_gt, low_memory=False)
            chk("姿态表含共晶接触集 Jaccard 列",
                "contact_jaccard_cocrystal" in _g.columns,
                f"列数={len(_g.columns)}")
    else:
        chk("pose_gate_summary.json 存在", False, "缺失")

    print("\n[17] 确定性判据的写入语义（评审 P1）")
    _dc = ROOT / "05_robustness/reports/determinism_compare.json"
    if _dc.is_file():
        _d = json.loads(_dc.read_text(encoding="utf-8"))
        _need = ("mode", "baseline_sha256", "determinism_verified")
        _miss = [k for k in _need if k not in _d]
        chk("判定 JSON 含 mode/baseline_sha256/determinism_verified", not _miss,
            f"缺 {_miss}" if _miss else f"mode={_d['mode']}")
        if not _miss:
            chk("写基线的一轮不得记为确定性通过",
                (_d["mode"] != "baseline_written")
                or (_d["determinism_verified"] is False and _d["identical"] is None),
                f"mode={_d['mode']}, verified={_d['determinism_verified']}, "
                f"identical={_d['identical']}")
            chk("mode=compared 时 determinism_verified=true",
                _d["mode"] != "compared" or _d["determinism_verified"] is True)
    _dets = ROOT / "05_robustness/scripts/12_determinism.py"
    if _dets.is_file():
        _src = _dets.read_text(encoding="utf-8")
        chk("确定性哈希范围含 07_md_batch", "\"07_md_batch\"" in _src,
            "此前该目录不在比对范围内")
    # 交付层在 G 臂之前执行：验收表不得把「上一轮链末读数」记为本轮的确定性通过。
    _accp = ROOT / "06_deliverables/acceptance_report.json"
    if _accp.is_file():
        _a = json.loads(_accp.read_text(encoding="utf-8"))
        _detrows = [c for c in _a["checks"] if "确定性" in c["criterion"]]
        chk("验收表含两条确定性判据", len(_detrows) == 2, f"n={len(_detrows)}")
        chk("验收表的确定性行在交付层恒为 pending（不得自证）",
            all(c["pass"] is None and str(c.get("status", "")).startswith("pending")
                for c in _detrows),
            "、".join(str(c.get("status")) for c in _detrows))
        _dv = _a.get("determinism_verdict", {})
        chk("验收表登记链末判据来源", "determinism_compare.json" in str(_dv.get("source", ""))
            and "12_determinism" in str(_dv.get("judged_by", "")),
            f"source={_dv.get('source')}")
        chk("验收表登记上一轮链末读数而非本轮结论",
            isinstance(_dv.get("last_recorded"), dict)
            and "last_recorded" in _dv,
            f"last_recorded={_dv.get('last_recorded')}")
    _repp = ROOT / "06_deliverables/SOP-05_report.md"
    if _repp.is_file():
        _rt = _repp.read_text(encoding="utf-8")
        _sec7 = _rt.split("## 7. 步骤 5")[-1].split("## 8.")[0]
        chk("报告 §7 区分上一轮链末读数与本轮判据",
            ("上一轮链末" in _sec7) and ("本步之后" in _sec7), "已加时点声明")
        chk("报告 §7 指回链末判据文件与机器断言",
            ("determinism_compare.json" in _sec7) and ("[17]" in _sec7), "文件 + 自检组号")
        _sec8 = _rt.split("## 8. 验收对账")[-1].split("## 9.")[0]
        chk("报告 §8 的确定性行标 pending 而非通过",
            "pending（由链末 G 臂判定" in _sec8, "两行均 pending")

    print("\n[18] 报告 §10 与 MD 批次口径一致（评审 P0-2b / P0-3）")
    _mbj = ROOT / "07_md_batch/reports/md_batch.json"
    _rep = ROOT / "06_deliverables/SOP-05_report.md"
    if _mbj.is_file() and _rep.is_file():
        _mbp = json.loads(_mbj.read_text(encoding="utf-8")).get("pool", {})
        _pol = str(_mbp.get("hERG_cross_model_policy", ""))
        _rt = _rep.read_text(encoding="utf-8")
        _sec10 = _rt.split("## 10. MD 批次")[-1].split("## 11.")[0] if "## 10. MD 批次" in _rt else ""
        chk("md_batch.json 登记 hERG 跨模型策略", _pol != "", f"policy={_pol}")
        if _pol == "retain_and_annotate":
            chk("报告 §10 不再断言「剔除第二 hERG 模型」",
                "剔除第二 hERG 模型" not in _sec10,
                "已改述为保留+标注")
            chk("报告 §10 转述保留与标注口径",
                ("保留入池" in _sec10) and ("模型依赖" in _sec10),
                f"n_retained={_mbp.get('n_retained_hERG_cross_model_fail')}")
            chk("报告 §10 列举被保留并标注的分子",
                all(c in _sec10 for c in
                    pd.read_csv(ROOT / "07_md_batch/tables/md_batch20.csv",
                                low_memory=False)
                    .query("md_eligibility_note == md_eligibility_note")
                    .loc[lambda d: d.md_eligibility_note.astype(str).str.len() > 0, "cid"]),
                "2 个")
        _exq = pd.read_csv(ROOT / "07_md_batch/tables/md_batch_excluded.csv")
        _dupc = _exq.loc[_exq.reason.str.contains("骨架"), "cid"].tolist()
        chk("报告 §10 列举同系物骨架去重名单",
            all(c in _sec10 for c in _dupc) if _dupc else True, f"{len(_dupc)} 个")
        chk("报告 §10 转述新一级判据（晶体样参照系）",
            ("晶体样" in _sec10) and ("≥ 70%" in _sec10),
            "RMSD/指纹/锚点/方向性四条")
        chk("报告 §10 指回 MD 计划判据与对照",
            "SOP-05_downstream_MD_plan.md" in _sec10, "§4–§6")
    else:
        chk("报告 §10 一致性核验前置文件", False, "缺 md_batch.json 或报告")

    print("\n[19] 步序拓扑：读取者必须晚于写出者（评审 D1）")
    import ast as _ast
    _ord = json.loads((ROOT / "config/run_order.json").read_text(encoding="utf-8"))["steps"]
    _idx = {st["id"]: i for i, st in enumerate(_ord)}
    _writer = {}
    for i, st in enumerate(_ord):
        for out in (st.get("outputs") or []):
            _writer[out] = (i, st["id"])
    # 允许的例外（必须写明理由，否则本组断言失败）：
    _allow = {
        ("13_assemble", "05_robustness/reports/determinism_compare.json"):
            "交付层位于 G 臂之前，只能登记上一轮读数（已在验收表标 pending）",
        ("14_report", "05_robustness/reports/determinism_compare.json"):
            "同上：报告 §7 明示「读数取自上一轮链末」",
    }
    _viol, _seen = [], []
    for i, st in enumerate(_ord):
        src = (ROOT / st["script"]).read_text(encoding="utf-8")
        lits = {n.value for n in _ast.walk(_ast.parse(src))
                if isinstance(n, _ast.Constant) and isinstance(n.value, str)}
        for lit in lits:
            hit = _writer.get(lit)
            if hit is None:
                continue
            j, wid = hit
            if j <= i:
                continue                      # 正常方向
            _seen.append(f"{st['id']}→{wid}")
            if (st["id"], lit) not in _allow:
                _viol.append(f"{st['id']}(第{i + 1}步) 读 {wid}(第{j + 1}步) 的 {lit}")
    chk("无未登记的「读下游产物」边", not _viol,
        "；".join(_viol) if _viol else f"已知例外 {len(_seen)} 条均已登记")
    chk("已知例外均有理由", all(v for v in _allow.values()), f"{len(_allow)} 条")

    # D1 回归：08b 必须晚于选名单，且不得从交付层（06_deliverables）取输入
    _08b = ROOT / "03_gates/scripts/08b_herg_cross_model.py"
    if _08b.is_file() and "08b_herg_cross_model" in _idx:
        _src08b = _08b.read_text(encoding="utf-8")
        chk("08b 位于 10_select_lists 之后",
            _idx["08b_herg_cross_model"] > _idx["10_select_lists"],
            f"第 {_idx['08b_herg_cross_model'] + 1} 步 vs 名单第 {_idx['10_select_lists'] + 1} 步")
        _read_deliv = [l for l in _src08b.splitlines()
                       if "06_deliverables" in l and not l.strip().startswith("#")]
        chk("08b 不从交付层取输入（单一数据源）", not _read_deliv,
            "；".join(_read_deliv)[:80] if _read_deliv else "只读 04_selection")
        chk("08b 读上游名单产物",
            "04_selection/tables/list_main_QED0.50_anchor4.csv" in _src08b, "04_selection 主名单")
    chk("run_order 步骤序号连续且与数组顺序一致",
        [st.get("order") for st in _ord] == list(range(1, len(_ord) + 1)),
        f"{len(_ord)} 步")
    chk("每个分析步都声明 role（不得冒充门/名单变更）",
        all(st.get("role") for st in _ord if st["id"] in
            {"08b_herg_cross_model", "10b_admet_annotation", "10c_herg_threshold",
             "15_md_batch_select"}),
        "4/4")

    print("\n[20] 重跑鲁棒性：产物目录可整体删除后重建（收尾项）")
    import importlib.util as _ilu
    import shutil as _sh
    _rasrc = (ROOT / "run_all.py").read_text(encoding="utf-8")
    chk("run_all 定义 ensure_output_dirs（按 run_order.outputs 预建目录）",
        "def ensure_output_dirs" in _rasrc, "起链前预建父目录")
    chk("run_all 在起链前调用预建（而非只打印失败提示）",
        "ensure_output_dirs(steps)" in _rasrc, "在 _run 内调用")
    _no_out = [st["id"] for st in _ord if not (st.get("outputs") or [])]
    chk("run_order 每步都声明 outputs（预建覆盖可与步序对齐）", not _no_out,
        "；".join(_no_out) if _no_out else f"{len(_ord)}/{len(_ord)} 步")
    _miss = sorted({str((ROOT / o).parent.relative_to(ROOT))
                    for st in _ord for o in (st.get("outputs") or [])
                    if not (ROOT / o).parent.is_dir()})
    chk("当前包内全部声明产物的父目录存在", not _miss,
        "；".join(_miss) if _miss else f"{len({str((ROOT / o).parent) for st in _ord for o in (st.get('outputs') or [])})} 个目录")
    # 功能测试：删除整目录后，预建必须把它建回来（只建目录、不建文件）
    _spec = _ilu.spec_from_file_location("_sop5_run_all_probe", ROOT / "run_all.py")
    _mod = _ilu.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    _probe_rel = "logs/_dirtest_probe/deep/out.csv"
    _probe = ROOT / _probe_rel
    _sh.rmtree(ROOT / "logs/_dirtest_probe", ignore_errors=True)
    _created = _mod.ensure_output_dirs([{"id": "_probe", "outputs": [_probe_rel]}])
    chk("整目录缺失时预建可重建（功能测试）",
        _probe.parent.is_dir(), "；".join(_created) or "未建")
    chk("预建只建目录、不生成占位文件", not _probe.is_file(), "不会污染产物")
    _sh.rmtree(ROOT / "logs/_dirtest_probe", ignore_errors=True)
    chk("功能测试残留已清理", not (ROOT / "logs/_dirtest_probe").exists(), "临时目录已删")

    # 冻结输入层（raw/）不在 run_order.outputs 内，但被步骤读取；被误删必须能立刻发现
    _raw_layer = ["03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv",
                  "03_gates/herg/cross_model/raw/chembl_herg_kcnh2_activities.csv",
                  "03_gates/herg/cross_model/raw/fetch_meta.json"]
    _raw_missing = [r for r in _raw_layer if not (ROOT / r).is_file()]
    chk("冻结输入层 raw/ 齐备（不随产物目录清空而消失）", not _raw_missing,
        "；".join(_raw_missing) if _raw_missing else f"{len(_raw_layer)} 个文件")
    _08bsrc = (ROOT / "03_gates/scripts/08b_herg_cross_model.py").read_text(encoding="utf-8")
    chk("08b 对 raw/ 缺失给出可操作提示（不抛裸异常）",
        "缺少冻结输入" in _08bsrc, "SystemExit 带修复指引")

    print("\n[21] MD 批次最优性证书（v1.7.0 目标口径；评审「MD 候选应是最优子集」）")
    from rdkit import DataStructs as _QDS
    from rdkit.DataStructs.cDataStructs import ExplicitBitVect as _QEBV
    _q_js = ROOT / "07_md_batch/reports/md_batch.json"
    _q_csv = ROOT / "07_md_batch/tables/md_batch20.csv"
    _q_exq = ROOT / "07_md_batch/tables/md_batch_excluded.csv"
    _q_rule = rule.get("step6_md_batch", {})
    chk("规则登记 MD 目标函数块 objective",
        bool(_q_rule.get("objective")) and bool(_q_rule.get("expected_certificate")),
        f"轴={len(_q_rule.get('objective', {}).get('coverage_axes', []))} "
        f"K={_q_rule.get('objective', {}).get('n_bins_per_axis')} "
        f"阶段={len(_q_rule.get('objective', {}).get('lexicographic', []))}")
    if not (_q_js.is_file() and _q_csv.is_file() and _q_exq.is_file()):
        chk("最优性证书前置文件齐备", False,
            "缺 md_batch.json / md_batch20.csv / md_batch_excluded.csv")
    else:
        _QJ = json.loads(_q_js.read_text(encoding="utf-8"))
        _QOBJ, _QOPT = _QJ.get("objective", {}), _QJ.get("optimality", {})
        _QHR = _QJ.get("heuristic_reference", {})
        _QEXP = _q_rule.get("expected_certificate", {})
        _QAX = [(a["column"], a["label"]) for a in _QOBJ.get("coverage_axes", [])]
        _QK = int(_QOBJ.get("n_bins_per_axis", 0))
        chk("证书存在且声明为证书", _QOPT.get("certificate") is True,
            f"t*={_QOPT.get('t_star')} span={_QOPT.get('span')} Σ={_QOPT.get('sum_CNN_VS')}")
        _QLEX = _QOBJ.get("lexicographic", [])
        chk("字典序四阶段均 status=0 且 mip_gap=0",
            len(_QLEX) == len(_q_rule["objective"]["lexicographic"])
            and all(x["status"] == 0 and float(x["mip_gap"]) <= 1e-9 for x in _QLEX),
            "；".join(f"{x['name']}: status={x['status']} gap={x['mip_gap']}" for x in _QLEX))
        chk("重解名册一致（确定性）", _QOPT.get("resolved_roster_consistent") is True,
            "同一输入两次求解名册逐位相同")

        # ---- 独立重建资格池与分位：主名单候选 − 剔除表 → 资格池自分布分位（不复用步骤脚本中间量）
        _q_main = pd.read_csv(ROOT / "04_selection/tables/list_main_QED0.50_anchor4.csv",
                              low_memory=False)
        _q_cand = _q_main[_q_main.candidate_type == "candidate"]
        _q_excl_ids = set(pd.read_csv(_q_exq).cid)
        _q_pool = _q_cand[~_q_cand.cid.isin(_q_excl_ids)].reset_index(drop=True)
        chk("资格池可由主名单−剔除表独立重建",
            len(_q_pool) == int(_QJ["pool"]["n_eligible"]),
            f"{len(_q_cand)} − {len(_q_excl_ids)} = {len(_q_pool)}")
        _q_edges = {nm: np.percentile(_q_pool[col].to_numpy(float),
                                      [100 * i / _QK for i in range(1, _QK)])
                    for col, nm in _QAX}
        # 证书按 6 位小数存储 ⇒ 容差 1e-6（仍足以暴露任何上游漂移）
        _q_edge_bad = [nm for _c, nm in _QAX
                       if max(abs(a - b) for a, b in
                              zip(_q_edges[nm], _QOBJ["bin_edges"][nm])) > 1e-6]
        chk("分位边界与证书一致（独立复算）", not _q_edge_bad,
            "；".join(_q_edge_bad) if _q_edge_bad else f"{len(_QAX)} 轴 × {_QK} 层")

        def _q_bin_id(_v, _edges):
            return int(np.searchsorted(np.asarray(_edges), float(_v), side="left"))

        _q_roster = pd.read_csv(_q_csv, low_memory=False)
        _q_roster = _q_roster[_q_roster.candidate_type == "candidate"]
        _q_ids = sorted(str(c) for c in _q_roster.cid)
        _q_pool_ix = _q_pool.set_index("cid")
        _q_binof: dict[tuple[str, int], set] = {}
        for _col, _nm in _QAX:
            _s = _q_pool[["cid", _col]].copy()
            _s["_b"] = [_q_bin_id(v, _q_edges[_nm]) for v in _s[_col].to_numpy(float)]
            for _b in range(_QK):
                _q_binof[(_nm, _b)] = set(_s.loc[_s._b == _b, "cid"].astype(str))
        _q_counts = {k: len(v & set(_q_ids)) for k, v in _q_binof.items()}
        _q_cnt_k = {f"{k[0]}|Q{k[1] + 1}": v for k, v in _q_counts.items()}
        _q_live = [v for k, v in _q_counts.items() if _q_binof[k]]
        _q_t = min(_q_live) if _q_live else 0
        _q_z = np.array([float(_q_pool_ix.loc[c, "CNN_VS"]) for c in _q_ids])
        _q_span, _q_sum = float(_q_z.max() - _q_z.min()), float(_q_z.sum())
        chk("从批次表独立复算 t / 跨度 / Σ 与证书一致",
            _q_t == int(_QOPT["t_star"])
            and abs(_q_span - float(_QOPT["span"])) < 1e-6
            and abs(_q_sum - float(_QOPT["sum_CNN_VS"])) < 1e-6,
            f"t={_q_t}（证书 {_QOPT['t_star']}）、跨度 {round(_q_span, 6)}"
            f"（{_QOPT['span']}）、Σ {round(_q_sum, 6)}（{_QOPT['sum_CNN_VS']}）")
        chk("证书值与规则登记值一致（漂移检测）",
            _q_t == int(_QEXP["t_star"])
            and abs(round(_q_span, 4) - float(_QEXP["span"])) < 1e-9
            and abs(round(_q_sum, 4) - float(_QEXP["sum_CNN_VS"])) < 1e-6
            and int(_QHR["n_swaps"]) == int(_QEXP["n_swaps_vs_heuristic"]),
            f"规则登记 t*={_QEXP['t_star']} 跨度={_QEXP['span']} "
            f"Σ={_QEXP['sum_CNN_VS']} 互换={_QEXP['n_swaps_vs_heuristic']}")

        # ---- 四条硬约束独立实测
        _q_chems = _q_roster.chemotype.value_counts().to_dict()
        _q_risk = int(_q_roster.risk_am1bcc.sum())
        _q_med = float(_q_cand.CNN_VS.median())
        _q_below = int((_q_roster.CNN_VS < _q_med).sum())
        _q_tc_cap = float(_q_rule["caps"]["pairwise_tc_max"])
        _q_bad = []
        if len(_q_ids) != int(_q_rule["n_candidates"]):
            _q_bad.append(f"席位数 {len(_q_ids)}")
        if max(_q_chems.values()) > int(_q_rule["caps"]["chem_cap_value"]):
            _q_bad.append(f"化学型 {_q_chems}")
        if _q_risk > int(_q_rule["caps"]["risk_cap_value"]):
            _q_bad.append(f"风险档 {_q_risk}")
        if _q_below < int(_q_rule["stratification"]["min_below_median_candidates"]):
            _q_bad.append(f"下半区 {_q_below}")
        chk("四条硬约束独立实测均满足且与证书一致",
            not _q_bad and bool(_QOPT["compliance"]["pass"])
            and all(_QOPT["compliance"]["pass"].values()),
            "；".join(_q_bad) if _q_bad else
            f"化学型上限命中 {max(_q_chems.values())}、风险档 {_q_risk}、下半区 {_q_below}、"
            f"批内 Tc 上限 {_QOPT['compliance']['max_pairwise_tc']} < {_q_tc_cap}")

        # ---- 独立重解：自写一遍「max t」的 MILP，验证不存在覆盖下限更高的解
        try:
            from scipy.optimize import (Bounds as _QB, LinearConstraint as _QLC,
                                        milp as _Qmilp)
            _q_n = len(_q_pool)
            _q_all = [str(c) for c in _q_pool.cid]
            _q_pos = {c: i for i, c in enumerate(_q_all)}
            _q_rows, _q_lb, _q_ub = [], [], []
            _r = np.zeros((1, _q_n + 1)); _r[0, :_q_n] = 1.0
            _q_rows.append(_r); _q_lb.append(float(_q_rule["n_candidates"]))
            _q_ub.append(float(_q_rule["n_candidates"]))
            for _ch in sorted(_q_pool.chemotype.unique()):
                _r = np.zeros((1, _q_n + 1))
                _r[0, [_q_pos[str(c)] for c in _q_pool.loc[_q_pool.chemotype == _ch, "cid"]]] = 1.0
                _q_rows.append(_r); _q_lb.append(-np.inf)
                _q_ub.append(float(_q_rule["caps"]["chem_cap_value"]))
            _r = np.zeros((1, _q_n + 1))
            _r[0, :_q_n] = _q_pool.risk_am1bcc.to_numpy(float)
            _q_rows.append(_r); _q_lb.append(-np.inf)
            _q_ub.append(float(_q_rule["caps"]["risk_cap_value"]))
            _r = np.zeros((1, _q_n + 1))
            _r[0, :_q_n] = (_q_pool.CNN_VS.to_numpy(float) < _q_med).astype(float)
            _q_rows.append(_r)
            _q_lb.append(float(_q_rule["stratification"]["min_below_median_candidates"]))
            _q_ub.append(np.inf)
            _q_fp = {}
            _qz = np.load(ROOT / "02_druglikeness/tables/fingerprints_ecfp4.npz",
                          allow_pickle=True)
            for _c, _bits in zip(_qz["cid"], _qz["bits"]):
                _bv = _QEBV(2048)
                _bv.SetBitsFromList(np.flatnonzero(_bits).tolist())
                _q_fp[str(_c)] = _bv
            _q_conf = 0
            for _i in range(_q_n):
                for _j in range(_i + 1, _q_n):
                    if _QDS.TanimotoSimilarity(_q_fp[_q_all[_i]],
                                               _q_fp[_q_all[_j]]) >= _q_tc_cap:
                        _r = np.zeros((1, _q_n + 1)); _r[0, _i] = _r[0, _j] = 1.0
                        _q_rows.append(_r); _q_lb.append(-np.inf); _q_ub.append(1.0)
                        _q_conf += 1
            for _c2, _nm2 in _QAX:
                for _b in range(_QK):
                    _mem = sorted(_q_binof[(_nm2, _b)])
                    if not _mem:
                        continue
                    _r = np.zeros((1, _q_n + 1))
                    _r[0, [_q_pos[c] for c in _mem]] = 1.0
                    _r[0, _q_n] = -1.0
                    _q_rows.append(_r); _q_lb.append(0.0); _q_ub.append(np.inf)
            _qc = np.zeros(_q_n + 1); _qc[_q_n] = -1.0
            _qres = _Qmilp(c=_qc, integrality=np.array([1] * _q_n + [0]),
                           bounds=_QB(np.zeros(_q_n + 1),
                                      np.concatenate([np.ones(_q_n),
                                                      [float(_QJ["objective"]["n_bins"])]])),
                           constraints=_QLC(np.vstack(_q_rows), np.array(_q_lb),
                                            np.array(_q_ub)),
                           options={"presolve": True, "disp": False})
            _q_tmax = round(-float(_qres.fun))
            chk("独立重解：不存在覆盖下限更高的解（t 上界 = t*）",
                int(_qres.status) == 0 and _q_tmax == int(_QOPT["t_star"]) == _q_t,
                f"重解 t_max={_q_tmax}（status={_qres.status}）；证书 t*="
                f"{_QOPT['t_star']}、批次表实算 t={_q_t}；Tc 冲突对 {_q_conf}（非绑定）")
        except ImportError as _e:
            chk("独立重解可用（scipy.optimize.milp）", False, f"import 失败：{_e}")

        # ---- 对照实现留痕与席位角色标签
        _q_hr_ids = set(str(c) for c in _QHR.get("roster", []))
        _q_hz = np.array([float(_q_pool_ix.loc[c, "CNN_VS"]) for c in sorted(_q_hr_ids)])
        _q_hc = {k: len(v & _q_hr_ids) for k, v in _q_binof.items()}
        _q_hlive = [v for k, v in _q_hc.items() if _q_binof[k]]
        _q_sym = set(_q_ids) ^ _q_hr_ids
        chk("对照（v1.6.0 贪心）留痕且目标值可独立复算",
            len(_q_hr_ids) == int(_q_rule["n_candidates"])
            and (min(_q_hlive) if _q_hlive else 0) == int(_QHR["t"])
            and abs(float(_q_hz.max() - _q_hz.min()) - float(_QHR["span"])) < 1e-6
            and abs(float(_q_hz.sum()) - float(_QHR["sum_CNN_VS"])) < 1e-6
            and len(_q_sym) % 2 == 0 and len(_q_sym) // 2 == int(_QHR["n_swaps"]),
            f"t={_QHR['t']} 跨度={_QHR['span']} Σ={_QHR['sum_CNN_VS']}；"
            f"与最优解重叠 {_QHR.get('overlap_with_optimal')}、互换 {len(_q_sym) // 2}")
        _q_roles = _q_roster.selection_reason.value_counts().to_dict()
        _q_roster_ids = set(str(c) for c in _q_roster.cid)
        chk("席位角色标签为 v1.7.0 三类且计数与证书自洽",
            sum(_q_roles.values()) == len(_q_ids)
            and set(_q_roles).issubset({"跨度端点席位（移除此席跨度收窄）",
                                        "覆盖下限约束席位（移除此席 t 会下降）",
                                        "质量兜底填充（移除后 t 与跨度均不变）"})
            and {k: v for k, v in _QOPT["selection_role_counts"].items() if v}
            == dict(_q_roles),
            f"{len(_q_roster_ids)} 席：{dict(_q_roles)}")
        chk("证书分位席位表与独立复算逐格一致（含池内格容量）",
            all(int(v["optimal"]) == _q_cnt_k[k] for k, v in _QOPT["bin_counts"].items()),
            f"{len(_QOPT['bin_counts'])} 格（池容量 " +
            "/".join(str(v["pool_n"]) for v in _QOPT["bin_counts"].values()) + "）；"
            f"最优解占据 {sum(1 for v in _QOPT['bin_counts'].values() if v['optimal'] > 0)} 格")
        _q_deliv = ROOT / "06_deliverables/SOP-05_md_batch.json"
        if _q_deliv.is_file():
            _QD = json.loads(_q_deliv.read_text(encoding="utf-8"))
            chk("交付层证书副本与包内同源（13_assemble 复制链）",
                _QD.get("optimality", {}).get("roster") == _QOPT.get("roster")
                and str(_QD.get("rule_version")) == str(rule["rule_version"]),
                f"rule={_QD.get('rule_version')}；名册 "
                f"{len(_QD.get('optimality', {}).get('roster', []))} 席")
        else:
            chk("交付层证书副本存在", False,
                "缺 06_deliverables/SOP-05_md_batch.json（需先跑 13_assemble）")

    print(f"\n== 结论：{'全部通过' if not FAIL else '未通过 ' + str(len(FAIL)) + ' 项'} ==")
    if FAIL:
        for f in FAIL:
            print(f"  - {f}")
        return 7
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
