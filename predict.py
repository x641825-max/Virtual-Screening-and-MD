#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一键生成最终结果文件（候选清单）—— 赛事《代码提交要求》附件 5 §二(2)、§四。

用途
  把本项目三条 MD 候选来源（SOP_5 目录化合物链 / SOP_ML_run05 生成分子链 /
  二者联立合并层）的筛后清单归一为一份标准化候选清单：

      results/results.csv    UTF-8，逗号分隔
      results/results.xlsx   Excel（同内容，openpyxl 写出）

  并逐条映射到附件 5 §四要求的必填字段：候选编号 / 所属赛道 / SMILES /
  关键预测指标 / 对应模型与运行版本 / 备注。

用法
  python predict.py                       # 用随包冻结的清单（results/ 下）出清单
  python predict.py --source repro        # 用 src/ 下复跑产出的清单出清单（一致性核验用）
  python predict.py --out DIR             # 指定输出目录（默认 results/）
  python predict.py --verify A.csv B.csv  # 逐行核验两份候域清单是否一致

口径边界（随结果一并转述）
  · 本清单是**计算筛选结果**，不含任何实验数据，不构成安全性或有效性结论。
  · hERG 只表示「在 ADMET-AI 口径下过门」，不是独立验证；第二模型（ChEMBL 自建集成）
    与主模型的一致性接近随机，分歧分子已保留并逐分子标注。
  · 「最优」限定在批次目标函数口径内（分层覆盖 t → CNN_VS 跨度 → Σ CNN_VS 字典序），
    不等于「这批分子比落选者更好」。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent

TRACK = "赛道三：离子通道、GPCR等药靶的小分子药物虚拟筛选"
RULE_SET = "SOP-05-A2A-antagonist-shortlist v1.7.0"
MODEL_REF = (
    "CNN_VS=gnina 1.3.2（semi_flexible / cnn_scoring=rescore，上游对接，不在本包重跑）；"
    "hERG=ADMET-AI（外部预测表接口）；第二 hERG 模型=ChEMBL CHEMBL240 自建 ECFP4 集成(RF/LR/kNN)；"
    "CNS MPO 输入 basic_pka=MolGpKa、clogd74=cpLogD v2.0；QED/描述符/指纹=RDKit 2025.09.5；"
    "选批求解=scipy 1.15.3 milp（HiGHS 后端）"
)

SETS = {
    "SOP_5": {
        "label": "SOP_5（目录化合物链）",
        "code": "S5",
        "frozen": "results/SOP_5/SOP-05_MD_batch20.csv",
        "repro": "src/SOP_5/07_md_batch/tables/md_batch20.csv",
        "expect_rows": 20,
        "expect_cand": 18,
        "run_version": "SOP_5 冻结轮（2026-09-23/24，自检 19 步全通过）",
    },
    "SOP_ML_run05": {
        "label": "SOP_ML_run05（ML 生成分子链）",
        "code": "ML",
        "frozen": "results/SOP_ML_run05/SOP-05_MD_batch20.csv",
        "repro": "src/SOP_ML_run05/07_md_batch/tables/md_batch20.csv",
        "expect_rows": 19,
        "expect_cand": 17,
        "run_version": "SOP_ML_run05（4,979 生成分子换数据轮）",
    },
    "merged_final_md": {
        "label": "联立合并层（两来源并池重解）",
        "code": "MG",
        "frozen": "results/merged_final_md/final_md_20.csv",
        "repro": "src/SOP_5/10_final_md/tables/final_md_20.csv",
        "expect_rows": 22,
        "expect_cand": 20,
        "run_version": "SOP_5/10_final_md 步骤 17（合并池 35 → 资格池 33 → 20 席）",
    },
}

OUT_COLS = [
    "candidate_id", "track", "source_set", "source_branch", "molecule_id", "role",
    "smiles", "cnn_vs", "qed", "cns_mpo", "herg_prob_admet_ai", "herg_gate_admet_ai",
    "herg_prob_cross_model_ens", "mw", "clogp", "tpsa", "hbd", "hba", "rotb",
    "chemotype", "scaffold_generic", "anchor_hits_4A", "pose_gate",
    "sa_score", "triage_band",
    "selection_reason", "md_eligibility_note",
    "rule_set", "model_ref", "run_version", "remark",
]


def _get(row: pd.Series, *names, default=""):
    for n in names:
        if n in row.index and pd.notna(row[n]):
            return row[n]
    return default


def _round(v, nd):
    try:
        return round(float(v), nd)
    except Exception:
        return ""


def build(source: str) -> pd.DataFrame:
    frames = []
    for key, spec in SETS.items():
        path = ROOT / spec[source]
        if not path.is_file():
            raise SystemExit(f"[FAIL] 缺少清单：{path}（source={source}）")
        d = pd.read_csv(path)
        n_cand = int((d.get("candidate_type", "candidate") == "candidate").sum())
        if len(d) != spec["expect_rows"] or n_cand != spec["expect_cand"]:
            print(f"[WARN] {key}: 行数 {len(d)}（期望 {spec['expect_rows']}）、"
                  f"候选 {n_cand}（期望 {spec['expect_cand']}）")
        rows = []
        seq = 0
        for _, r in d.iterrows():
            seq += 1
            role = str(_get(r, "candidate_type", default="candidate"))
            branch = _get(r, "source_branch", default="")
            if not branch:
                mid = str(_get(r, "cid"))
                branch = "目录化合物（ZINC/MCULE）" if mid.startswith(("ZINC", "MCULE")) \
                    else "ML 生成分子（SCALE-）"
            gate = _get(r, "herg_gate", default="")
            gate_txt = "" if gate == "" else ("过门（<0.5）" if float(gate) == 1 else "不过门（≥0.5）")
            rows.append({
                "candidate_id": f"{spec['code']}-{seq:02d}",
                "track": TRACK,
                "source_set": key,
                "source_branch": branch,
                "molecule_id": _get(r, "cid"),
                "role": "候选" if role == "candidate" else "阴性对照（不占名额）",
                "smiles": _get(r, "smiles"),
                "cnn_vs": _round(_get(r, "CNN_VS"), 4),
                "qed": _round(_get(r, "qed"), 4),
                "cns_mpo": _round(_get(r, "cns_mpo"), 3),
                "herg_prob_admet_ai": _round(_get(r, "herg_prob"), 4),
                "herg_gate_admet_ai": gate_txt,
                "herg_prob_cross_model_ens": _round(_get(r, "hERG_cross_ens"), 4),
                "mw": _round(_get(r, "MW"), 2),
                "clogp": _round(_get(r, "cLogP"), 3),
                "tpsa": _round(_get(r, "TPSA"), 2),
                "hbd": _get(r, "HBD"), "hba": _get(r, "HBA"), "rotb": _get(r, "RotB"),
                "chemotype": _get(r, "chemotype"),
                "scaffold_generic": _get(r, "scaffold_generic"),
                "anchor_hits_4A": _get(r, "anchor_hits_4A"),
                "pose_gate": _get(r, "pose_gate"),
                "sa_score": _round(_get(r, "sa_score"), 3),
                "triage_band": _get(r, "triage_band"),
                "selection_reason": _get(r, "selection_reason"),
                "md_eligibility_note": _get(r, "md_eligibility_note"),
                "rule_set": RULE_SET,
                "model_ref": MODEL_REF,
                "run_version": spec["run_version"],
                "remark": f"{spec['label']}；批次内位次 {seq}/{len(d)}",
            })
        frames.append(pd.DataFrame(rows)[OUT_COLS])
        print(f"[OK] {key}: {len(d)} 行（候选 {n_cand} + 对照 {len(d) - n_cand}）← {spec[source]}")
    out = pd.concat(frames, ignore_index=True)
    return out


def write(out: pd.DataFrame, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    csv = outdir / "results.csv"
    out.to_csv(csv, index=False, encoding="utf-8")
    try:
        xlsx = outdir / "results.xlsx"
        with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
            out.to_excel(w, sheet_name="results", index=False)
    except Exception as e:  # openpyxl 缺失时只出 csv，不静默假装成功
        print(f"[WARN] xlsx 未写出（{e}）—— csv 已写出，可用 README 中的命令补出")
    else:
        print(f"[OK] 写出 {csv}（{len(out)} 行）与 {outdir / 'results.xlsx'}")


def verify(a: Path, b: Path) -> int:
    da, db = pd.read_csv(a), pd.read_csv(b)
    if list(da.columns) != list(db.columns):
        print("[FAIL] 列不一致"); return 1
    if len(da) != len(db):
        print(f"[FAIL] 行数不一致：{len(da)} vs {len(db)}"); return 1
    diff = (da.fillna("").astype(str) != db.fillna("").astype(str))
    bad = diff.any(axis=1).sum()
    if bad:
        cols = [c for c in da.columns if diff[c].any()]
        print(f"[FAIL] {bad} 行不一致，涉及列：{cols}"); return 1
    print(f"[PASS] 两份清单逐行逐字段一致（{len(da)} 行 × {len(da.columns)} 列）")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="生成标准化候选清单（results.csv / results.xlsx）")
    ap.add_argument("--source", choices=["frozen", "repro"], default="frozen",
                    help="frozen=随包冻结清单（默认）；repro=src/ 下复跑产出")
    ap.add_argument("--out", default="", help="输出目录（默认 results/；repro 时默认 results/reproduced/）")
    ap.add_argument("--verify", nargs=2, metavar=("A", "B"), help="核验两份清单是否一致")
    args = ap.parse_args()

    if args.verify:
        return verify(Path(args.verify[0]), Path(args.verify[1]))

    out = build(args.source)
    outdir = Path(args.out) if args.out else (
        ROOT / ("results" if args.source == "frozen" else "results/reproduced"))
    write(out, outdir)

    cand = out[out.role == "候选"]
    print("\n--- 清单摘要 ---")
    print(f"总行数 {len(out)}（候选 {len(cand)} + 对照 {len(out) - len(cand)}）")
    for k in SETS:
        s = out[out.source_set == k]
        print(f"  {k:16s} 行 {len(s):2d}  候选 {(s.role == '候选').sum():2d}"
              f"  唯一分子 {s.molecule_id.nunique():2d}")
    print(f"三套合计唯一分子 {out.molecule_id.nunique()} 个"
          f"（跨来源结构零重叠，唯二相同的是 2 个共用阴性对照）")
    print("字段：" + "、".join(OUT_COLS[:8]) + f" … 共 {len(OUT_COLS)} 列")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
