#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""输入换版校验：对「全库打分/排名表」做结构与取值契约校验（换数据版专用）。

背景：本包的输入层只有一个文件是「换数据」对象——
  00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv
它的 35 列里有 31 列被下游 11 个脚本按列名直接消费，且第 1 步会做
「SDF 属性 vs 上游表逐值比对」与「|CNNscore × CNNaffinity − CNN_VS| ≤ 5e-4」两条硬校验。
换表前用本脚本把契约一次性跑清楚，比让链在第 3 步炸掉便宜得多。

用法：
  python3 tests/check_swap_csv.py --file <新表.csv>
  python3 tests/check_swap_csv.py --file <新表.csv> \\
      --sdf  00_inputs/raw/inputs/top5000_best_pose.sdf \\
      --herg 03_gates/herg/herg_predictions.csv \\
      --admet 03_gates/admet/raw/CNN_VS_ranking_49963_ADMET.csv \\
      --json reports/swap_check.json
可选：
  --tol-identity 5e-4   恒等式容差（与 config/frozen_rule.json 的口径一致）
  --min-rows N          覆盖 E1 的行数下限（缺省：给了 --sdf 就取该 SDF 的记录数，否则 5,000）
  --allow-identity-break 恒等式不成立时降级为警告（必须同时登记口径变更，否则不要用）
  --allow-new-status     允许 pka_status 出现 {ok, no_basic_center} 之外的新取值（下游会按
                         imputed_0.5 静默处理，用前请确认这是想要的）
  --strict              把 WARN 也当作失败
本脚本只读；不写输入、不改规则、不落包外路径。
退出码：0 = 通过；7 = 有 FAIL；8 = 文件缺失/无法读取。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# ---- 契约常量（取自 SOP_5 输入层的实测：49,963 行 × 35 列）----
COLS = ["rank", "cid", "smiles", "CNN_VS", "CNNscore", "CNNaffinity", "CNNaffinity_variance",
        "minimizedAffinity", "part", "mw", "logp", "hbd", "hba", "tpsa", "rotb", "pains",
        "brenk", "p_RF", "rank_score", "basic_pka", "clogd74", "pka_provider",
        "pka_provider_version", "pka_status", "pka_basic_center_count",
        "pka_basic_center_atom_indices", "pka_basic_center_values", "pka_error",
        "clogd_provider", "clogd_provider_version", "clogd_confidence", "clogd74_lower_80",
        "clogd74_upper_80", "clogd_status", "clogd_error"]
# 下游按列名直接读取、且必须非空的数值列
NUMERIC_REQUIRED = ["rank", "CNN_VS", "CNNscore", "CNNaffinity", "CNNaffinity_variance",
                    "minimizedAffinity", "mw", "logp", "hbd", "hba", "tpsa", "rotb",
                    "p_RF", "rank_score", "basic_pka", "clogd74", "clogd_confidence",
                    "clogd74_lower_80", "clogd74_upper_80", "pka_basic_center_count"]
# 允许整列为空（下游不读其值：pains 由 RDKit 重算；*_error 为上游诊断列）
ALLOWABLE_EMPTY = ["pains", "pka_error", "clogd_error"]
# 允许部分空
ALLOWABLE_PARTIAL_NULL = ["brenk"]
PKA_STATUS_OK = {"ok", "no_basic_center"}
SDF_PROPS = ["CNN_VS", "CNNscore", "CNNaffinity", "CNNaffinity_variance",
             "minimizedAffinity", "rank_score"]

R = []          # 结果行：(级别, 编号, 断言, 详情)


def add(level: str, code: str, msg: str, detail: str = "") -> None:
    R.append({"level": level, "code": code, "check": msg, "detail": detail})


def main() -> int:
    ap = argparse.ArgumentParser(description="SOP_ML 换表契约校验（只读）")
    ap.add_argument("--file", required=True, help="待校验的排名/打分表 CSV")
    ap.add_argument("--sdf", help="本批姿态 SDF（校验 cid 覆盖与属性同源性；也用于定 E1 行数下限）")
    ap.add_argument("--min-rows", type=int, default=None,
                    help="E1 行数下限；缺省 = SDF 记录数（给了 --sdf）或 5000")
    ap.add_argument("--herg", help="hERG 落门表 herg_predictions.csv（校验 cid 覆盖）")
    ap.add_argument("--admet", help="ADMET 冻结表（校验 cid 集合是否一致）")
    ap.add_argument("--json", help="把结果写成 JSON 报告")
    ap.add_argument("--tol-identity", type=float, default=5e-4)
    ap.add_argument("--allow-identity-break", action="store_true")
    ap.add_argument("--allow-new-status", action="store_true")
    ap.add_argument("--strict", action="store_true")
    a = ap.parse_args()

    p = Path(a.file)
    if not p.is_file():
        print(f"[ERR] 文件不存在：{p}")
        return 8
    try:
        df = pd.read_csv(p, low_memory=False)
    except Exception as exc:                                    # noqa: BLE001
        print(f"[ERR] 无法按 CSV 读取：{exc}")
        return 8

    add("INFO", "E0", "读取成功", f"{p.name}：{df.shape[0]} 行 × {df.shape[1]} 列")

    # E2 列集合
    missing = [c for c in COLS if c not in df.columns]
    extra = [c for c in df.columns if c not in COLS]
    if missing:
        add("FAIL", "E2", "35 列齐备", f"缺列 {missing}")
    else:
        add("PASS", "E2", "35 列齐备", "列名逐一匹配（顺序不限）")
    if extra:
        add("WARN", "E2b", "无多余列", f"多出 {extra}（下游忽略，但请确认不是改名残留）")
    # 缺列时仍继续跑「只依赖 cid / smiles / rank」的检查与耦合覆盖检查：
    # 这样一次运行就能同时给出「缺什么列」和「姿态/hERG/ADMET 覆盖如何」，便于一次性判断工作量
    has_cols = not missing

    # E1 行数：下限 = 本批姿态数（表必须能容纳每条姿态的主人）。
    # 原口径写死 5,000（SOP_5 的姿态条数）；换数据后姿态条数随批次变化，
    # 判据改为「≥ 姿态条数」，SDF 未给时退回 5,000。2026-09-24 登记的口径变更。
    n = len(df)
    floor = a.min_rows
    if floor is None:
        if a.sdf and Path(a.sdf).is_file():
            floor = Path(a.sdf).read_text(encoding="utf-8", errors="replace").count("$$$$")
        else:
            floor = 5000
    if n < floor:
        add("FAIL", "E1", f"行数 ≥ {floor}",
            f"仅 {n} 行——本批 {floor} 条姿态的主人必须能在表中找到")
    else:
        add("PASS", "E1", f"行数 ≥ {floor}", f"{n} 行（下限 {floor}）")

    # E3 cid
    cid = df["cid"].astype(str)
    if cid.isna().any() or (cid.str.strip() == "").any() or cid.str.lower().eq("nan").any():
        add("FAIL", "E3", "cid 非空", "存在空值")
    elif not cid.is_unique:
        dup = cid[cid.duplicated()].unique()[:5].tolist()
        add("FAIL", "E3", "cid 唯一", f"重复 {int(cid.duplicated().sum())} 行，例 {dup}")
    else:
        add("PASS", "E3", "cid 非空且唯一", f"{len(cid)} 个")

    # E4 smiles
    s = df["smiles"].astype(str)
    bad_smi = []
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        for i, v in enumerate(s):
            if Chem.MolFromSmiles(v) is None:
                bad_smi.append((i + 2, v[:40]))
            if len(bad_smi) >= 5:
                break
    except Exception as exc:                                     # noqa: BLE001
        add("WARN", "E4", "smiles 可解析（RDKit）", f"RDKit 不可用：{exc}")
    if bad_smi:
        add("FAIL", "E4", "smiles 可解析（RDKit）", f"不可解析行号（含表头）{bad_smi}")
    else:
        add("PASS", "E4", "smiles 全部可解析（RDKit）", "抽查至首个坏例为止")
    dup_s = int(s.duplicated().sum())
    if dup_s:
        add("WARN", "E4b", "smiles 无重复", f"重复 {dup_s} 行；10b0_admet_regenerate 会就此告警")
    else:
        add("PASS", "E4b", "smiles 无重复", "0")

    # E5 rank 为 1..N 的排列
    if "rank" not in df.columns:
        add("WARN", "E5", "rank 为 1..N 的排列", "缺 rank 列，跳过")
        rk = None
    else:
        rk = pd.to_numeric(df["rank"], errors="coerce")
    if rk is not None:
        if rk.isna().any():
            add("FAIL", "E5", "rank 为整数", f"非数值 {int(rk.isna().sum())} 行")
        elif sorted(rk.astype(int).tolist()) != list(range(1, n + 1)):
            add("FAIL", "E5", "rank 为 1..N 的排列", f"min={rk.min()} max={rk.max()} n={n}")
        else:
            add("PASS", "E5", "rank 为 1..N 的排列", f"1..{n}")

    if not has_cols:
        add("WARN", "E6", "列级取值检查", "缺列，E6–E11 跳过（先补齐 35 列再跑）")
    else:
        # E6 必需数值列非空且可解析
        bad_cols = []
        for c in NUMERIC_REQUIRED:
            v = pd.to_numeric(df[c], errors="coerce")
            k = int(v.isna().sum())
            if k:
                bad_cols.append(f"{c}(空/非数值 {k})")
        if bad_cols:
            add("FAIL", "E6", "必需数值列非空可解析", "；".join(bad_cols))
        else:
            add("PASS", "E6", f"{len(NUMERIC_REQUIRED)} 个数值列非空可解析", "全部通过")

        # E6b 允许整列为空的列
        filled = [c for c in ALLOWABLE_EMPTY if int(df[c].notna().sum())]
        add("INFO", "E6b", "允许为空的列", f"其中已有值的：{filled or '无'}")
        partial = {c: int(df[c].isna().sum()) for c in ALLOWABLE_PARTIAL_NULL}
        add("INFO", "E6c", "允许部分空的列", f"{partial}")

        # E7 恒等式 |CNNscore × CNNaffinity − CNN_VS| ≤ tol（包内第 1 步硬校验）
        d = (pd.to_numeric(df["CNNscore"], errors="coerce")
             * pd.to_numeric(df["CNNaffinity"], errors="coerce")
             - pd.to_numeric(df["CNN_VS"], errors="coerce")).abs()
        mx = float(d.max())
        if mx <= a.tol_identity:
            add("PASS", "E7", "|CNNscore×CNNaffinity − CNN_VS| ≤ 容差", f"max={mx:.2e} ≤ {a.tol_identity:g}")
        else:
            lvl = "WARN" if a.allow_identity_break else "FAIL"
            add(lvl, "E7", "|CNNscore×CNNaffinity − CNN_VS| ≤ 容差",
                f"max={mx:.2e}，越界 {int((d > a.tol_identity).sum())} 行——该恒等式是 gnina 输出定义，"
                f"新模型若不满足，第 1 步属性一致性门会失败；确需放行请用 --allow-identity-break 并在 SWAP_SPEC 登记")

        # E8 clogd 区间自洽
        q = df.dropna(subset=["clogd74", "clogd74_lower_80", "clogd74_upper_80"])
        if len(q):
            bad = int((~((q["clogd74_lower_80"] <= q["clogd74"])
                         & (q["clogd74"] <= q["clogd74_upper_80"]))).sum())
            if bad:
                add("FAIL", "E8", "clogd74_lower_80 ≤ clogd74 ≤ clogd74_upper_80", f"越界 {bad} 行")
            else:
                add("PASS", "E8", "cliogd 区间自洽".replace("cliogd", "clogd"), f"{len(q)} 行")
        else:
            add("WARN", "E8", "clogd 区间自洽", "三列全空，跳过")

        # E9 pka_status 取值域（下游只识别 ok / no_basic_center）
        st = set(df["pka_status"].dropna().astype(str).unique())
        unknown = sorted(st - PKA_STATUS_OK)
        if unknown and not a.allow_new_status:
            add("FAIL", "E9", "pka_status ∈ {ok, no_basic_center}", f"未知取值 {unknown}")
        elif unknown:
            add("WARN", "E9", "pka_status 取值域", f"未知取值 {unknown}（下游按 imputed_0.5 静默处理）")
        else:
            add("PASS", "E9", "pka_status ∈ {ok, no_basic_center}", f"分布 {df['pka_status'].value_counts().to_dict()}")

        # E10 provider / confidence
        conf = pd.to_numeric(df["clogd_confidence"], errors="coerce")
        if conf.isna().any() or ((conf <= 0) | (conf > 1)).any():
            add("WARN", "E10", "clogd_confidence ∈ (0,1]", f"range=[{conf.min()},{conf.max()}]")
        else:
            add("PASS", "E10", "clogd_confidence ∈ (0,1]", f"const={sorted(conf.unique())[:3]}")
        for c in ("pka_provider", "pka_provider_version", "clogd_provider", "clogd_provider_version",
                  "clogd_status"):
            k = int(df[c].isna().sum())
            if k:
                add("FAIL", "E10b", f"{c} 非空", f"空 {k} 行")
        add("INFO", "E10c", "来源标识", f"pka={df['pka_provider'].iloc[0]} / clogd={df['clogd_provider'].iloc[0]}")

        # E11 part 命名
        bad_part = [x for x in df["part"].astype(str).unique() if not (x.startswith("part_"))]
        if bad_part:
            add("WARN", "E11", "part 命名 part_NN", f"异常取值 {bad_part[:5]}")
        else:
            add("PASS", "E11", "part 命名 part_NN", f"{df['part'].nunique()} 个分组")

    tc = set(cid)

    # E12 姿态 SDF 覆盖 + 属性同源性（E13 需要属性列）
    if a.sdf:
        sp = Path(a.sdf)
        if not sp.is_file():
            add("WARN", "E12", "姿态 SDF 覆盖", f"未找到 {sp}")
        else:
            from rdkit import Chem, RDLogger
            RDLogger.DisableLog("rdApp.*")
            sup = Chem.SDMolSupplier(str(sp), sanitize=False, removeHs=False)
            rows, missing_cid, mism = [], [], 0
            for m in sup:
                if m is None:
                    continue
                c = m.GetProp("cid") if m.HasProp("cid") else m.GetProp("_Name")
                rows.append(c)
                if c not in tc:
                    missing_cid.append(c)
                    continue
                ref = df.loc[cid[cid == c].index[0]]
                for pr in SDF_PROPS:
                    if m.HasProp(pr):
                        try:
                            if abs(float(m.GetProp(pr)) - float(ref[pr])) > 1e-3:
                                mism += 1
                                break
                        except ValueError:
                            mism += 1
                            break
            if missing_cid:
                add("FAIL", "E12", "SDF 的 cid ⊆ 表的 cid",
                    f"{len(missing_cid)} 个不在表中，例 {missing_cid[:5]}——第 1 步会判为「上游缺失」")
            else:
                add("PASS", "E12", "SDF 的 cid ⊆ 表的 cid", f"{len(rows)} 个姿态全部命中")
            if not has_cols:
                add("WARN", "E13", "SDF 内嵌属性 vs 新表同源", "缺属性列，跳过比对")
            elif mism:
                add("WARN", "E13", "SDF 内嵌属性 vs 新表同源",
                    f"{mism}/{len(rows)} 个姿态的属性与表值不一致（>1e-3）——姿态属性来自旧打分轮，"
                    f"第 1 步的「SDF 属性 vs 上游表逐值比对」会据此判违规；见 SWAP_SPEC §4 的三种处置")
            else:
                add("PASS", "E13", "SDF 内嵌属性 vs 新表同源", "逐值一致（≤1e-3）")

    # E14 hERG 落门表覆盖
    if a.herg:
        hp = Path(a.herg)
        if not hp.is_file():
            add("WARN", "E14", "hERG 落门表覆盖", f"未找到 {hp}")
        else:
            hc = set(pd.read_csv(hp, usecols=["cid"], low_memory=False).cid.astype(str))
            miss = sorted(hc - tc)
            if miss:
                add("FAIL", "E14", "hERG 落门表 cid ⊆ 表的 cid",
                    f"{len(miss)} 个不在表中，例 {miss[:5]}——第 3 步按 cid 左连接会静默丢行")
            else:
                add("PASS", "E14", "hERG 落门表 cid ⊆ 表的 cid", f"{len(hc)} 个全部命中")

    # E15 ADMET 冻结表 cid 集合
    if a.admet:
        dp = Path(a.admet)
        if not dp.is_file():
            add("WARN", "E15", "ADMET 冻结表覆盖", f"未找到 {dp}")
        else:
            ac = set(pd.read_csv(dp, usecols=["cid"], low_memory=False).cid.astype(str))
            only_t, only_a = len(tc - ac), len(ac - tc)
            if only_t or only_a:
                add("WARN", "E15", "ADMET 冻结表 cid 集合一致",
                    f"表有而 ADMET 无 {only_t} 个；ADMET 有而表无 {only_a} 个"
                    f"——需按新表重生成 03_gates/admet/raw 的冻结表，否则 10b/10c 的覆盖会掉")
            else:
                add("PASS", "E15", "ADMET 冻结表 cid 集合一致", f"{len(ac)} 个")

    _dump(a.json)
    return _report(strict=a.strict)


def _dump(path: str | None) -> None:
    if not path:
        return
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"results": R,
                               "n_pass": sum(1 for r in R if r["level"] == "PASS"),
                               "n_fail": sum(1 for r in R if r["level"] == "FAIL"),
                               "n_warn": sum(1 for r in R if r["level"] == "WARN")},
                              ensure_ascii=False, indent=2), encoding="utf-8")


def _report(strict: bool = False) -> int:
    order = {"FAIL": 0, "WARN": 1, "PASS": 2, "INFO": 3}
    for r in sorted(R, key=lambda x: order[x["level"]]):
        if r["level"] == "INFO":
            continue
        print(f'[{r["level"]:4s}] {r["code"]:5s} {r["check"]}' + (f' — {r["detail"]}' if r["detail"] else ""))
    n_fail = sum(1 for r in R if r["level"] == "FAIL")
    n_warn = sum(1 for r in R if r["level"] == "WARN")
    n_pass = sum(1 for r in R if r["level"] == "PASS")
    print(f"\n合计：PASS {n_pass}｜WARN {n_warn}｜FAIL {n_fail}")
    if n_fail:
        print("结论：**不通过** —— 先把 FAIL 项修掉（或按 SWAP_SPEC 显式登记口径变更）再接入链路")
        return 7
    if n_warn and strict:
        print("结论：**严格模式下不通过**（有 WARN）")
        return 7
    print("结论：通过（WARN 项请逐条确认是否属预期口径）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
