"""CNS MPO（Wager TT, Hou X, Verhoest PR, Villalobos A. ACS Chem Neurosci 2010;1(6):435-449,
DOI 10.1021/cn100008c, Table 1）的确定性实现。

六个单属性转换函数 T0（0~1），等权求和 → 0~6：
  cLogP  单调降：≤3 → 1.0；>5 → 0.0；3~5 线性
  cLogD  单调降：≤2 → 1.0；>4 → 0.0；2~4 线性
  MW     单调降：≤360 → 1.0；>500 → 0.0；360~500 线性
  TPSA   峰形  ：40 < TPSA ≤ 90 → 1.0；(20,40] 线性升；(90,120] 线性降；≤20 或 >120 → 0.0
  HBD    单调降：≤0.5 → 1.0；>3.5 → 0.0；0.5~3.5 线性
  pKa    单调降：≤8 → 1.0；>10 → 0.0；8~10 线性
原文 Table 2/3 的示例值作为回归夹具（tests/fixtures/cns_mpo_reference.json）。
"""
from __future__ import annotations

import math

# 原文 Table 3 三个 Pfizer 候选（属性值取自正文）与其报告分数
REFERENCE_CASES = [
    {"name": "PF-02545920 (PDE10)", "clogp": 3.8, "clogd": 3.5, "mw": 392.5,
     "tpsa": 52.8, "hbd": 0, "pka": 4.3, "reported_mpo": 4.6},
    {"name": "PF-03654746 (H3)", "clogp": 2.4, "clogd": 0.0, "mw": 322.4,
     "tpsa": 32.3, "hbd": 1, "pka": 9.2, "reported_mpo": 4.9},
    {"name": "PF-04447943 (PDE9)", "clogp": -1.5, "clogd": -0.7, "mw": 395.4,
     "tpsa": 101.9, "hbd": 1, "pka": 7.9, "reported_mpo": 5.2},
]


def _down(x, lo, hi):
    """单调降：x ≤ lo → 1.0；x > hi → 0.0；区间内线性。"""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 0.5
    if x <= lo:
        return 1.0
    if x > hi:
        return 0.0
    return (hi - x) / (hi - lo)


def _hump(x, lo_bad, lo_good, hi_good, hi_bad):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return 0.5
    if x <= lo_bad or x > hi_bad:
        return 0.0
    if x <= lo_good:
        return (x - lo_bad) / (lo_good - lo_bad)
    if x <= hi_good:
        return 1.0
    return (hi_bad - x) / (hi_bad - hi_good)


def t_clogp(v):
    return _down(v, 3.0, 5.0)


def t_clogd(v):
    return _down(v, 2.0, 4.0)


def t_mw(v):
    return _down(v, 360.0, 500.0)


def t_tpsa(v):
    return _hump(v, 20.0, 40.0, 90.0, 120.0)


def t_hbd(v):
    return _down(v, 0.5, 3.5)


def t_pka(v):
    return _down(v, 8.0, 10.0)


def cns_mpo(clogp, clogd, mw, tpsa, hbd, pka, pka_status: str | None = None) -> dict:
    """返回总分与六个 T0 分项。

    pKa 缺失（如 pka_status == 'no_basic_center' 且无数值）时的处置：
      记 T0=1.0 并在 pka_source 标注 "no_basic_center"（无碱性中心 → 该轴不构成风险）。
      若数值缺失但状态未知 → T0=0.5，pka_source="imputed_0.5"。
    """
    if pka is None or (isinstance(pka, float) and math.isnan(pka)):
        pka_val, pka_src = 1.0, "no_basic_center" if pka_status == "no_basic_center" else "imputed_0.5"
    else:
        pka_val, pka_src = t_pka(pka), "measured_input"
    parts = {
        "T0_clogp": t_clogp(clogp), "T0_clogd": t_clogd(clogd), "T0_mw": t_mw(mw),
        "T0_tpsa": t_tpsa(tpsa), "T0_hbd": t_hbd(hbd), "T0_pka": pka_val,
    }
    total = sum(parts.values())
    return {"cns_mpo": round(total, 4), "pka_source": pka_src,
            **{k: round(v, 4) for k, v in parts.items()}}


def cns_mpo_bounded(clogp, clogd, clogd_lo80, clogd_hi80, mw, tpsa, hbd, pka,
                    pka_status: str | None = None) -> dict:
    """点估计 + **cLogD 80% 预测区间**界定下的 CNS MPO 上下界（标注用，不作门）。

    口径（与 Wager 2010 Table 1 一致）：cLogD 越大该轴越差，故
      cns_mpo_lo80 = 悲观下界：cLogD 取区间上界 clogd_hi80
      cns_mpo_hi80 = 乐观上界：cLogD 取区间下界 clogd_lo80
    pKa 轴本身没有区间（MolGpKa 只给点值），按点值处理。
    """
    # 点估计一律用上游 clogd74 点值（不用区间中点），保证与既有口径逐位一致
    base = cns_mpo(clogp, clogd, mw, tpsa, hbd, pka, pka_status)
    if clogd_lo80 is None or clogd_hi80 is None:
        return {**base, "cns_mpo_lo80": None, "cns_mpo_hi80": None,
                "cns_mpo_pi80_span": None, "pi80_source": "missing_clogd_interval"}
    t_hi, t_lo = t_clogd(clogd_hi80), t_clogd(clogd_lo80)
    tot_lo = base["cns_mpo"] - base["T0_clogd"] + t_hi      # 悲观：取区间上界
    tot_hi = base["cns_mpo"] - base["T0_clogd"] + t_lo      # 乐观：取区间下界
    return {**base,
            "cns_mpo_lo80": round(tot_lo, 4), "cns_mpo_hi80": round(tot_hi, 4),
            "cns_mpo_pi80_span": round(tot_hi - tot_lo, 4),
            "T0_clogd_lo80": round(t_hi, 4), "T0_clogd_hi80": round(t_lo, 4),
            "pi80_source": "cpLogD_v2.0_80pct_interval"}


def regression_check(verbose: bool = False) -> dict:
    """对原文 Table 2/3 的示例回算，返回最大偏差。"""
    diffs = []
    for c in REFERENCE_CASES:
        got = cns_mpo(c["clogp"], c["clogd"], c["mw"], c["tpsa"], c["hbd"], c["pka"])
        d = abs(got["cns_mpo"] - c["reported_mpo"])
        diffs.append({"name": c["name"], "recomputed": got["cns_mpo"],
                      "reported": c["reported_mpo"], "abs_diff": round(d, 4)})
        if verbose:
            print(f"  {c['name']:22s} 回算 {got['cns_mpo']:.4f} vs 原文 {c['reported_mpo']} "
                  f"(|Δ|={d:.4f})")
    return {"cases": diffs, "max_abs_diff": round(max(x["abs_diff"] for x in diffs), 4)}


if __name__ == "__main__":
    print(regression_check(verbose=True))
