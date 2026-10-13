#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 3b：hERG 门的外部接口。

本包**不含** hERG 预测实现（本地无 ADMET-AI / chemprop 运行时，权重不可得）。
因此把该门做成显式接口：
  ① 生成「四臂并集」补算清单 → 03_gates/herg/SOP-05_herg_input_union.csv
     （四臂 = {QED≥0.6, ≥0.5, ≥0.4}×anchor=4 ∪ {QED≥0.6}×anchor≥3，各先过
       PAINS0/Ro5/Veber0 与姿态闸门 pass、maxTc<0.9）
  ② 若外部把预测结果写到 03_gates/herg/herg_predictions.csv（列 cid,hERG），
     本步自动落门（hERG < 0.5）并登记覆盖度；
  ③ 否则记为 pending：不删分子，但名单与报告全部标注 herg_gate=pending。
输出：03_gates/herg/SOP-05_herg_input_union.csv、HERG_INTERFACE.md
      03_gates/reports/herg_status.json
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

from vendor.io_utils import StepLog, read_json, write_csv, write_json, write_text
from vendor.master import build_master

ARMS = [("QED0.60_anchor4", 0.60, 4), ("QED0.50_anchor4", 0.50, 4), ("QED0.40_anchor4", 0.40, 4),
        ("QED0.60_anchor3", 0.60, 3)]

TEMPLATE = """# hERG 门外部接口（SOP_5）

## 为什么需要外部接口

本工程包内**没有** hERG 预测实现：本地环境未安装 `admet_ai` / `chemprop` / `torch`
等运行时，且 ADMET-AI 权重不在包内。为保持"包内自包含、可复现"，本步只生成
**待补算清单**并预留落门位置，不伪造预测值。

## 需要外部产出的文件

路径（包内相对）：`03_gates/herg/herg_predictions.csv`

| 列 | 类型 | 说明 |
|---|---|---|
| `cid` | str | 与清单 `SOP-05_herg_input_union.csv` 的 `cid` 一致 |
| `hERG` | float | hERG 阻断概率，0–1，越大越危险（ADMET-AI 口径） |
| `source`（可选） | str | 模型名与版本，例如 `ADMET-AI v1.x / chemprop-2.x` |
| `weight_sha256`（可选） | str | 权重文件哈希，用于固化 |

参考实现（外部执行，命令示例）：

```bash
# ADMET-AI（Chemprop 后端）；需联网下载权重，建议在包外独立环境执行
pip install admet-ai
admet_predict --input-file SOP-05_herg_input_union.csv --output-file herg_predictions.csv
# 该命令输出的 hERG 列名为 "hERG"（阻断概率）
```

## 落门规则

`hERG < 0.5` 通过（阈值见 `config/frozen_rule.json → step3_gates.herg.threshold`）。
缺值**不得**当作「判为高风险」：状态记为 `not_evaluated`（未评估），并在 `herg_status` 与
报告中同时报出**作用域覆盖率**（`scope_coverage = 并集内分子数 / 池规模`）。
把缺值写成「保守不过门」会把「没测」冒名为「不合格」，属口径混淆。

## 外部产出溯源（登记于 `config/frozen_rule.json → step3_gates.herg`）

- 落门表：`{pred_path}`（md5 `{pred_md5}`、sha256 `{pred_sha256}`）
- 来源：{source}
- 外部队列：{external_run}

## 当前状态

- 清单规模：**{n_scope}** 个分子（四臂并集；规则期望 {expected}）
- 预测表是否存在：**{exists}**
- 状态：**{status}**
- 作用域覆盖率：**{scope_coverage}**（{n_scope}/{n_pool}；池内未评估 {n_not_evaluated} 个）
- 为什么作用域小于池：hERG 外部补算的范围是四臂并集（代价控制），而非整个池。
  池内并集外的分子**没有预测值**，其门状态是 `not_evaluated`，不是「高风险」。
  当前并集外的分子均不在任何交付名单中，但下游若扩大候选池必须同步扩大补算范围，
  否则会把「未评估」当作「不合格」使用。

## 模型依赖度（第二模型交叉验证 · 分析步 08b）

本接口的 `hERG` 列来自**单一模型**（ADMET-AI）。为量化该判定的模型依赖度，链上分析步
`08b_herg_cross_model` 用第二个**独立训练**的 hERG 模型（ChEMBL `CHEMBL240` 活性数据自建，
特征仅 ECFP4）复算同一批四臂并集分子，给出阈值一致性、与阈值无关的排序一致性，
以及对池/名单的实测影响。产物见 `03_gates/herg/cross_model/`（`HERG_CROSS_MODEL.md`、
`reports/herg_cross_model.json`、`tables/`、`figures/`）。

要点（数字见上述产物与最终报告 §4b）：两模型在 0.5 阈值上的门决策一致率仅三分之一到一半，
且**与阈值无关**的排序交集同样低，说明差异不只是标定问题；因此
「hERG 过门」应转述为「在 ADMET-AI 下过门」，不宜当作可迁移的硬判据。
独立性的诚实边界：ADMET-AI 的 hERG 头同样基于公开 hERG 数据（含 TDC 集，多源自 ChEMBL）训练，
故该步量化的是**算法/训练管道依赖度**，不是数据独立的外部验证。

## 降级声明（表缺失时）

`hERG` 门记为 `pending`：名单**不因 hERG 被删减**，所有名单与决策卡带
`herg_gate=pending` 标注。此时：

1. 名单在"结构与活性相似度"两轴上仍是严格达标的（anchor/姿态/maxTc 均实测）；
2. hERG 风险需在 MD 前置的人工复核环节补齐，不得视为已通过；
3. 报告中 `hERG 区间` 一栏显示 `pending（未补算）` 而非区间值。
"""


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    herg = rule["step3_gates"]["herg"]
    thr = herg["threshold"]
    scope_mode = str(herg.get("scope", "union4"))
    with StepLog(ROOT, "08_herg") as log:
        m = build_master(ROOT, log=log)
        base = m[(m.druglike_pass == 1) & (m.pose_gate == "pass") & (m.maxTc_actives < 0.9)]
        scopes = {}
        for name, q, a in ARMS:
            sub = base[(base.qed >= q) & (base.anchor_hits_4A >= a)]
            scopes[name] = set(sub.cid)
            log.info(f"臂 {name}：池 {len(sub)}")
        union = sorted(set().union(*scopes.values()))
        expected = herg.get("expected_scope_size")
        log.info(f"四臂并集（审计件）= {len(union)}（规则期望 {expected}）")

        flags = {name: set(s) for name, s in scopes.items()}
        scope = m[m.cid.isin(union)].copy()
        scope["arm_membership"] = scope.cid.map(
            lambda c: ";".join(name for name, s in flags.items() if c in s))
        write_csv(ROOT / "03_gates/herg/SOP-05_herg_input_union.csv",
                  scope[["cid", "smiles", "qed", "anchor_hits_4A", "maxTc_actives",
                         "MW", "cLogP", "TPSA", "chemotype", "arm_membership"]])

        # 门作用域口径：scope=pool → 全池落门（依据包内冻结表，零新计算）；
        # scope=union4 → 旧口径（只覆盖四臂并集，并集外为未评估）。
        n_pool = int(len(m))
        gate_scope_n = n_pool if scope_mode == "pool" else len(union)
        n_outside = n_pool - gate_scope_n
        scope_meta = {
            "scope_mode": scope_mode,
            "n_pool": n_pool,
            "gate_scope_n": gate_scope_n,
            "n_in_scope": gate_scope_n,
            "n_not_evaluated": n_outside,
            "scope_coverage": round(gate_scope_n / max(n_pool, 1), 4),
            "audit_union_n": len(union),
            "not_evaluated_policy": "缺预测值 → herg_gate_state='not_evaluated'（不当作不合格）",
        }
        pred = ROOT / herg["external_table"]
        status = {
            "status": "applied" if bool(m.herg_prob.notna().any()) else "pending",
            "threshold": thr, "scope_n": len(union), "expected_scope_n": expected,
            "scope_size_match": bool(expected is None or len(union) == expected),
            "arm_pool_sizes": {name: len(s) for name, s in scopes.items()},
            "predictions_file": str(pred.relative_to(ROOT)),
            "predictions_present": bool(pred.is_file()),
            "coverage_in_scope": None,
            "n_pass_in_scope": None,
            **scope_meta,
            "note": ("hERG 表存在 → 已落门" if pred.is_file() else
                     "hERG 表缺失 → pending：不删分子，名单/报告显式标注偏差"),
        }
        # 全池口径下：落门依据是包内冻结表，外部表作对账件（逐值比对）
        _ov = (m.dropna(subset=["herg_prob", "herg_prob_external"])
               if "herg_prob_external" in m.columns else m.iloc[0:0])
        if len(_ov):
            status["external_vs_gate_overlap_n"] = int(len(_ov))
            status["external_vs_gate_max_abs_diff"] = float(
                (_ov.herg_prob - _ov.herg_prob_external).abs().max())
            log.info(f"落门依据 vs 外部表对账：交集 {len(_ov)}，max|Δ| = "
                     f"{status['external_vs_gate_max_abs_diff']:.2e}")
        status["gate_source"] = ("包内冻结 ADMET 全量表（全池口径）" if scope_mode == "pool"
                                 else "外部 hERG 落门表（四臂并集口径）")
        if pred.is_file():
            hp = pd.read_csv(pred, low_memory=False)
            col = "hERG" if "hERG" in hp.columns else [c for c in hp.columns
                                                       if "herg" in c.lower()][0]
            keep = ["cid", col] + [c for c in ("source", "weight_sha256", "model_version")
                                   if c in hp.columns]
            hp = hp[keep].rename(columns={col: "hERG"})
            def _first(col):
                if col not in hp.columns:
                    return None
                v = hp[col].iloc[0]
                return None if pd.isna(v) or str(v).strip() == "" else str(v)

            status["source"] = _first("source")
            status["weight_sha256"] = _first("weight_sha256")   # 应用接口不回传权重哈希 → None
            joined = scope[["cid"]].merge(hp[["cid", "hERG"]], on="cid", how="left")
            if scope_mode == "pool":
                # 全池口径：门内覆盖按**落门依据**（冻结表）统计；
                # 外部落门表对审计用四臂并集的覆盖另列，避免把「外部表没覆盖」误读成「门内有未评估」。
                _gate_v = m.herg_prob
                status["coverage_in_scope"] = round(float(_gate_v.notna().mean()), 4)
                status["n_pass_in_scope"] = int((_gate_v < thr).sum())
                _un_ext = m[m.cid.isin(union)].herg_prob_external \
                    if "herg_prob_external" in m.columns else m.iloc[0:0].herg_prob
                status["audit_union_n"] = int(len(union))
                status["external_table_coverage_of_audit_union"] = round(
                    float(_un_ext.notna().mean()), 4)
                status["n_audit_union_outside_external"] = int(len(union) - _un_ext.notna().sum())
                log.info(f"全池口径：门内覆盖 {status['coverage_in_scope']:.2%}；"
                         f"外部落门表对四臂并集覆盖 "
                         f"{status['external_table_coverage_of_audit_union']:.2%}"
                         f"（并集 {len(union)}，其中 "
                         f"{status['n_audit_union_outside_external']} 个外部表未覆盖，"
                         f"其值取自包内冻结表）")
            else:
                status["coverage_in_scope"] = round(float(joined.hERG.notna().mean()), 4)
                status["n_pass_in_scope"] = int((joined.hERG < thr).sum())
            if scope_mode == "pool":
                status["n_fail_in_scope"] = int((m.herg_prob >= thr).sum())
            else:
                status["n_fail_in_scope"] = int((joined.hERG >= thr).sum())
            # 池层覆盖：并集之外的分子没有预测值；它们本来就进不了任何臂/名单，
            # 但必须显式报数，避免把「保守判为不过门」误读成「被判为高风险」
            in_pool_missing = m[m.herg_prob.isna()]
            status["n_pool_total"] = int(len(m))
            status["n_pool_outside_scope"] = int(len(in_pool_missing))
            _slab = ("门作用域（全池）内覆盖" if scope_mode == "pool" else "并集内覆盖")
            status["note"] = (f"hERG 表存在 → 已落门（作用域 {scope_mode}）；{_slab} "
                              f"{status['coverage_in_scope']:.2%}，通过 "
                              f"{status['n_pass_in_scope']}、不通过 {status['n_fail_in_scope']}；"
                              f"池内 {len(in_pool_missing)} 个分子没有预测值，"
                              f"其门状态为 **not_evaluated（未评估）**，不得转述为「高风险」或"
                              f"「不合格」；作用域覆盖率 "
                              f"{status['scope_coverage']:.2%}")
            log.info(f"hERG 覆盖 {status['coverage_in_scope']:.2%}；范围内通过 "
                     f"{status['n_pass_in_scope']}、不通过 {status['n_fail_in_scope']}；"
                     f"并集外池内分子 {len(in_pool_missing)}（not_evaluated，非「判为高风险」）；"
                     f"作用域覆盖率 {status['scope_coverage']:.2%}")
            log.info(f"hERG 来源：{status.get('source')}｜权重哈希 "
                     f"{str(status.get('weight_sha256'))[:16] + '…' if status.get('weight_sha256') else '未提供（应用接口不回传）'}")
        write_json(ROOT / "03_gates/reports/herg_status.json", status)

        write_text(ROOT / "03_gates/herg/HERG_INTERFACE.md",
                   TEMPLATE.format(n_scope=len(union), expected=expected,
                                   exists=pred.is_file(), status=status["status"],
                                   scope_coverage=f"{status['scope_coverage']:.2%}",
                                   n_pool=status["n_pool"],
                                   n_not_evaluated=status["n_not_evaluated"],
                                   pred_path=str(rule["step3_gates"]["herg"]["external_table"]),
                                   pred_md5=str(rule["step3_gates"]["herg"].get("external_md5") or "未登记"),
                                   pred_sha256=str(rule["step3_gates"]["herg"].get("external_sha256") or "未登记"),
                                   source=str(status.get("source") or "未标注"),
                                   external_run=str(rule["step3_gates"]["herg"].get("external_run") or "未登记")))
        log.info(f"补算清单 → 03_gates/herg/SOP-05_herg_input_union.csv；"
                 f"状态 = {status['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
