"""候选主表构建：把描述符、姿态闸门、活性相似度、上游 pKa/cLogD、hERG 合到一张表。

所有下游脚本（08/09/10/11/13）都只通过本模块取数，口径单点维护。
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import chem_core as cc
from .io_utils import read_json

UPSTREAM_COLS = ["cid", "basic_pka", "clogd74", "pka_status", "pka_basic_center_count",
                 "clogd_status", "clogd_confidence", "clogd74_lower_80", "clogd74_upper_80",
                 "pka_provider", "pka_provider_version",
                 "clogd_provider", "clogd_provider_version"]


def build_master(root: Path | str, log=None) -> pd.DataFrame:
    root = Path(root)
    rule = read_json(root / "config/frozen_rule.json")
    des = pd.read_csv(root / "02_druglikeness/tables/descriptors.csv", low_memory=False)
    gate = pd.read_csv(root / "01_pose_qc/tables/pose_gate.csv", low_memory=False)
    mt = pd.read_csv(root / "03_gates/tables/maxTc_actives.csv", low_memory=False)
    pool = pd.read_csv(root / "01_pose_qc/tables/pose_pool.csv", low_memory=False)
    up = pd.read_csv(root / "00_inputs/raw/inputs/CNN_VS_ranking_49963_with_pKa_cLogD74.csv",
                     low_memory=False)

    # 上游表是 pKa/cLogD 及其型号/区间的唯一来源：先去掉描述符表里的同名列，
    # 否则二次合并会产生 basic_pka_x/_y 之类重复列，下游按「列存在才写入」的取数
    # 会静默丢列（曾导致交付决策卡缺 basic_pka/clogd74）。
    dup = [c for c in UPSTREAM_COLS if c != "cid" and c in des.columns]
    if dup:
        des = des.drop(columns=dup)

    m = (des.merge(gate[["cid", "pose_gate", "anchor_hits_4A", "outside_box",
                         "n_contacts_4A", "warn_or_reject_reason"]], on="cid", how="left")
            .merge(mt[["cid", "maxTc_actives", "nearest_active_id",
                       "pass_maxTc_lt_09", "pass_maxTc_lt_07"]], on="cid", how="left")
            # stereo_source 同样已由 05 写入描述符表，只保留一列，避免 _x/_y 重复列
            .merge(pool[["cid"] + [c for c in ("protonation_note", "stereo_source")
                                   if c not in des.columns]], on="cid", how="left")
            .merge(up[UPSTREAM_COLS], on="cid", how="left"))

    herg_path = root / rule["step3_gates"]["herg"]["external_table"]
    hspec = rule["step3_gates"]["herg"]
    thr = float(hspec["threshold"])
    scope_mode = str(hspec.get("scope", "union4"))
    m["herg_gate_sources"] = ""

    if scope_mode == "pool":
        # 门作用域 = 全池：依据是**包内已冻结**的 ADMET-AI 全量表（49,963 行，覆盖池 100%），
        # 外部落门表作为对账件保留。此前作用域是四臂并集，池内并集外的分子没有预测值 →
        # 门状态是「未评估」，却被写成「保守判为不过门」（把「没测」当成「测出不合格」），
        # 且把两条 anchor≥3 敏感性臂池静默裁小。全池口径不需要新计算（表已在包内并登记哈希）。
        layer = rule["step3_gates"].get("admet_layer", {})
        frozen_p = root / str(layer.get("raw_table", ""))
        if not frozen_p.is_file():
            raise FileNotFoundError(
                f"hERG scope=pool 需要包内冻结的 ADMET 全量表：{frozen_p} 缺失 —— "
                f"该文件属输入层（03_gates/admet/raw/，非本链产物），不得随产物目录一并删除；"
                f"请从输入副本恢复后重跑（本链不重算 ADMET 全量表）。")
        adt = pd.read_csv(frozen_p, low_memory=False)
        if not {"cid", "ADMET_hERG"}.issubset(adt.columns):
            raise RuntimeError("冻结 ADMET 表缺少 cid / ADMET_hERG 列")
        m = m.merge(adt[["cid", "ADMET_hERG"]].rename(columns={"ADMET_hERG": "herg_prob"}),
                    on="cid", how="left")
        if "herg_prob" not in m.columns:      # 冻结表缺列时显式失败，不静默降级
            raise RuntimeError("冻结 ADMET 表合并后缺少 herg_prob 列")
        m["herg_gate_sources"] = "frozen_admet_table"
        if log:
            log.info(f"hERG 门作用域 = 全池（依据包内冻结表 {frozen_p.name}）：覆盖 "
                     f"{int(m.herg_prob.notna().sum())}/{len(m)}；阈值 <{thr}")
    elif herg_path.is_file():
        hp = pd.read_csv(herg_path, low_memory=False)
        col = "hERG" if "hERG" in hp.columns else [c for c in hp.columns if "herg" in c.lower()][0]
        m = m.merge(hp[["cid", col]].rename(columns={col: "herg_prob"}), on="cid", how="left")
        m["herg_gate_sources"] = "external_table"
    else:
        m["herg_prob"] = pd.NA      # pending：不删分子
        m["herg_gate_sources"] = "absent"
        if log:
            log.info(f"hERG 表已载入（作用域 union4）：{herg_path.name}，覆盖 "
                     f"{int(m.herg_prob.notna().sum())}/{len(m)}；阈值 <{thr}")

    # 外部落门表：无论主口径如何都并入，用于口径对账（逐值比对）
    if herg_path.is_file():
        hp = pd.read_csv(herg_path, low_memory=False)
        col = "hERG" if "hERG" in hp.columns else [c for c in hp.columns if "herg" in c.lower()][0]
        m = m.merge(hp[["cid", col]].rename(columns={col: "herg_prob_external"}),
                    on="cid", how="left")
        ov = m.dropna(subset=["herg_prob", "herg_prob_external"])
        if log and len(ov):
            log.info(f"hERG 落门表与外部表对账：交集 {len(ov)}，max|Δ| = "
                     f"{float((ov.herg_prob - ov.herg_prob_external).abs().max()):.2e}")

    # 三态门状态：pass / fail / not_evaluated（**不得**把缺值写成 fail）
    _has = m.herg_prob.notna()
    m["herg_gate_state"] = [
        "not_evaluated" if not h else ("pass" if float(v) < thr else "fail")
        for h, v in zip(_has, m.herg_prob.fillna(0.0))]
    m["herg_gate"] = pd.array(
        [pd.NA if not h else int(float(v) < thr)
         for h, v in zip(_has, m.herg_prob.fillna(0.0))], dtype="Int64")
    if _has.any():
        m["herg_status"] = "applied"
        if scope_mode != "pool" and log:
            log.info(f"hERG 缺值 {int((~_has).sum())} 个 → herg_gate_state=not_evaluated"
                     f"（不当作不合格）")
    else:
        m["herg_gate"] = pd.array([1] * len(m), dtype="Int64")   # pending：不删分子
        m["herg_status"] = "pending"
        if log:
            log.info(f"hERG 表缺失（{herg_path.name}）→ herg_gate=pending（不删分子，报告标注偏差）")
    return m
