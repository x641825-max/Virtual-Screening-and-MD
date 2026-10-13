#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
一键模型重训与验证入口 (train.py)
基于人源 A2AR 活性数据与实验非活性数据重新训练随机森林模型并导出纯 numpy 推理权重
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)

if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

from src.features import batch_smiles_to_ecfp4
from src.metrics import calc_bedroc, calc_enrichment_factor

PROJECT_SEED = 20260905
np.random.seed(PROJECT_SEED)

def export_portable_npz(rf, npz_path):
    """将 sklearn 随机森林导出为纯 numpy 便携包 (与 PortableRF 严格兼容)"""
    n_classes = int(rf.n_classes_)
    n_features = int(rf.n_features_in_)

    cl, cr, feat, thr, val, offs = [], [], [], [], [], [0]
    for est in rf.estimators_:
        t = est.tree_
        cl.append(t.children_left.astype(np.int32))
        cr.append(t.children_right.astype(np.int32))
        feat.append(t.feature.astype(np.int32))
        thr.append(t.threshold.astype(np.float64))
        # 分类树 value 形状 [n_nodes, 1, n_classes]，取第 0 个 output
        v = t.value[:, 0, :].astype(np.float64)
        s = v.sum(axis=1, keepdims=True)
        v = np.divide(v, s, out=np.full_like(v, 1.0 / n_classes), where=s > 0)
        val.append(v)
        offs.append(offs[-1] + len(t.children_left))

    np.savez_compressed(
        npz_path,
        children_left=np.concatenate(cl),
        children_right=np.concatenate(cr),
        feature=np.concatenate(feat),
        threshold=np.concatenate(thr),
        value=np.concatenate(val),
        offsets=np.array(offs, dtype=np.int64),
        n_classes=np.array([n_classes]),
        n_features=np.array([n_features]),
        n_trees=np.array([len(rf.estimators_)]),
    )
    print(f"[+] 纯 numpy 模型导出成功: {npz_path} ({os.path.getsize(npz_path)/1024/1024:.2f} MB)")

def main():
    train_path = os.path.join(BASE_DIR, "data", "dataset_train.csv")
    ext_path = os.path.join(BASE_DIR, "data", "dataset_external.csv")
    output_npz = os.path.join(BASE_DIR, "models", "a2a_rf_trained_portable.npz")
    
    print("==================================================================")
    print("第一届全球大学生生命科学挑战赛 · A2AR 初筛模型一键训练")
    print(f"训练集: {train_path}")
    print(f"独立外部集: {ext_path}")
    print("==================================================================")
    
    df_train = pd.read_csv(train_path)
    X_train, valid_tr = batch_smiles_to_ecfp4(df_train["smiles"].tolist())
    y_train = df_train["label"].values[valid_tr]
    X_train = X_train[valid_tr]
    
    print(f"[+] 提取训练集 ECFP4 特征完成: {len(X_train)} 个样本 (活性 {y_train.sum()} 个)")
    
    rf = RandomForestClassifier(
        n_estimators=400,
        min_samples_leaf=2,
        class_weight="balanced_subsample",
        n_jobs=-1,
        random_state=PROJECT_SEED
    )
    
    print("[*] 正在训练 400 棵决策树随机森林模型...")
    rf.fit(X_train, y_train)
    print("[+] 模型训练完成！")
    
    if os.path.exists(ext_path):
        df_ext = pd.read_csv(ext_path)
        X_ext, valid_ext = batch_smiles_to_ecfp4(df_ext["smiles"].tolist())
        y_ext = df_ext["label"].values[valid_ext]
        X_ext = X_ext[valid_ext]
        
        preds_ext = rf.predict_proba(X_ext)[:, 1]
        auc = roc_auc_score(y_ext, preds_ext)
        pr_auc = average_precision_score(y_ext, preds_ext)
        bedroc = calc_bedroc(y_ext, preds_ext, alpha=20.0)
        ef1 = calc_enrichment_factor(y_ext, preds_ext, fraction=0.01)
        
        print(f"\n================ 独立外部测试集泛化验证指标 ================")
        print(f"ROC-AUC:       {auc:.4f}")
        print(f"PR-AUC:        {pr_auc:.4f}")
        print(f"BEDROC(α=20):  {bedroc:.4f}")
        print(f"EF@1%:         {ef1:.2f}")
        print(f"============================================================")

    export_portable_npz(rf, output_npz)

if __name__ == "__main__":
    main()
