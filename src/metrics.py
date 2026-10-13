# -*- coding: utf-8 -*-
"""
虚拟筛选效能评价指标计算：ROC-AUC, PR-AUC, BEDROC, Enrichment Factor (EF)
"""
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

def calc_bedroc(y_true, scores, alpha=20.0):
    """
    计算玻尔兹曼增强早期区分度 (BEDROC, alpha=20.0 对应关注前 8% 排名)
    """
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(scores, dtype=float)
    N = len(y)
    n = int(np.sum(y))
    if n == 0 or n == N:
        return 0.0
    
    # 按照分数从大到小排序
    order = np.argsort(-s)
    y_sorted = y[order]
    
    # 获得阳性化合物的 1-based 排名
    ranks = np.where(y_sorted == 1)[0] + 1
    
    # 计算 BEDROC
    s_sum = np.sum(np.exp(-alpha * ranks / N))
    ra = n / N
    rand_sum = (ra * (1.0 - np.exp(-alpha))) / (np.exp(alpha / N) - 1.0)
    max_sum = (1.0 - np.exp(-alpha * ra)) / (np.exp(alpha / N) - 1.0)
    
    bedroc = (s_sum - rand_sum) / (max_sum - rand_sum + 1e-12)
    return float(np.clip(bedroc, 0.0, 1.0))

def calc_enrichment_factor(y_true, scores, fraction=0.01):
    """
    计算富集因子 EF@x%
    """
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(scores, dtype=float)
    N = len(y)
    n = int(np.sum(y))
    if n == 0:
        return 0.0
    k = max(1, int(round(fraction * N)))
    order = np.argsort(-s)
    top_hits = np.sum(y[order[:k]])
    ef = (top_hits / k) / (n / N)
    return float(ef)
