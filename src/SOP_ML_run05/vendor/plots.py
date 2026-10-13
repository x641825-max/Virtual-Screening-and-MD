"""绘图统一风格：CJK 字体、无衬线、统一色板、PNG+SVG 双出。"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

CJK_CANDIDATES = ["Noto Sans CJK SC", "Noto Sans CJK JP", "Noto Sans SC",
                  "Source Han Sans SC", "WenQuanYi Zen Hei", "DejaVu Sans"]
BLUE, GREY, RED, GREEN, AMBER = "#2f6fb3", "#8a939e", "#c0504d", "#4f8a63", "#d08c1e"


def setup_cjk() -> str:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in CJK_CANDIDATES:
        if name in available:
            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams["font.sans-serif"] = [name] + CJK_CANDIDATES
            plt.rcParams["axes.unicode_minus"] = False
            return name
    plt.rcParams["axes.unicode_minus"] = False
    return "DEJAVU(无 CJK 字体)"


def save(fig, path, dpi: int = 200) -> list[str]:
    """写 PNG + SVG。SVG 需要显式去元数据才可复现：

    matplotlib 默认会把生成日期与随机 hashsalt 写进 SVG，导致同一份数据两次运行
    SVG 字节不同（确定性 G 臂因此失败）。此处固定 hashsalt 并关闭 Date 元数据。
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams["svg.hashsalt"] = "SOP-05"
    out = []
    fig.savefig(p.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    out.append(str(p.with_suffix(".png")))
    fig.savefig(p.with_suffix(".svg"), bbox_inches="tight",
                metadata={"Date": None, "Creator": "SOP-05/matplotlib"})
    out.append(str(p.with_suffix(".svg")))
    plt.close(fig)
    return out


def funnel_barh(steps: dict, title: str, path, xlabel: str = "分子数"):
    labels = list(steps.keys())
    vals = [steps[k] for k in labels]
    fig, ax = plt.subplots(figsize=(8.2, 0.55 * len(labels) + 1.5))
    y = range(len(labels))[::-1]
    colors = [BLUE] * len(labels)
    colors[-1] = GREEN
    ax.barh(list(y), vals, color=colors, height=0.62)
    for yy, v in zip(y, vals):
        ax.text(v + max(vals) * 0.012, yy, f"{v:,}", va="center", fontsize=9)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.set_xlim(0, max(vals) * 1.12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    return save(fig, path)


def grouped_barh(cats, series: dict, title, path, xlabel="分子数"):
    """series: {名称: [每类数值]}"""
    import numpy as np
    n = len(cats)
    k = len(series)
    fig, ax = plt.subplots(figsize=(8.4, 0.5 * n + 1.8))
    idx = np.arange(n)
    h = 0.8 / k
    palette = [BLUE, AMBER, GREEN, RED, GREY]
    for j, (name, vals) in enumerate(series.items()):
        ax.barh(idx + (j - (k - 1) / 2) * h, vals, height=h * 0.92,
                color=palette[j % len(palette)], label=name)
    ax.set_yticks(idx)
    ax.set_yticklabels(cats, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    return save(fig, path)


def heatmap(mat, row_labels, col_labels, title, path, cbar_label="overlap", fmt="{:.2f}"):
    import numpy as np
    m = np.asarray(mat, dtype=float)
    fig, ax = plt.subplots(figsize=(1.0 * len(col_labels) + 3.2, 0.55 * len(row_labels) + 2.0))
    im = ax.imshow(m, cmap="Blues", vmin=0, vmax=max(1.0, float(np.nanmax(m))))
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=30, ha="right", fontsize=9)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, fmt.format(m[i, j]), ha="center", va="center", fontsize=8,
                    color="white" if m[i, j] > 0.6 * np.nanmax(m) else "#22303f")
    ax.set_title(title, fontsize=11)
    fig.colorbar(im, ax=ax, label=cbar_label, shrink=0.8)
    fig.tight_layout()
    return save(fig, path)


def line_ci(x, ys: dict, title, path, xlabel, ylabel, refs: dict | None = None):
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    palette = [BLUE, AMBER, GREEN, RED]
    for j, (name, y) in enumerate(ys.items()):
        ax.plot(x, y, marker="o", markersize=4, color=palette[j % len(palette)], label=name)
    if refs:
        for name, v in refs.items():
            ax.axhline(v, color=GREY, linestyle="--", linewidth=1)
            ax.text(x[-1], v, f" {name} {v}", fontsize=8, color=GREY, va="bottom")
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    return save(fig, path)


def hist(values, title, path, xlabel, bins=30, color=BLUE):
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    ax.hist(values, bins=bins, color=color, alpha=0.9)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("分子数", fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    return save(fig, path)
