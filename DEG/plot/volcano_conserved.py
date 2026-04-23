import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text

# ====================== 最终完整函数 ======================
def plot_volcano_conserved(deg_merge, padj_thresh=0.05, logfc_thresh=1, top_n=10):
    df = deg_merge.dropna(subset=["logFC_human", "padj_human", "logFC_mouse", "padj_mouse"]).copy()
    
    # 1. 筛选保守基因（完全用你已有的列）
    conserved = df[df["same_direction"] & df["significant_both"]].copy()

    # 2. 上调 TOP 基因（你写的公式 1:1 放进函数）
    conserved_up = conserved[conserved["logFC_human"] > 0].copy()
    conserved_up["mean_logFC"] = (conserved_up["logFC_human"] + conserved_up["logFC_mouse"]) / 2
    top_up = conserved_up.nlargest(top_n, "mean_logFC")

    # 3. 下调 TOP 基因（你写的公式 1:1 放进函数）
    conserved_down = conserved[conserved["logFC_human"] < 0].copy()
    conserved_down["mean_logFC"] = (conserved_down["logFC_human"] + conserved_down["logFC_mouse"]) / 2
    top_down = conserved_down.nsmallest(top_n, "mean_logFC")

    # ====================== 绘图 ======================
    plt.figure(figsize=(4, 4))

    # 背景灰色点
    plt.scatter(df["logFC_human"], -np.log10(df["padj_human"].clip(1e-300)),
                c="#bdbdbd", s=2, alpha=0.6, label="ns")

    # 保守上调
    plt.scatter(conserved_up["logFC_human"], -np.log10(conserved_up["padj_human"].clip(1e-300)),
                c="#2ca02c", s=5, alpha=0.9, label="conserved_up")

    # 保守下调
    plt.scatter(conserved_down["logFC_human"], -np.log10(conserved_down["padj_human"].clip(1e-300)),
                c="#ff7f0e", s=5, alpha=0.9, label="conserved_down")

    # 标注 TOP 基因
    top_labels = pd.concat([top_up, top_down])
    texts = []
    for _, row in top_labels.iterrows():
        texts.append(plt.text(
            row["logFC_human"],
            -np.log10(max(row["padj_human"], 1e-300)),
            str(row["gene_human"]),
            fontsize=7,
            bbox=dict(boxstyle="round,pad=0.2", edgecolor="k", facecolor="none", linewidth=0.4)
        ))

    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.3))

    # 阈值线
    plt.axvline(-logfc_thresh, c="gray", ls="--", lw=1)
    plt.axvline(logfc_thresh, c="gray", ls="--", lw=1)
    plt.axhline(-np.log10(padj_thresh), c="gray", ls="--", lw=1)

    plt.xlabel("log2FC (human)")
    plt.ylabel("-log10(adjusted pvalue)")
    plt.legend(frameon=False, fontsize=8)
    plt.gca().spines[['top','right']].set_visible(False)
    plt.tight_layout()
    plt.savefig("figs/volcano_conserved.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()