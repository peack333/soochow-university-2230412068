import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

def plot_volcano(
    df,
    suffix="human",    # 只需要改这里！！！
    lfc_cut=1,
    diff_cut=0.2,
    padj_cut=0.05,
    figsize=(6, 5)
):
    # 自动拼接列名（不用手动写了）
    gene_col = f"gene_{suffix}"
    logFC_col = f"logFC_{suffix}"
    padj_col = f"padj_{suffix}"
    pct1_col = f"pct.1_{suffix}"
    pct2_col = f"pct.2_{suffix}"
    
    df = df.copy().dropna()
    df["Difference"] = df[pct1_col] - df[pct2_col]

    # 阈值分组
    def threshold(row):
        if row[logFC_col] > lfc_cut and row["Difference"] > diff_cut and row[padj_col] < padj_cut:
            return "Up"
        elif row[logFC_col] < -lfc_cut and row["Difference"] < -diff_cut and row[padj_col] < padj_cut:
            return "Down"
        else:
            return "NS"

    df["threshold"] = df.apply(threshold, axis=1)

    # 绘图
    plt.figure(figsize=figsize)
    color_map = {"Down": "#1f77b4", "NS": "#bdbdbd", "Up": "#d62728"}

    for group in ["Down", "NS", "Up"]:
        sub = df[df["threshold"] == group]
        plt.scatter(sub["Difference"], sub[logFC_col], c=color_map[group], 
                    s=2, alpha=0.8, label=group)

    # 上调、下调各取10个
    up = df[df["threshold"] == "Up"].sort_values(logFC_col, ascending=False).head(10)
    down = df[df["threshold"] == "Down"].sort_values(logFC_col, ascending=True).head(10)

    # 标签（无灰底）
    texts = []
    for _, row in pd.concat([up, down]).iterrows():
        texts.append(plt.text(
            row["Difference"], row[logFC_col],
            row[gene_col],
            fontsize=8,
            color="black",
            bbox=dict(boxstyle="round,pad=0.2", edgecolor="black", facecolor="white", alpha=0.9, linewidth=0.4)
        ))
        
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.3))

    # 参考线
    plt.axvline(0, linestyle='--', c='k', lw=1)
    plt.axhline(0, linestyle='--', c='k', lw=1)

    # 样式
    plt.gca().spines[['top','right']].set_visible(False)
    plt.xlabel("Difference (pct.1 - pct.2)")
    plt.ylabel("log2FC")
    plt.legend(frameon=False)  # 图例无灰底无边框
    plt.tight_layout()

    plt.savefig(
    f"figs/T_cell_volcano_{suffix}.png",
    bbox_inches="tight",
    facecolor="white"
)
    
    plt.show()

    up_genes = df[df["threshold"] == "Up"].sort_values(logFC_col, ascending=False)[gene_col].tolist()
    down_genes = df[df["threshold"] == "Down"].sort_values(logFC_col, ascending=True)[gene_col].tolist()
    all_significant_genes = up_genes + down_genes

    return {
        "up_genes": up_genes,
        "down_genes": down_genes,
        "all_significant_genes": all_significant_genes
    }