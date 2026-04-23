import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os
os.makedirs("figs", exist_ok=True)  # 自动创建文件夹

# ====================== 优化版火山图 ======================
def plot_volcano_specific(deg_merge, padj_thresh=0.05, logfc_thresh=1, top_n=5):
    df = deg_merge.copy().dropna(subset=["logFC_human", "padj_human", "logFC_mouse", "padj_mouse"])

    # 分开绘制：人类 + 小鼠（不再合并，逻辑更清晰）
    human_df = df[["gene_human", "logFC_human", "padj_human", "deg_class"]].copy()
    human_df.columns = ["gene", "logFC", "padj", "deg_class"]
    human_df["species"] = "human"

    mouse_df = df[["Mouse gene name", "logFC_mouse", "padj_mouse", "deg_class"]].copy()
    mouse_df.columns = ["gene", "logFC", "padj", "deg_class"]
    mouse_df["species"] = "mouse"

    plot_df = pd.concat([human_df, mouse_df], ignore_index=True)
    plot_df["-log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))

    # 配色（更柔和、更专业）
    COLORS = {
        "human_specific": "#FF5252",    # 人特异：亮红
        "mouse_specific": "#42A5F5",    # 鼠特异：亮蓝
        "opposite": "#AB47BC",          # 相反趋势：紫色
        "other": "#E0E0E0"              # 其他：浅灰
    }

    plt.figure(figsize=(4, 4))

    # 1. 背景点（全部不显著/保守）
    bg = plot_df[~plot_df["deg_class"].isin(["human_specific", "mouse_specific", "opposite"])]
    plt.scatter(bg["logFC"], bg["-log10_padj"],
                c=COLORS["other"], s=1.5, alpha=0.5, zorder=1, label="Other")

    # 2. 人类特异
    hum = plot_df[plot_df["deg_class"] == "human_specific"]
    plt.scatter(hum["logFC"], hum["-log10_padj"],
                c=COLORS["human_specific"], s=2.5, alpha=0.9, zorder=3, label="Human-specific")

    # 3. 小鼠特异
    mou = plot_df[plot_df["deg_class"] == "mouse_specific"]
    plt.scatter(mou["logFC"], mou["-log10_padj"],
                c=COLORS["mouse_specific"], s=2.5, alpha=0.9, zorder=3, label="Mouse-specific")

    # 4. 相反趋势
    opp = plot_df[plot_df["deg_class"] == "opposite"]
    plt.scatter(opp["logFC"], opp["-log10_padj"],
                c=COLORS["opposite"], s=3, alpha=1, zorder=4, label="Opposite")

    # ====================== 智能标注（按你之前的排序规则） ======================
    texts = []

    # 人特异：按人类 logFC 绝对值排序
    if not hum.empty:
        top_hum = hum.assign(abs_lfc=hum["logFC"].abs()).sort_values("abs_lfc", ascending=False).head(top_n)
        for _, r in top_hum.iterrows():
            texts.append(plt.text(r["logFC"], r["-log10_padj"], str(r["gene"]),
                                  fontsize=7, weight='bold', 
                                  bbox=dict(boxstyle="round,pad=0.2", edgecolor="k", facecolor="none", linewidth=0.4)
                                 )
                        )

    # 鼠特异：按小鼠 logFC 绝对值排序
    if not mou.empty:
        top_mou = mou.assign(abs_lfc=mou["logFC"].abs()).sort_values("abs_lfc", ascending=False).head(top_n)
        for _, r in top_mou.iterrows():
            texts.append(plt.text(r["logFC"], r["-log10_padj"], str(r["gene"]),
                                  fontsize=7, weight='bold',
                                 bbox=dict(boxstyle="round,pad=0.2", edgecolor="k", facecolor="none", linewidth=0.4)
                                 )
                        )

    # 相反：按人和小鼠 logFC 平均值排序
    if not opp.empty:
        # 匹配你之前的定制排序：平均值
        opp = opp.copy()
        opp["avg_abs"] = (abs(df["logFC_human"]) + abs(df["logFC_mouse"])).reindex(opp.index) / 2
        top_opp = opp.sort_values("avg_abs", ascending=False).head(top_n)
        for _, r in top_opp.iterrows():
            texts.append(plt.text(r["logFC"], r["-log10_padj"], str(r["gene"]),
                                  fontsize=7, weight='bold',
                                  bbox=dict(boxstyle="round,pad=0.2", edgecolor="k", facecolor="none", linewidth=0.4)
                                 )
                        )

    # 自动调整文字，不重叠
    adjust_text(texts,
                arrowprops=dict(arrowstyle="-", color='black', lw=0.4),
                expand_points=(1.5, 1.8),
                force_points=0.2)

    # 阈值线
    plt.axvline(-logfc_thresh, c="gray", linestyle="--", lw=1, zorder=0)
    plt.axvline(logfc_thresh, c="gray", linestyle="--", lw=1, zorder=0)
    plt.axhline(-np.log10(padj_thresh), c="gray", linestyle="--", lw=1, zorder=0)

    # 样式美化
    plt.xlabel("log2FC", fontsize=11)
    plt.ylabel("-log10(adjusted p-value)", fontsize=11)
    plt.legend(frameon=False, fontsize=7)
    plt.gca().spines[['top', 'right']].set_visible(False)
    plt.tight_layout()

    # 保存
    plt.savefig("figs/volcano_specific_opposite.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()