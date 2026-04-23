## function
import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt

def run_species_gsea(deg_merge, species, top_n=5, seed=6, gene_sets = "MSigDB_Hallmark_2020"):
    """
    最终版 GSEA 函数
    自动匹配：gene_{species} + logFC_{species}
    数据库：MSigDB_Hallmark_2020
    输出：显著通路(FDR<0.05) + GSEA图 + 气泡图
    """
    # ===================== 自动匹配列名 =====================
    gene_col = f"gene_{species}"
    logfc_col = f"logFC_{species}"
    
    # 构建排序基因集
    rnk = deg_merge[[gene_col, logfc_col]].dropna()
    rnk = rnk.sort_values(logfc_col, ascending=False)

    # ===================== 运行 GSEA =====================
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        seed=seed
    )
    
    # ===================== 结果处理 =====================
    res = pre_res.res2d

    # ✅ 关键：按 FDR 显著性从小到大排序
    res_sorted = res.sort_values("FDR q-val", ascending=True)

    # ✅ 筛选：只保留显著通路 FDR < 0.05
    res_sig = res_sorted[res_sorted["FDR q-val"] < 0.05].copy()

    # ===================== 绘图并保存 =====================
    if len(res_sig) > 0:
        top_terms = res_sig.Term[:top_n].tolist()
        print(f"✅ {species} 找到 {len(res_sig)} 条显著通路，绘制并保存图片")

        # 1. GSEA 曲线图
        axes = pre_res.plot(terms=top_terms)
        plt.savefig(f"{species}_GSEA_plot.png", dpi=300, bbox_inches="tight")
        plt.close()

        # 2. 气泡图
        gp.dotplot(
            res_sig,
            column="FDR q-val",
            title=f"{species} GSEA (FDR < 0.05)",
            cmap="RdBu_r",
            top_n=15
        )
        plt.savefig(f"{species}_bubble_plot.png", dpi=dpi, bbox_inches="tight")
        plt.close()

    else:
        print(f"⚠️ {species} 无显著通路（FDR < 0.05）")
    return pre_res, res_sig, res_sorted