import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import textwrap
from pathlib import Path
from typing import List, Union


def TEfamily(
    DEG_file: str,
    rmsk_file: str,
    out: str
):
    df_TE = pd.read_csv(DEG_file, sep="\t")
    df_TE.reset_index(names="subfamily", inplace=True)
    df_rmsk = pd.read_csv(rmsk_file, sep="\t", header=None)
    
    df_annotation = df_rmsk[[10, 11]].drop_duplicates(keep="first")
    df_annotation.columns = ["subfamily", "family"]
    
    new_df = pd.merge(df_TE, df_annotation, on="subfamily", how="left")
    
    # 筛选显著
    significant_df = new_df[
        (abs(new_df["logFC"]) > 0.58) & 
        (-np.log10(new_df["PValue"]) > 1.3)
    ].copy()
    
    if significant_df.empty:
        raise ValueError("The filtered DataFrame of significant TEs is empty.")
    
    # 直接分组计数：同时保留方向信息
    significant_df["regulation"] = np.where(significant_df["logFC"] > 0, "Up", "Down")
    result = (
        significant_df.groupby(["family", "regulation"])
        .size()
        .reset_index(name="count")
    )
    
    # 保存 TSV
    result_pivot = result.pivot(index="family", columns="regulation", values="count").fillna(0).astype(int)
    result_pivot.to_csv(f"{out}_TEfamily.tsv", sep="\t", header=True)

    # --- 绘图 ---
    sns.set_style("white")  # 去掉网格线
    fig, ax = plt.subplots(figsize=(12, 6))

    # 分开取调色板：深色给Up，浅色给Down
    up_families = result[result["regulation"] == "Up"]["family"].unique()
    down_families = result[result["regulation"] == "Down"]["family"].unique()
    up_colors = sns.color_palette("dark", len(up_families))
    down_colors = sns.color_palette("pastel", len(down_families))

    # 上调在左侧
    for i, (fam, cnt) in enumerate(result[result["regulation"] == "Up"][["family", "count"]].values):
        ax.bar(i, cnt, color=up_colors[i])
        ax.text(i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    # 下调在右侧
    offset = len(up_families) + 2  # 空2列作为分隔
    for i, (fam, cnt) in enumerate(result[result["regulation"] == "Down"][["family", "count"]].values):
        ax.bar(offset + i, cnt, color=down_colors[i])
        ax.text(offset + i, cnt + 0.1, str(cnt), ha="center", va="bottom", fontsize=9)

    # 设置 x 轴刻度
    xticks = list(range(len(up_families))) + list(range(offset, offset + len(down_families)))
    xlabels = list(up_families) + list(down_families)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right")

    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Top TE Families: Upregulated vs Downregulated", fontsize=16)

    # 重新计算 Y 轴上限，横线位置稍高于最大值

    y_max = max(result["count"])
    # 在上方标注 Up 和 Down
    ax.hlines(y=y_max*1.1, xmin=0, xmax=len(up_families)-1, color="black", linewidth=1.5)
    ax.text((len(up_families)-1)/2, y_max*1.13, "Upregulated", ha="center", va="bottom", fontsize=12)

    ax.hlines(y=y_max*1.1, xmin=offset, xmax=offset+len(down_families)-1, color="black", linewidth=1.5)
    ax.text(offset + (len(down_families)-1)/2, y_max*1.13, "Downregulated", ha="center", va="bottom", fontsize=12)

    y_max = max(result["count"]) * 1.25
    ax.set_ylim(0, y_max)
    plt.tight_layout()
    plt.savefig(f"{out}_TEfamily.png", dpi=300)
    plt.close()

def GSEA_BarPlot_Enhanced(
    gsea_result: pd.DataFrame,
    outplot: str,
    pval_cutoff: float = 0.05,
    sort_by: str = "ES",    # "ES" or "pval"
    top_flag: bool = False,
    top_number: int = 5,
    label_wrap: int = 40,
    cmap_name: str = "viridis"
):
    """
    增强版 GSEA barplot：
    - 中心对称：负 ES 左、正 ES 右
    - 按 ES 或 pval 自动排序
    - 颜色按显著性渐变（-log10(pval)）
    - ggplot 风格
    """

    # ------------------------
    # 1. 过滤和预处理
    # ------------------------
    df = gsea_result[gsea_result["pval"] < pval_cutoff].copy()
    df = df.rename(columns={"score": "ES"}).sort_values("ES")
    if top_flag:
        df = pd.concat([df.head(top_number), df.tail(top_number)])
    df["logP"] = -np.log10(df["pval"].replace(0, 1e-7))

    # 名称换行
    df["pathway_wrapped"] = df["source"].apply(
        lambda s: "\n".join(textwrap.wrap(s, label_wrap))
    )


    # ------------------------
    # 4. 图尺寸
    # ------------------------
    height = max(6, len(df) * 0.45)
    fig, ax = plt.subplots(figsize=(12, height))

    # ------------------------
    # 5. 颜色映射
    # ------------------------
    cmap = sns.color_palette(cmap_name, as_cmap=True)
    norm = plt.Normalize(df["logP"].min(), df["logP"].max())
    df["color"] = df["logP"].apply(lambda x: cmap(norm(x)))

    # ------------------------
    # 6. 中心对称 barplot
    # ------------------------
    ax.barh(
        y=df["pathway_wrapped"],
        width=df["ES"],
        color=df["color"],
        edgecolor="black",
        linewidth=0.5
    )

    # ------------------------
    # 7. 中心线（x=0）
    # ------------------------
    ax.axvline(0, color="black", linewidth=1.2)

    # ------------------------
    # 8. 坐标轴与标签
    # ------------------------
    ax.set_xlabel("Enrichment Score (ES)")
    ax.set_ylabel("Pathway")

    ax.set_title(
        "GSEA Enrichment Barplot\n"
        "Left: Negative ES | Right: Positive ES",
        fontsize=14,
        pad=15
    )

    # ------------------------
    # 9. Colorbar (p-value 显著性)
    # ------------------------
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("-log10(P-value)")

    # ------------------------
    # 10. 美化
    # ------------------------
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    ax.grid(False, axis="y")

    plt.tight_layout()
    plt.savefig(outplot, dpi=300, bbox_inches="tight")
    plt.close()

# -----------------------------
# 函数3：堆叠柱状图（支持 x 轴彩色标签 + 样本分组图例）
# -----------------------------
def plot_mutation_distribution_multi(
        df_counts: pd.DataFrame,
        xlabels: List[str] = None,
        groups: List[str] = None,
        group_colors: dict = None,
        title: str = "Func.refGene Mutation Distribution (Proportion)",
        save_path: Union[str, Path] = None,
        legend_width:float=0.2,
        figsize:tuple = (12, 6)):

    df_prop = df_counts.div(df_counts.sum(axis=0), axis=1)

    fig, ax = plt.subplots(figsize=figsize)

    # ======== 核心：右侧预留空间给图例 ========
    fig.subplots_adjust(right=1 - legend_width)

    # 堆叠柱状图
    df_prop.T.plot(
        kind="bar",
        stacked=True,
        colormap="tab20",
        width=0.8,
        ax=ax,
        legend=False
    )

    n_samples = df_counts.shape[1]

    # -------------------------------
    # X 轴标签设置
    # -------------------------------
    if xlabels is None:
        xlabels = df_counts.columns.tolist()

    if len(xlabels) != n_samples:
        raise ValueError("xlabels 长度必须与样本数量一致")

    xlabel_colors = ["black"] * n_samples

    # -------------------------------
    # 根据分组上色
    # -------------------------------
    if groups is not None:
        if len(groups) != n_samples:
            raise ValueError("groups 长度必须与样本数量一致")

        if group_colors is None:
            unique_groups = list(dict.fromkeys(groups))
            cmap = plt.get_cmap("tab10")
            group_colors = {g: cmap(i) for i, g in enumerate(unique_groups)}

        xlabel_colors = [group_colors[g] for g in groups]

    # 设置 x 轴颜色
    for label, c in zip(ax.get_xticklabels(), xlabel_colors):
        label.set_color(c)

    # -------------------------------
    # 第一个图例：突变类型
    # -------------------------------
    mutation_types = df_prop.index.tolist()
    cmap = plt.get_cmap("tab20")
    mutation_patches = [
        mpatches.Patch(color=cmap(i), label=mutation_types[i])
        for i in range(len(mutation_types))
    ]

    legend1 = ax.legend(
        handles=mutation_patches,
        title="Mutation Type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )
    ax.add_artist(legend1)

    # -------------------------------
    # 第二个图例：样本分组（绘制在右下）
    # -------------------------------
    if groups is not None:
        unique_groups = sorted(set(groups))
        group_patches = [
            mpatches.Patch(color=group_colors[g], label=g)
            for g in unique_groups
        ]

        ax.legend(
            handles=group_patches,
            title="Sample Group",
            bbox_to_anchor=(1.02, 0.25),
            loc="upper left"
        )

    ax.set_xlabel("Sample")
    ax.set_ylabel("Proportion")
    ax.set_title(title)
    ax.set_xticklabels(xlabels, rotation=45, ha='right')

    if save_path:
        plt.savefig(save_path, dpi=300)


if __name__ == '__main__':
    # DEG_file = "/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell_DEG.tsv"
    # rmsk_file = "/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz"
    # TEfamily(DEG_file,rmsk_file,"a")

    gsea_result = pd.read_csv("/disk5/luosg/scRNAseq/output/combine/Intestine/GSEA/table/combined_groupE2_chang_10XSC3_Intestinal_stem_cell-combined_groupCKO_chang_10XSC3_Intestinal_stem_cell_DEG_gene_gsea.tsv",sep="\t")
    GSEA_BarPlot_Enhanced(
        gsea_result=gsea_result,
        outplot="/disk5/luosg/scRNAseq/output/combine/Intestine/GSEA/plotTop/Plasmacytoid_dendritic_cells_E2_CKO_gsea.png",
        sort_by="pval",
        top_flag=True,
        top_number=5
    )
