import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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




if __name__ == '__main__':
    DEG_file = "/disk5/luosg/scRNAseq/output/result/DEG/Intestine/table/combined_groupCKO_chang_10XSC3_B_cell-combined_groupWT_chang_10XSC3_B_cell_DEG.tsv"
    rmsk_file = "/ChIP_seq_2/Data/index/Mus_musculus/UCSC/mm39/rmsk_mm39.txt.gz"
    TEfamily(DEG_file,rmsk_file,"a")