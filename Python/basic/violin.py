from itertools import combinations
from statannotations.Annotator import Annotator
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
import itertools
import numpy as np

def adjust_group_color(
    color,
    group_idx,
    n_groups,
    lightness_strength: float = 0.25,
    saturation_strength: float = 0.35,
    max_lightness: float = 0.75,
):
    """
    Adjust group color while keeping
    mutation hue stable.

    Parameters
    ----------
    color :
        Base RGB tuple.

    group_idx :
        Group index.

    n_groups :
        Total group count.

    lightness_strength :
        How much brighter later groups become.

        Larger value:
            stronger lightness contrast.

    saturation_strength :
        How much saturation decreases.

        Larger value:
            stronger pale effect.

    max_lightness :
        Upper limit of brightness.

        Prevents colors becoming nearly white.
    """

    r, g, b = color

    h, l, s = colorsys.rgb_to_hls(
        r,
        g,
        b,
    )

    if n_groups == 1:
        return color

    ratio = group_idx / (n_groups - 1)

    # increase lightness
    l = min(
        max_lightness,
        l + ratio * lightness_strength,
    )

    # reduce saturation
    s = max(
        0.05,
        s - ratio * saturation_strength,
    )

    return colorsys.hls_to_rgb(
        h,
        l,
        s,
    )


def violin_plot_advanced(
    df: pd.DataFrame,
    key: str,
    values: List[str],
    outfile: str,
    xlabel: str = "SV Type",
    ylabel: str = "Delta Frequency",
    order: Optional[List[str]] = None,
    sort_keys: bool = True,
    group_names: Optional[List[str]] = None,
    group_order: Optional[List[str]] = None,
    base_palette: str = "Set2",
    figsize: tuple = (16, 8),
):
    """
    Draw violin + strip plot.

    Features
    --------
    - Different mutations:
        different base colors

    - Different groups:
        different saturation/lightness
        within same mutation

    - Black solid strip points

    - Publication-style aesthetics
    """

    # ==========================================================
    # matplotlib settings
    # ==========================================================

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    sns.set_style("whitegrid")

    # ==========================================================
    # validation
    # ==========================================================

    if key not in df.columns:
        raise ValueError(
            f"x-axis column not found: {key}"
        )

    if not values:
        raise ValueError(
            "`values` cannot be empty"
        )

    missing_cols = [
        v for v in values
        if v not in df.columns
    ]

    if missing_cols:
        raise ValueError(
            f"y-axis columns not found: {missing_cols}"
        )

    if (
        group_names is not None
        and len(group_names) != len(values)
    ):
        raise ValueError(
            "`group_names` length must equal "
            "`values` length"
        )

    # ==========================================================
    # build long dataframe
    # ==========================================================

    long_parts = []

    for i, v_col in enumerate(values):

        group_name = (
            group_names[i]
            if group_names is not None
            else v_col
        )

        sub_df = (
            df[[key, v_col]]
            .dropna(subset=[key, v_col])
            .copy()
        )

        sub_df = sub_df.rename(
            columns={v_col: "_y"}
        )

        sub_df["group"] = group_name

        long_parts.append(sub_df)

    if not long_parts:
        raise ValueError(
            "No valid plotting data"
        )

    df_plot = pd.concat(
        long_parts,
        ignore_index=True,
    )

    # ==========================================================
    # x-axis order
    # ==========================================================

    if order is None:

        observed_order = (
            df_plot[key]
            .dropna()
            .unique()
            .tolist()
        )

        if sort_keys:
            order = sorted(observed_order)
        else:
            order = observed_order

    # ==========================================================
    # group order
    # ==========================================================

    if group_order is None:

        group_order = (
            df_plot["group"]
            .dropna()
            .unique()
            .tolist()
        )

    # ==========================================================
    # base colors
    # ==========================================================

    base_colors = sns.color_palette(
        base_palette,
        n_colors=len(order),
    )

    mutation_color_map = {
        mutation: color
        for mutation, color in zip(
            order,
            base_colors,
        )
    }

    # ==========================================================
    # color adjustment
    # ==========================================================

    

    # ==========================================================
    # plotting
    # ==========================================================

    fig, ax = plt.subplots(
        figsize=figsize
    )

    # ----------------------------------------------------------
    # violin plot
    # ----------------------------------------------------------

    sns.violinplot(
        data=df_plot,
        x=key,
        y="_y",
        hue="group",
        order=order,
        hue_order=group_order,
        dodge=True,
        inner="quartile",
        linewidth=1.2,
        saturation=1,
        bw_adjust=1.2,
        density_norm="area",
        width=0.95,
        ax=ax,
    )

    # ==========================================================
    # recolor violins
    # ==========================================================

    violin_bodies = [
        c for c in ax.collections
        if isinstance(
            c,
            mcollections.PolyCollection
        )
    ]

    violin_idx = 0

    for mutation in order:

        base_color = mutation_color_map[
            mutation
        ]

        for group_idx, group in enumerate(
            group_order
        ):

            if violin_idx >= len(
                violin_bodies
            ):
                break

            color = adjust_group_color(
                base_color,
                group_idx,
                len(group_order),
            )

            violin = violin_bodies[
                violin_idx
            ]

            violin.set_facecolor(color)

            violin.set_edgecolor(
                "black"
            )

            violin.set_alpha(0.95)

            violin.set_linewidth(1)

            violin_idx += 1

    # ==========================================================
    # black solid strip points
    # ==========================================================

    sns.stripplot(
        data=df_plot,
        x=key,
        y="_y",
        hue="group",
        order=order,
        hue_order=group_order,
        dodge=True,
        jitter=0.12,
        color="black",
        size=2.5,
        alpha=0.9,
        linewidth=0,
        ax=ax,
    )

    # ==========================================================
    # remove duplicated legends
    # ==========================================================

    if ax.get_legend() is not None:
        ax.get_legend().remove()

    # ==========================================================
    # explanatory legend
    # ==========================================================

    example_mutation = order[0]

    example_base_color = mutation_color_map[
        example_mutation
    ]

    legend_handles = []

    n_groups = len(group_order)

    if n_groups == 2:

        shade_words = [
            "Dark",
            "Light",
        ]

    else:

        shade_words = [
            f"Shade {i + 1}"
            for i in range(n_groups)
        ]

    for group_idx, group in enumerate(
        group_order
    ):

        legend_color = adjust_group_color(
            example_base_color,
            group_idx,
            n_groups,
        )

        label = (
            f"{shade_words[group_idx]} "
            f"{group}"
        )

        legend_handles.append(
            mpatches.Patch(
                facecolor=legend_color,
                edgecolor="black",
                label=label,
            )
        )

    ax.legend(
        handles=legend_handles,
        title="Color meaning",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=11,
        title_fontsize=12,
    )

    # ==========================================================
    # labels
    # ==========================================================

    ax.set_xlabel(
        xlabel,
        fontsize=14,
    )

    ax.set_ylabel(
        ylabel,
        fontsize=14,
    )

    # ==========================================================
    # ticks
    # ==========================================================

    ax.tick_params(
        axis="x",
        rotation=45,
        labelsize=11,
    )

    ax.tick_params(
        axis="y",
        labelsize=11,
    )

    # ==========================================================
    # style
    # ==========================================================

    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    # ==========================================================
    # layout
    # ==========================================================

    fig.tight_layout()

    # ==========================================================
    # save
    # ==========================================================

    fig.savefig(
        outfile,
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

def violin_plot(
    df: pd.DataFrame,
    key: str,
    value: str,
    outfile:str,
    xlabel: str = "SV Type",
    ylabel: str = "Delta Frequency (Freq - ddPCR_AF)",
):
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.violinplot(
        x=key,
        y=value,
        data=df,
        inner='box',
        linewidth=1.5,
        hue=key,
        palette='Set2',
        ax=ax
    )
    sns.stripplot(
        x=key,
        y=value,
        data=df,
        color='black',
        size=3,
        alpha=0.5,
        ax=ax
    )
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis='x', labelrotation=45)
    fig.savefig(outfile, dpi=300, bbox_inches='tight')

def violin(
        df:pd.DataFrame,
        outplot:str,
        key:str,
        value:str
):
    unique_cell_types = df[key].unique().tolist()
    box_pairs = list(combinations(unique_cell_types, 2))
    # combinations(unique_cell_types, 2) 会生成所有不重复的 (组A, 组B) 组合

    plt.figure(figsize=(10, 8))
    ax = sns.violinplot(
        x=key,
        y=value,
        data=df,
        palette='Set2',
        inner='box',
        linewidth=1.5
    )
    annotator = Annotator(
        ax,
        box_pairs,
        data=df,
        x=key,
        y=value,
        order=unique_cell_types # 确保 order 与图中的顺序一致
    )
    # 注意：进行多重比较时，通常需要进行多重检验校正 (Multiple Testing Correction)。
    annotator.configure(
        test='Mann-Whitney',
        text_format='star',
        loc='inside',
        # 多重检验校正：例如使用 Benjamini/Hochberg FDR
        comparisons_correction='fdr_bh', # <--- 关键参数！
        pvalue_thresholds=[
            [0.05, "*"],
            [0.01, "**"],
            [0.001, "***"],
            [0.0001, "****"]
        ]
    )

    annotator.apply_and_annotate()
    ax.set_xlabel('') 
    plt.ylabel('GSI')
    plt.tight_layout()
    
    plt.savefig(outplot,dpi=300)




def violin_by_method_celltype_step(df_dict: dict, outplot: str):
    all_dfs = []
    for method, df in df_dict.items():
        tmp = df.copy()
        tmp['Method'] = method
        all_dfs.append(tmp)
    combined_df = pd.concat(all_dfs, axis=0)

    methods = combined_df['Method'].unique().tolist()
    cell_types = combined_df['Cell_type'].unique().tolist()
    palette = sns.color_palette("Set2", n_colors=len(cell_types))
    cell_color_map = dict(zip(cell_types, palette))

    plt.figure(figsize=(12, 8))
    ax = sns.violinplot(
        x='Method',
        y='LV1',
        hue='Cell_type',
        data=combined_df,
        palette=cell_color_map,
        inner='box',
        linewidth=1.5,
        dodge=True
    )

    # 计算每个方法+cell_type的小提琴中心点
    violin_positions = {}
    for method_idx, method in enumerate(methods):
        cells_in_method = combined_df[combined_df['Method']==method]['Cell_type'].unique()
        n_cell_types = len(cells_in_method)
        offsets = np.linspace(-0.2, 0.2, n_cell_types)
        for i, ct in enumerate(cells_in_method):
            violin_positions[(method, ct)] = method_idx + offsets[i]

    # 遍历每个方法做两两检验
    for method in methods:
        subset = combined_df[combined_df['Method'] == method]
        cells_in_method = subset['Cell_type'].unique().tolist()
        if len(cells_in_method) < 2:
            continue
        pairs = list(itertools.combinations(cells_in_method, 2))
        step = subset['LV1'].max() * 0.05

        for i, (c1, c2) in enumerate(pairs):
            data1 = subset[subset['Cell_type']==c1]['LV1']
            data2 = subset[subset['Cell_type']==c2]['LV1']
            stat, p = mannwhitneyu(data1, data2, alternative='two-sided')

            x1 = violin_positions[(method, c1)]
            x2 = violin_positions[(method, c2)]
            y_base = max(data1.max(), data2.max())
            y = y_base + (i+1)*step
            y_top = y + step*0.3  # 梯形顶部高度

            # 画梯形：竖直-水平-竖直
            ax.plot([x1, x1], [y, y_top], color='black', lw=2)
            ax.plot([x1, x2], [y_top, y_top], color='black', lw=2)
            ax.plot([x2, x2], [y_top, y], color='black', lw=2)

            # 显示星号
            ax.text((x1+x2)/2, y_top + step*0.05,
                    "****" if p<0.0001 else
                    "***" if p<0.001 else
                    "**" if p<0.01 else
                    "*" if p<0.05 else "ns",
                    ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Method')
    ax.set_ylabel('LV1 score')
    plt.tight_layout()
    plt.savefig(outplot, dpi=300)
    plt.close()
    print(f"图已保存到 {outplot}")

def combine_run():
    Latent = {
        "pca": "/home/luosg/Data/genomeStability/output/result/latent/fa_factors.csv",
        "nmf": "/home/luosg/Data/genomeStability/output/result/latent/nmf_factors.csv",
        "factor": "/home/luosg/Data/genomeStability/output/result/latent/pca_factors.csv"
    }
    # Latent = {
    #     "pca": "/home/luosg/Data/genomeStability/output/result1/latent/pca_factors.csv",
    #     "nmf": "/home/luosg/Data/genomeStability/output/result1/latent/nmf_factors.csv",
    #     "factor": "/home/luosg/Data/genomeStability/output/result1/latent/fa_factors.csv" 
    # }
    outdir="/home/luosg/Data/genomeStability/output/result/latent"
    # outdir = "/home/luosg/Data/genomeStability/output/result1/latent"
    df_dict = {}
    for method, infile in Latent.items():
        df = pd.read_csv(infile, index_col=0)
        df_annotation = pd.read_csv("/home/luosg/Data/genomeStability/data/target_fq.tsv", sep="\t")
        df_annotation.drop_duplicates(subset=["Sample_id"], inplace=True)
        Sample_to_Status_map = df_annotation.set_index('Sample_id')['Status']
        new_index = df.index.map(Sample_to_Status_map)
        df.index = new_index
        df_scaled = (df - df.min()) / (df.max() - df.min())
        df_scaled = df_scaled.reset_index(names="Cell_type")
        df_dict[method] = df_scaled

    outplot = f"{outdir}/all_methods_violin.png"
    violin_by_method_celltype_step(df_dict, outplot)
if __name__ == "__main__":
    # single_run()
    combine_run()
