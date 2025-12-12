import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import matplotlib.gridspec as gridspec

def plot_heatmap(
        heatmap_data: pd.DataFrame,
        outplot:str,
        title: str = "Heatmap of Positions",
        xlabel: str = "Positions",
        ylabel: str = "Samples"
):
    """
    heatmap_data: 2D DataFrame where rows are samples and columns are positions, with values to be represented in the heatmap.
    """
    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, annot=True, cmap='coolwarm', cbar_kws={'label': 'Values'}, linewidths=0.5)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(outplot,dpi=300)
    plt.close()

def plot_heatmap_line():
    ### 数据详情 WT-1,WT-2,OE-1,OE-2, score
    dfnew = pd.read_csv("/disk5/luosg/BETA20250928/output/heatmap/Spz1_heamtmap_final.txt",sep="\t")
    df_plot = dfnew[["WT-1","WT-2","OE-1","OE-2"]]
    df_zscore = df_plot.sub(df_plot.mean(axis=1), axis=0).div(df_plot.std(axis=1), axis=0)
    fig = plt.figure(figsize=(5, 9))
    gs = GridSpec(
        nrows=20, 
        ncols=5, 
        figure=fig,
        wspace=0.05, # 水平间距（与原代码 wspace 匹配）
        hspace=0.2   # 垂直间距，可以根据需要调整
    )
    ax_heatmap = fig.add_subplot(gs[:18, :4])
    ax_curve = fig.add_subplot(gs[:18, 4], sharey=ax_heatmap)
    cax = fig.add_subplot(gs[19, :4])
    im = ax_heatmap.imshow(df_zscore.values, cmap="RdBu_r", aspect='auto', interpolation='nearest')
    ax_heatmap.invert_yaxis()
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xticks(np.arange(df_zscore.shape[1]))
    ax_heatmap.set_xticklabels(df_zscore.columns)
    ax_heatmap.xaxis.tick_top() 
    ax_heatmap.xaxis.set_label_position('top')

    # 定义注释的位置（Axes 坐标）
    axes_y_line = 1.07 # 横线位置（略微在顶部边界 1.0 上方）
    axes_y_text = 1.09  # 文本位置（在横线上方）
    line_color = 'black'

    # 获取 X 轴坐标的映射（data 坐标） imshow: 中心点坐标默认 0,1,2……,边界默认加减0.5
    x_min_ctrl = -0.1
    x_max_ctrl = 1.1
    x_center_ctrl = 0.5

    x_min_spz1 = 1.9
    x_max_spz1 = 3.1
    x_center_spz1 = 2.5

    # a. "Ctrl" 分组
    # 绘制横线：X 轴使用 Data 坐标，Y 轴使用 Axes 坐标
    ax_heatmap.plot(
        [x_min_ctrl, x_max_ctrl], 
        [axes_y_line, axes_y_line], 
        color=line_color, 
        linewidth=0.8,
        transform=ax_heatmap.get_xaxis_transform(), # 关键：X轴使用Data，Y轴使用Axes
        clip_on=False
    ) 
    # 放置文本：X 轴使用 Data 坐标，Y 轴使用 Axes 坐标
    ax_heatmap.text(
        x_center_ctrl, 
        axes_y_text, 
        "Ctrl", 
        ha='center', 
        va='center', 
        fontsize=10, 
        weight='normal',
        transform=ax_heatmap.get_xaxis_transform() # 关键：X轴使用Data，Y轴使用Axes
    ) 

    # b. "Spz1" 分组
    ax_heatmap.plot(
        [x_min_spz1, x_max_spz1], 
        [axes_y_line, axes_y_line], 
        color=line_color, 
        linewidth=0.8,
        transform=ax_heatmap.get_xaxis_transform(),
        clip_on=False
    ) 
    ax_heatmap.text(
        x_center_spz1, 
        axes_y_text, 
        "Spz1", 
        ha='center', 
        va='center', 
        fontsize=10, 
        weight='normal',
        transform=ax_heatmap.get_xaxis_transform()
    ) 

    y = np.arange(df_plot.shape[0])
    x = dfnew["score"]
    window_size = 50
    x_smooth = x.rolling(window=window_size, center=True).mean().fillna(x.rolling(window=window_size, center=True).mean().iloc[0])
    ax_curve.plot(
        x_smooth, y,
        color='red',
        linewidth=1.5
    )
    ax_curve.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    ax_curve.xaxis.tick_top()
    ax_curve.xaxis.set_label_position('top')
    # ax_curve.set_title("Moving average of BETA score", fontsize=10, pad=15) # 使用 pad 调整标题与刻度的距离
    ax_curve.text(
        0.5, 
        1.05, 
        "Moving average of\nBETA score", 
        linespacing=1.5,
        fontsize=10, 
        ha='center', 
        va='bottom', 
        transform=ax_curve.transAxes # 👈 使用曲线轴的 Axes 坐标
    )

    fig.colorbar(
        im, 
        cax=cax, # 👈 使用 GridSpec 创建的专用 cax
        label=None, 
        orientation='horizontal'
    )
    fig.savefig("/disk5/luosg/BETA20250928/output/heatmap/minus.png", dpi=300, bbox_inches='tight')

def plot_heatmap_with_group_colors(
        heatmap_data:pd.DataFrame,
        outplot:str,
        row_rename_map: dict,
        column_rename_map: dict,
        show_group_labels:bool = True,   # 是否显示分组文字
        title:str = "Heatmap of Positions",
        xlabel:str = "Positions",
        colorbar_width:float = 0.02,
        x_annotation_height:float = 0.2
):
    ### data preprocessing
    heatmap_data.rename(index=row_rename_map, inplace=True)
    heatmap_data.rename(columns=column_rename_map, inplace=True)
    heatmap_data = heatmap_data.sort_index()
    new_col_order = [new_name for new_name in column_rename_map.values()]
    heatmap_data = heatmap_data[new_col_order]
    group_colors = {}
    for group in heatmap_data.index:
        if group not in group_colors:
            group_colors[group] = sample(
                [
                    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                    '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
                ],
                1
            )[0]
    
    n_rows = heatmap_data.shape[0]
    fig = plt.figure(figsize=(7, 11))
    gs = gridspec.GridSpec(2, 2, width_ratios=[colorbar_width, 1],height_ratios = [1,x_annotation_height],wspace=0,hspace=0.3)
    # ---------------------
    # 颜色条
    ax_bar = fig.add_subplot(gs[0, 0])
    group_positions = {}
    for idx, group in enumerate(heatmap_data.index):
        ax_bar.add_patch(
            plt.Rectangle((0.6, idx), 0.3, 1, facecolor=group_colors[group])
        )
        if group not in group_positions:
            group_positions[group] = []
        group_positions[group].append(idx)
    # 添加分组标签
    if show_group_labels:
        for group, positions in group_positions.items():
            mid_idx = positions[len(positions)//2]
            ax_bar.text(
                0, mid_idx + 0.5, group,  # 文字在颜色条左边
                ha='right', va='center',
                fontsize=8
            )
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, n_rows)
    ax_bar.invert_yaxis()
    ax_bar.axis('off')
    # ---------------------

    # ---------------------
    # 热图
    ax_heat = fig.add_subplot(gs[0, 1])
    sns.heatmap(
        heatmap_data,
        ax=ax_heat,
        cmap='coolwarm',
        linewidths=0.01,
    )
    for label in ax_heat.get_xticklabels():
        label.set_rotation(60)
        label.set_fontsize(8)   # ★ 设置大小
        label.set_ha('center')
    ax_heat.set_yticks([])  # 去掉 y 轴
    ax_heat.set_ylabel('')
    ax_heat.set_xlabel(xlabel)
    ax_heat.set_title(title)
    cbar = ax_heat.collections[0].colorbar
    cbar.ax.text(
        0.5, 1.02,  # x=0.5 colorbar 中心, y=1.05 颜色条上方
        'ES',
        ha='center', va='bottom',
        rotation=0, fontsize=8,
        transform=cbar.ax.transAxes
    )
    # ---------------------

    # ---------------------
    # 在热图下方添加 x 轴分组图例
    ax_legend = fig.add_subplot(gs[1, :]) # 下方图例，与热图同列
    ax_legend.axis('off')
    short_names = list(heatmap_data.columns)
    pathway_rename_map = {v:k for k,v in column_rename_map.items()}
    legend_texts = [
        f"{abbr}: {pathway_rename_map[abbr]}"
        for abbr in short_names
        if abbr in pathway_rename_map
    ]

    max_chars_per_row = 60  # 每行最大字符数
    rows = []
    current_row = []
    current_len = 0

    for text in legend_texts:
        text_len = len(text)
        if current_len + text_len > max_chars_per_row and current_row:
            # 当前行满了，换行
            rows.append(current_row)
            current_row = [text]
            current_len = text_len
        else:
            current_row.append(text)
            current_len += text_len + 2  # 预留间隔
    if current_row:
        rows.append(current_row)
    n_rows = len(rows)
    # 绘制文字，左对齐
    for row_idx, row in enumerate(rows):
        x_pos = 0  # 每行从左边开始
        y_pos = n_rows - row_idx
        for text in row:
            ax_legend.text(
                x_pos, y_pos,
                text,
                ha='left', va='center',  # 左对齐
                rotation=0,
                fontsize=8
            )
            x_pos += len(text) + 2  # 累加下一个文字起始位置

    # 设置坐标范围，让文字显示完整
    ax_legend.set_xlim(0, max_chars_per_row)
    ax_legend.set_ylim(0, len(rows))
    # ---------------------
    plt.savefig(outplot, dpi=300)
    plt.close()

def plot_heatmap_with_colgroup(
        data: pd.DataFrame,
        col_group_map: dict,
        row_labels: dict = None,
        outplot: str = "heatmap.png"
    ):
    """
    不聚类热图，x轴列名染色显示分组，列分组图例在颜色条下方，每个分组独占一行
    """
    if row_labels:
        data = data.rename(index=row_labels)

    cmap = sns.color_palette("RdBu_r", as_cmap=True)
    vmin, vmax = np.nanmin(data.values), np.nanmax(data.values)

    fig, ax = plt.subplots(figsize=(max(8, data.shape[1]*0.6), max(6, data.shape[0]*0.4)))

    # 右侧颜色条
    cbar_ax = fig.add_axes([0.88, 0.25, 0.03, 0.7])
    sns.heatmap(
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        ax=ax,
        cbar_ax=cbar_ax,
        linewidths=0.5,
        linecolor='gray'
    )

    # 列分组颜色
    col_groups = pd.Series(col_group_map)
    unique_groups = col_groups.unique()
    palette = sns.color_palette("Set2", n_colors=len(unique_groups))
    group_colors = dict(zip(unique_groups, palette))

    # x轴标签染色
    for xtick_label in ax.get_xticklabels():
        label = xtick_label.get_text()
        if label in col_groups:
            xtick_label.set_color(group_colors[col_groups[label]])
            xtick_label.set_rotation(45)
            xtick_label.set_horizontalalignment('right')

    # 调整子图区域，留出下边距和右边距
    fig.subplots_adjust(bottom=0.2, right=0.78, top=0.95)

    # 绘制分组图例，位于右侧颜色条下方，每个分组占一行
    legend_handles = [plt.Line2D([0], [0], color=color, lw=6) for color in group_colors.values()]
    fig.legend(
        handles=legend_handles,
        labels=list(group_colors.keys()),
        loc='center right',
        bbox_to_anchor=(1.0, 0.15),  # 下方偏移
        ncol=1,
        frameon=False,
        title=""
    )

    plt.savefig(outplot, dpi=300)
    plt.close()

def plot_heatmap_with_group_colors1(
        heatmap_data,
        outplot,
        row_rename_map: dict,
        column_rename_map: dict,
        show_group_labels: bool = True,
        cmap: str = "coolwarm",

        # unified GridSpec weights
        y_colorbar_width: float = 1,
        heatmap_width: float = 50,
        right_cbar_width: float = 2,

        heatmap_height: float = 50,
        x_colorbar_height: float = 0.6,

        # ---- NEW ----
        min_group_fraction_for_label: float = 0.03   # 小组不显示文字
):
    """
    Only show text for large groups.
    Small groups appear as legend on the far right.
    """

    # ---- 1. rename ----
    heatmap_data = heatmap_data.copy()
    heatmap_data.rename(index=row_rename_map, inplace=True)
    heatmap_data.rename(columns=column_rename_map, inplace=True)

    col_order = [c for c in column_rename_map.values() if c in heatmap_data.columns]
    heatmap_data = heatmap_data.loc[:, col_order]

    # ---- 2. group list ----
    row_groups = list(heatmap_data.index)
    n_rows = len(row_groups)

    # ---- 3. 统计每个组行数 ----
    group_sizes = {}
    for g in row_groups:
        group_sizes[g] = group_sizes.get(g, 0) + 1

    total_n = float(n_rows)

    # ---- 4. 按 group 大小降序排序 ----
    sorted_groups = sorted(group_sizes.keys(), key=lambda g: group_sizes[g], reverse=True)

    # reorder rows
    ordered_index = []
    for g in sorted_groups:
        ordered_index.extend([idx for idx in heatmap_data.index if idx == g])

    heatmap_data = heatmap_data.loc[ordered_index, :]

    # update row groups
    row_groups = list(heatmap_data.index)
    n_rows, n_cols = heatmap_data.shape

    # ---- 5. define colors ----
    color_pool = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
        '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
        '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
    ]

    unique_row_groups = sorted_groups
    row_color_map = {g: color_pool[i % len(color_pool)] for i, g in enumerate(unique_row_groups)}

    # col group colors
    uniq_col = list(dict.fromkeys(list(heatmap_data.columns)))
    palette = sns.color_palette("hls", len(uniq_col))
    col_color_map = {g: palette[i] for i, g in enumerate(uniq_col)}

    # ---- 6. row group blocks ----
    blocks = []
    start = 0
    for i in range(1, n_rows):
        if row_groups[i] != row_groups[i - 1]:
            blocks.append((row_groups[i - 1], start, i))
            start = i
    blocks.append((row_groups[-1], start, n_rows))

    # ---- 6.1 识别大组、小组 ----
    big_groups = [g for g in sorted_groups if group_sizes[g] / total_n >= min_group_fraction_for_label]
    small_groups = [g for g in sorted_groups if g not in big_groups]

    # ---- 7. GridSpec ----
    fig = plt.figure(figsize=(8, 16))
    gs = gridspec.GridSpec(
        2, 3,
        width_ratios=[y_colorbar_width, heatmap_width, right_cbar_width],
        height_ratios=[heatmap_height, x_colorbar_height],
        wspace=0,
        hspace=0
    )

    ax_ybar = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_cbar = fig.add_subplot(gs[0, 2])
    ax_blank = fig.add_subplot(gs[1, 0])
    ax_xbar = fig.add_subplot(gs[1, 1])
    ax_blank_r = fig.add_subplot(gs[1, 2])

    ax_blank.axis("off")
    ax_blank_r.axis("off")

    vmin = np.percentile(heatmap_data.values, 5)
    vmax = np.percentile(heatmap_data.values, 95)
    # ---- 8. heatmap ----
    im = ax_heat.imshow(
        heatmap_data.values, 
        aspect='auto', 
        origin='upper', 
        cmap=cmap,
        vmin = vmin,
        vmax=vmax)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])

    # ---- 9. colorbar ----
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.ax.text(0.5, 1.02, "log2(cpm+1)", ha='center', va='bottom', transform=cbar.ax.transAxes)

    # ---- 10. y-axis color strip ----
    ax_ybar.set_xlim(0, 1)
    ax_ybar.set_ylim(0, n_rows)
    ax_ybar.invert_yaxis()

    for g, s, e in blocks:
        ax_ybar.add_patch(
            patches.Rectangle((0, s), 1, e - s, facecolor=row_color_map[g], edgecolor='none')
        )

    # ---- 大 group：左侧显示文字 ----
    if show_group_labels:
        for g, s, e in blocks:
            if g in big_groups:     # only large groups
                ax_ybar.text(
                    -0.1, (s + e) / 2, g,
                    ha="right", va="center", fontsize=8
                )

    ax_ybar.axis("off")

    # ---- 11. x-axis color strip ----
    col_groups = list(heatmap_data.columns)
    ax_xbar.set_xlim(0, n_cols)
    ax_xbar.set_ylim(0, 1)

    for i, g in enumerate(col_groups):
        ax_xbar.add_patch(
            patches.Rectangle((i, 0), 1, 1, facecolor=col_color_map[g], edgecolor='none')
        )

    if show_group_labels:
        for g in uniq_col:
            pos = [i for i, x in enumerate(col_groups) if x == g]
            mid = (pos[0] + pos[-1] + 1) / 2
            ax_xbar.text(mid, 0, g, ha="center", va="top", rotation=90, fontsize=8)

    ax_xbar.axis("off")

    # ---- 12. 小 group legend（最右边） ----
    if len(small_groups) > 0:
        handles = []
        labels = []
        for g in small_groups:
            handles.append(plt.Rectangle((0, 0), 1, 1, color=row_color_map[g]))
            labels.append(g)

        ax_cbar.legend(
            handles, labels,
            title="Groups",
            loc="upper left",
            bbox_to_anchor=(5, 1),
            fontsize=8,
            title_fontsize=9
        )

    # ---- SAVE ----
    plt.savefig(outplot, dpi=300, bbox_inches='tight')
    plt.close(fig)

def long_to_wide(
        df: pd.DataFrame,
        index_col: str,
        columns_col: str,
        values_col: str
) -> pd.DataFrame:
    """
    Convert long format DataFrame to wide format for heatmap plotting.
    """
    wide_df = df.pivot(index=index_col, columns=columns_col, values=values_col)
    return wide_df

if __name__ == "__main__":
    df = pd.read_csv("/home/luosg/Data/genomeStability/output/result/gsva/gseapy.gene_set.gsva.report.csv")
    df_wide = long_to_wide(df, index_col='Name', columns_col='Term', values_col='ES')
    plot_heatmap(df_wide,"/home/luosg/Data/genomeStability/output/result/gsva/heatmap.png")
