import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.interpolate import make_interp_spline
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec

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