import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#  python scripts/sweepfinder_plot.py 
Crowntail = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Crowntail/Crowntail_chr"
outfile_cr = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Crowntail/plot/Crowntail.jpg"
Halfmoon = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Halfmoon/Halfmoon_chr"
outfile_hn = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Halfmoon/plot/Halfmoon.jpg"
HMPk = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/HMPK/HMPK_chr"
outfile_hk = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/HMPK/plot/HMPK.jpg"
Veiltail = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Veiltail/Veiltail_chr"
outfile_vl = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Veiltail/plot/Veiltail.jpg"
Fighter = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Fighter/Fighter_chr"
outfile_fr = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Fighter/plot/Fighter.jpg"
domestication = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/domestication/domestication_chr"
outfile_dn = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/domestication/plot/domestication.jpg"
LongFin = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/LongFin/LongFin_chr"
outfile_Ln = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/LongFin/plot/LongFin.jpg"
Splendens = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Splendens/Splendens_chr"
outfile_ss = "/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Splendens/plot/Splendens.jpg"
xticks = [0]*21
xlab = ["Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9",
        "Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr16","Chr17","Chr18","Chr19","Chr20","Chr21"]
def sweepfinder_plot(file_path,ax,xticks):
    tmp = [0]
    pos = [0]*22 #防止tmp被更新
    chr_df = [0]*21
    for i in range(1,22):
        infile = file_path + str(i) + "_SF_5K.out"
        # infile = f"/home/lsg/GWAS/Betta_splendens/sweepfinder/output/Crowntail/Crowntail_chr{i}_SF_5K.out"/
        # print(infile)
        df = pd.read_table(infile)
        tmp.append(df["location"].max()) #每组染色体序列最大值
        pos[i] = sum(tmp)
    print(tmp)
    for i in range(1,22):
        infile = file_path + str(i) + "_SF_5K.out"
        df = pd.read_table(infile)
        df = df[df['LR'] > 0]
        df['location'] = df['location'] + pos[i-1]
        df = df.sort_values("location")#避免存在染色体位置不是按顺序排列的情况，以防杂线出现
        xticks[i-1] = df["location"].median()
        ax.scatter(df['location'], df['LR'])
# ########       
# fig1, ax1 = plt.subplots(figsize=(20, 6))        
# sweepfinder_plot(Crowntail,ax1,xticks)
# ax1.set_xticks(xticks)
# ax1.set_xticklabels(xlab)
# ax1.set_xlabel("Chromosome")
# ax1.set_ylabel("LR")
# ax1.set_title("Crowntail")
# fig1.savefig(outfile_cr, bbox_inches='tight')
# plt.close(fig1)
# ##########
# fig2, ax2 = plt.subplots(figsize=(20, 6))
# sweepfinder_plot(Halfmoon,ax2,xticks)
# #执行函数后xticks自然会更新
# ax2.set_xticks(xticks)
# ax2.set_xticklabels(xlab)
# ax2.set_xlabel("Chromosome")
# ax2.set_ylabel("LR")
# ax2.set_title("Halfmoon")
# fig2.savefig(outfile_hn, bbox_inches='tight')
# plt.close(fig2)
# ########       
# fig3, ax3 = plt.subplots(figsize=(20, 6))        
# sweepfinder_plot(HMPk,ax3,xticks)
# ax3.set_xticks(xticks)
# ax3.set_xticklabels(xlab)
# ax3.set_xlabel("Chromosome")
# ax3.set_ylabel("LR")
# ax3.set_title("HMPK")
# fig3.savefig(outfile_hk, bbox_inches='tight')
# plt.close(fig3)
# ########       
# fig4, ax4 = plt.subplots(figsize=(20, 6))        
# sweepfinder_plot(Veiltail,ax4,xticks)
# ax4.set_xticks(xticks)
# ax4.set_xticklabels(xlab)
# ax4.set_xlabel("Chromosome")
# ax4.set_ylabel("LR")
# ax4.set_title("Veiltail")
# fig4.savefig(outfile_vl, bbox_inches='tight')
# plt.close(fig4)
# fig5, ax5 = plt.subplots(figsize=(20, 6))        
# sweepfinder_plot(Fighter,ax5,xticks)
# ax5.set_xticks(xticks)
# ax5.set_xticklabels(xlab)
# ax5.set_xlabel("Chromosome")
# ax5.set_ylabel("LR")
# ax5.set_title("Fighter")
# fig5.savefig(outfile_fr, bbox_inches='tight')
# plt.close(fig5)

# fig6, ax6 = plt.subplots(figsize=(20, 6))        
# sweepfinder_plot(domestication,ax6,xticks)
# ax6.set_xticks(xticks)
# ax6.set_xticklabels(xlab)
# ax6.set_xlabel("Chromosome")
# ax6.set_ylabel("LR")
# ax6.set_title("domestication")
# fig6.savefig(outfile_dn, bbox_inches='tight')
# plt.close(fig6)
# fig7, ax7 = plt.subplots(figsize=(20, 6))        
# sweepfinder_plot(LongFin,ax7,xticks)
# ax7.set_xticks(xticks)
# ax7.set_xticklabels(xlab)
# ax7.set_xlabel("Chromosome")
# ax7.set_ylabel("LR")
# ax7.set_title("LongFin")
# fig7.savefig(outfile_Ln, bbox_inches='tight')
# plt.close(fig7)
fig8, ax8 = plt.subplots(figsize=(20, 6))        
sweepfinder_plot(Splendens,ax8,xticks)
ax8.set_xticks(xticks)
ax8.set_xticklabels(xlab)
ax8.set_xlabel("Chromosome")
ax8.set_ylabel("LR")
ax8.set_title("Splendens")
fig8.savefig(outfile_ss, bbox_inches='tight')
plt.close(fig8)