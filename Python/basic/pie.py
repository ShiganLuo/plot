import re
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
# reda file 
infile = "/home/lsg/GWAS/Betta_splendens240901/gene_Investication/output/snp/dip2a.txt"
outdir = "/home/lsg/GWAS/Betta_splendens240901/gene_Investication/output/plot/pie"
fo = open(infile,'r')
texts = fo.readlines()
fo.close()
def genotype(text):
# match genotype
# print(text)
    #有的基因型以|分割，且在相位中也可能展示基因型；所以第一个\t和:包括内容才为基因型
    pattern1 = re.compile(r'\t(\d.\d):')
    pattern2 = re.compile(r'(\d*\t)(\d*\t)(.\t)([a-zA-Z]*\t)([a-zA-Z]*\t)')
    genotype = pattern1.findall(text)
    match = pattern2.match(text)
    print(match.groups())
    w1 = match.groups()[3].replace("\t","")
    w2 = match.groups()[4].replace("\t","")
    print(w1)
    print(w2)
    # print(genotype)
    # print(len(genotype))
    # plot
    #w1:0；w2:1
    g1 = genotype.count('0/0') + genotype.count('0|0')
    g2 = genotype.count('0/1') + genotype.count('1/0') + genotype.count('0|1') + genotype.count('1|0')
    g3 = genotype.count('1/1') + genotype.count('1|1')
    print(f'{w1}{w1}:{g1}')
    print(f'{w1}{w2}:{g2}')
    print(f'{w2}{w2}:{g3}')
    sta = {f'{w1}/{w1}':g1,f'{w1}/{w2}':g2,f'{w2}/{w2}':g3}
    return sta
####饼图
fig,(ax1,ax2,ax3) = plt.subplots(1, 3, sharex=True,figsize=(10, 6))
Domestication = genotype(texts[1])
print(Domestication)
Dvalue = list(Domestication.values())
Dlabel = list(Domestication.keys())
Splendens = genotype(texts[3])
print(Splendens)
Svalue = list(Splendens.values())
Sibling = genotype(texts[5])
Sivalue = list(Sibling.values())
print(Sibling)
ax1.pie(Dvalue,colors=['red','green','blue'])
ax2.pie(Svalue,colors=['red','green','blue'])
ax3.pie(Sivalue,colors=['red','green','blue'])
fig.legend(Dlabel, loc="upper right")
#标记需手动修改
ax1.text(-0.4,1.2,"Domestication")
ax1.text(-0.1,1.4,"479")
ax2.text(-0.35,1.2,"Splendens")
ax2.text(-0.1,1.4,"27")
ax3.text(-0.35,1.1,"Relatives")
ax3.text(-0.1,1.4,"93")
fig.text(0.45,0.05,"Chr21:10201624")
fig.savefig(f'{outdir}/PLA2G1B_Chr21_10201624.jpg')
# ##柱状图
# category = list()
# geno = list()
# g_color = ['red','blue','green']
# g_label = []

# for text in texts:
#     if len(text) < 20:
#         category.append(text)
#         # continue
#     else:
#         mid_sta = genotype(text)
#         geno.append(mid_sta)
#         g_label = mid_sta.keys()
# print(category)
# print(geno)
# g_label = list(g_label)
# print(g_label)
# #绘图
# fig,(ax2,ax1) = plt.subplots(2, 1, sharex=True,figsize=(10, 6))
# ax1.set_ylim(0,20)#显示较小数据
# ax2.set_ylim(20,400)#显示较大数据
# ax2.set_xticks([])#隐藏第二个轴的x轴标签
# legend_elements = [Patch(facecolor=g_color[i], label=g_label[i]) for i in range(len(g_label))]
# for i in range(len(geno)):
#     ge1 = geno[i][g_label[0]]
#     ge2 = geno[i][g_label[1]]
#     ge3 = geno[i][g_label[2]]
#     print(ge1)
#     category[i] = category[i].replace("\n","")
#     print(category[i])
#     ax1.bar(category[i],ge1,color = g_color[0])
#     ax1.bar(category[i],ge2,bottom=ge1,color = g_color[1])
#     ax1.bar(category[i],ge3,bottom=ge1+ge2,color = g_color[2])
#     ax2.bar(category[i],ge1,color = g_color[0])
#     ax2.bar(category[i],ge2,bottom=ge1,color = g_color[1])
#     ax2.bar(category[i],ge3,bottom=ge1+ge2,color = g_color[2])
#     yval = ge1+ge2+ge3  # 获取柱子的高度
#     ax1.text(i, 15, yval, ha='center', va='bottom')
# ax1.set_xticks(ticks=range(len(category)), labels=category)
# # plt.tight_layout() 
# # d = .015  # 截断线的大小
# # kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
# # ax1.plot((-d, +d), (-d, +d), **kwargs)  # 左下角
# # ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # 右下角

# # kwargs.update(transform=ax2.transAxes)  # 切换到第二个轴的坐标系
# # ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # 左上角
# # ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # 右上角
# ax2.legend(handles=legend_elements, loc='upper right', frameon=False)
# fig.savefig('./bar.jpg')