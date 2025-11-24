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
