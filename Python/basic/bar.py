import pandas as pd
from matplotlib import pyplot as plt
import re
import sys
Domestication = sys.argv[1]
Splendens = sys.argv[2]
Sibling = sys.argv[3]
outfile = sys.argv[4]
annovar = sys.argv[5]
outsnp = f'{outfile}.snp'
# Domestication = "/home/lsg/GWAS/Betta_splendens/sweepfinder/data/vcf/Domestication.Chr21.1580811_1588574.vcf"
# Splendens = "/home/lsg/GWAS/Betta_splendens/sweepfinder/data/vcf/Splendens.Chr21.1580811_1588574.vcf"
# Sibling = "/home/lsg/GWAS/Betta_splendens/sweepfinder/data/vcf/Sibling.Chr21.1580811_1588574.vcf"
def genotype(text):
# match genotype
# print(text)
    #有的基因型以|分割，且在相位中也可能展示基因型；所以第一个\t和:包括内容才为基因型
    pattern1 = re.compile(r'\t(\d.\d):')
    pattern2 = re.compile(r'(Chr\d*\t)(\d*\t)(.\t)([a-zA-Z]*\t)([a-zA-Z]*\t)')
    genotype = pattern1.findall(text)
    match = pattern2.match(text)
    print(match.groups())
    w1 = match.groups()[3].replace("\t","")
    w2 = match.groups()[4].replace("\t","")
    Chr = match.groups()[0].replace("\t","")
    snp = match.groups()[1].replace("\t","")
    print(w1)
    print(w2)
    # print(snp)
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
    # sta = {f'{w1}/{w1}':g1,f'{w1}/{w2}':g2,f'{w2}/{w2}':g3}
    sta = [Chr,snp,w1,w2,g1,g2,g3]
    return sta
def read_vcf(infile):
    fo = open(infile,'r')
    texts = fo.readlines()
    fo.close()
    print(texts[360])
    return texts
#n代表从第几行开始
 
def analysis(texts,n):
    df = pd.DataFrame(columns=['Chr','snp','ref','alt','ref/ref','ref/alt','alt/alt'])
    a=1
    # print(df)
    print(genotype(texts[361]))
    for i in range(n,len(texts)):
        df.loc[i,] = genotype(texts[i])
        a = a+1
    df['ref/ref'] = df['ref/ref'].astype(int)
    df['ref/alt'] = df['ref/alt'].astype(int)
    df['alt/alt'] = df['alt/alt'].astype(int)
    df.reset_index(drop=True,inplace=True)
    return df
#注释非同义突变snp
p=1
def read_annovar(infile1,df):
    df1 = pd.read_table(infile1,sep=",",header=None)
    df1 = df1[df1[1] == 'nonsynonymous SNV'].reset_index(drop=True)
    for i in range(len(df1[1])):
        snp = str(df1.iloc[i,5])
        snp_index = df.index[df['snp'] == snp].tolist()
        df.loc[snp_index,'ref/ref'] = df.loc[snp_index,'ref/ref']*1.5
        df.loc[snp_index,'ref/alt'] = df.loc[snp_index,'ref/alt']*1.5
        df.loc[snp_index,'alt/alt'] = df.loc[snp_index,'alt/alt']*1.5
        # a=df.loc[snp_index,'ref/ref']+df.loc[snp_index,'ref/alt']+df.loc[snp_index,'alt/alt']
        # print(type(a))
        print(snp)
        print(snp_index)
        # print(type(snp))
        global p  
        if p <= len(df1[1]):
            p=p+1
            with open(outsnp, 'a') as f:
                print(snp, file=f)
                print(snp_index,file=f,end="")
                print(df.loc[snp_index,'ref/ref'].values,file=f,end=" ")
                print(df.loc[snp_index,'ref/alt'].values,file=f,end=" ")
                print(df.loc[snp_index,'alt/alt'].values,file=f)
    return df


lines_Dn = read_vcf(Domestication)
df_Dn = analysis(lines_Dn,361)

lines_Ss = read_vcf(Splendens)
df_Ss = analysis(lines_Ss,361)

lines_Sg = read_vcf(Sibling)
df_Sg = analysis(lines_Sg,361)
#控制输出野生基因型，而不是其它
df_Ss = read_annovar(annovar,df_Ss)
df_Dn = read_annovar(annovar,df_Dn)
df_Sg = read_annovar(annovar,df_Sg)
fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(18, 9))
# fig,ax1  = plt.subplots(figsize=(18, 9))
ax1.bar(x=range(1,len(df_Dn['ref'])+1),height=df_Dn['ref/ref'])
ax1.bar(x=range(1,len(df_Dn['ref'])+1),height=df_Dn['ref/alt'],bottom=df_Dn['ref/ref'])
ax1.bar(x=range(1,len(df_Dn['ref'])+1),height=df_Dn['alt/alt'],bottom=df_Dn['ref/ref']+df_Dn['ref/alt'])
ax1.set_ylabel('Domestication')
ax2.bar(x=range(1,len(df_Ss['ref'])+1),height=df_Ss['ref/ref'])
ax2.bar(x=range(1,len(df_Ss['ref'])+1),height=df_Ss['ref/alt'],bottom=df_Ss['ref/ref'])
ax2.bar(x=range(1,len(df_Ss['ref'])+1),height=df_Ss['alt/alt'],bottom=df_Ss['ref/ref']+df_Ss['ref/alt'])
ax2.set_ylabel('Splendens')
ax3.bar(x=range(1,len(df_Sg['ref'])+1),height=df_Sg['ref/ref'])
ax3.bar(x=range(1,len(df_Sg['ref'])+1),height=df_Sg['ref/alt'],bottom=df_Sg['ref/ref'])
ax3.bar(x=range(1,len(df_Sg['ref'])+1),height=df_Sg['alt/alt'],bottom=df_Sg['ref/ref']+df_Sg['ref/alt'])
ax3.set_ylabel('Sibling')
print(outfile)
fig.savefig(outfile) 
# print(df_vcf.head())