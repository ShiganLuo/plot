import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
infile = "/home/lsg/GWAS/Betta_splendens240901/haplotype/data/heatmap/haplo.txt"
ann = "/home/lsg/GWAS/Betta_splendens240901/haplotype/data/ID/Splendens615_Pop.txt"
df1 = pd.read_table(infile,header=0)
df2 = pd.read_table(ann,header=0)
df_merged = pd.merge(df1,df2,on='sample')
df_sorted = df_merged.sort_values(by='Pop')
heatmap_data = df_sorted.drop(columns=['sample', 'Pop']).set_index(df_sorted['sample'])
plt.figure(figsize=(10, 6))
sns.heatmap(heatmap_data, annot=True, cmap='coolwarm', cbar_kws={'label': 'Values'}, linewidths=0.5)
plt.title('Heatmap of Positions')
plt.xlabel('Positions')
plt.ylabel('Samples')
plt.savefig("heatmap.png")