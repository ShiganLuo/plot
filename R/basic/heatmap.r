ann="/home/lsg/GWAS/Betta_splendens240901/haplotype/data/ID/Splendens615_Pop.txt"
input="/home/lsg/GWAS/Betta_splendens240901/haplotype/data/heatmap/sort_haplot.csv"
#加载所需的软件包
library(pheatmap)
#加载注释内容
annotation_row <- read.table(ann,header=T,row.names = 1,stringsAsFactors = T)
head(annotation_row)
#加载haplotype文件
hap <- read.table(input, header=T,row.names = 1, stringsAsFactors = F)
# # # color # 根据自己的调整

group = list(Pop = c("Stiktos" = "#735380","Smaragdina_guitar" = "#660000","Smaragdina" = "#8A2BE2","Siamorientalis" = "#20B2AA",
"Mahachaiensis" = "#FA8072","Imbellis" = "#E1FF00","Splendens" = "#BC80BD","Yellow" = "#E6D214","White" = "#EFDFBB",
"Veiltail" = "#AED0CD","Turquoise" = "#2ED29F","Steel" = "#4682B4","Royal" = "#0F55F1","Red" = "#C63731","Orange" = "#F7953B",
"Koi" = "#DE77AE","Halfmoon" = "#A2B61E" ,"Giant" = "#6F5CCC","Fighter" = "#907C62","Dumbo_HMPK" = "#F6D5A4",
"Dumbo_Halfmoon" = "#C47664","Dragon" = "#177D7A","Crowntail" = "#91FF00","Copper" = "#00FFF7","Black" = "#BEBFC5"))


p <- pheatmap(hap, show_colnames=FALSE, 
                legend_breaks = -1:2, 
                legend_labels = c("./.", "Reference", "Heterozygote", "Altemative"), 
                # color = c("#FFFFFF","#F9EFE5","#F6C8E0","#F472B5"),
                color = c("#FFFFFF","#21E706","#0C13D7","#D7270C"),
                #cutree_cols=2, 
                cluster_row = FALSE, 
                cluster_col = FALSE,
                annotation_row = annotation_row, 
                show_rownames = FALSE,
                annotation_colors = group, 
                # annotation_names_row = TRUE,
                # annotation_legend = FALSE,
                # labels_row = labels_row,
                main="chr3:0-5.21Mb",
                fontsize = 20,
                fointsize= 30,
                )

#pdf(file = "AAA.pdf", width = 10, height = 5)
png(file = paste0("chr3:0-5.21Mb",".png"), width = 6000, height = 4000,res = 300)
p
dev.off()
print("结束")