tbl=read.table("/home/lsg/GWAS/Betta_splendens/group_structure/output/Fighter.5.Q")
jpeg(width = 150, height = 150, units = "mm", res = 300, file = "/home/lsg/GWAS/Betta_splendens/group_structure/output/Caudal_Dorsal.jpeg")
barplot(t(as.matrix(tbl)), col=rainbow(5),xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()