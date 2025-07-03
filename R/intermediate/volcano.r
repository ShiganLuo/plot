volcano <- function (res, outjpeg,mode = "gene",geneAnnotation = "", nlabel = 10, label.by = "padj"){
  # add gene_name column to the results table
  if ( mode == 'TE'){
    res = res %>% 
    rownames_to_column(var = "index") %>%
    mutate(gene_name = word(index, 1, sep = ":"))
  } else if ( mode == 'gene' && geneAnnotation != ""){
    df <- read.csv(geneAnnotation, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
    # it is gene mode as long as it contains ensembl_id ,this can resolve that both str1:str2:…… ensembl_id
    res <- res %>%
      rownames_to_column(var = "index") %>%
      mutate(gene_name = ifelse(
      str_count(index, ":") >= 2, 
      str_extract(index, "^[^:]+"),  # else extract gene_name from df
      df$gene_name[match(index, rownames(df))]))
  } else if ( mode == 'gene' && geneAnnotation == ""){
    stop("Invalid geneAnnotation argument. Please provide a gene annotation file, because you choose mode gene.")
  } else {
    stop("Invalid mode argument. Choose either gene or TE.")
  }

    # print(head(res))
  # assign significance to results based on padj
  res <- res %>%
  mutate(
    gene_status = case_when(
      log2FoldChange > 0.58 & -log10(padj) > 1.3 ~ "Upregulated genes",
      log2FoldChange < -0.58 & -log10(padj) > 1.3 ~ "Downregulated genes",
      TRUE ~ "Other genes"
    )
  )
  res = res[!is.na(res$padj),]
  # Calculate the number of each gene_status
  status_counts <- res %>%
    group_by(gene_status) %>%
    summarise(count = n())

  # Create new legend labels containing quantities
  new_labels <- setNames(
    sapply(status_counts$gene_status, function(status) {
      count <- status_counts$count[status_counts$gene_status == status]
      paste0(status, " (n = ", count, ")")
    }),
    status_counts$gene_status
  )
  # print(head(res))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange
  if (label.by == "padj") {
    up_genes <- res %>%
      filter(gene_status == "Upregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
    down_genes <- res %>%
      filter(gene_status == "Downregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
  } else if (label.by == "log2FoldChange") {
    up_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], desc(log2FoldChange)),nlabel)
    down_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], log2FoldChange),nlabel)
  } else
    stop ("Invalid label.by argument. Choose either padj or log2FoldChange.")
  
    p = ggplot(res, aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col=gene_status)) + 
      # scale_fill_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "black"), 
      #                       name = "gene classification", 
      #                       labels = c("Upregulated genes" = uplegendLabel, "Downregulated genes" = downlegendLabel,"Other genes" = otherlegendLabel))+
      scale_color_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "black"),
        labels = new_labels) + 
      ggrepel::geom_text_repel(data=up_genes, aes(label=head(gene_name,nlabel)), color = "#F59494", size = 3)+
      ggrepel::geom_text_repel(data=down_genes, aes(label=head(gene_name,nlabel)), color = "#93ACF6", size = 3)+
      labs ( x = expression(log[2]("FoldChange")), y = expression(-log[10]("adjusted p-value")))+
      geom_vline(xintercept = 0.58, linetype = "dotted")+
      geom_vline(xintercept = -0.58, linetype = "dotted")+
      geom_hline(yintercept = 1.3, linetype = "dotted")+
      theme_classic() +  
      theme(
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA,linewidth = 0.2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()   
      )
  jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
    print(p)
  dev.off()
}