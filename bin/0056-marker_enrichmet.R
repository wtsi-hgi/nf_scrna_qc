library(clusterProfiler)
library(dplyr)
library(DOSE)
library(optparse)
library(ggplot2)
library(org.Hs.eg.db)

option_list <- list(
  make_option(c("--markers_table"), type = "character")
)

opt <- parse_args(OptionParser(option_list=option_list))
basename <- gsub(".tsv.gz", "", opt$markers_table)
markers <- read.table(file = opt$markers_table, sep = '\t', header = TRUE)

cl <- sort(unique(markers$cluster))

for(i in 1:length(cl)){
  d <- markers[which(markers$cluster==cl[i]),]
  geneList <- d$logfoldchanges
  eg <- bitr(
      d$gene_symbols,
      fromType="SYMBOL",
      toType="ENTREZID",
      OrgDb="org.Hs.eg.db"
  )
  names(geneList) <- eg$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  gene <- names(geneList)


  kk <- enrichKEGG(gene         = names(geneList),
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)


  ggo1 <- groupGO(gene    = gene,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "CC",
                 level    = 3,
                 readable = TRUE)

  ggo2 <- groupGO(gene    = gene,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "MF",
                 level    = 3,
                 readable = TRUE)

  ggo3 <- groupGO(gene    = gene,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = 3,
                 readable = TRUE)

  #Disease association
  edo <- enrichDGN(gene)
  edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')


  pdf(paste0(basename, "-enrichment_cluster_", cl[i], ".pdf"), onefile = TRUE)

  p1 <- barplot(kk, showCategory=20) + ggtitle("KEGG terms")
  p2 <- barplot(ggo1, showCategory=20) + ggtitle("Cellular Composition")
  p3 <- barplot(ggo2, showCategory=20) + ggtitle("Molecular Function")
  p4 <- barplot(ggo3, showCategory=20) + ggtitle("Biological Process")
  p5 <- barplot(edo, showCategory=20) + ggtitle("Disease association")
  p6 <- cnetplot(edox, node_label="all") + ggtitle("Disease association")

  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)

  dev.off()


}
