#!/usr/bin/env Rscript

suppressMessages(library('ggplot2'))

setwd("/net/isi-dcnl/ifs/user_data/JJChen_Grp/zhoukr/project/zhhchen/tet1/rna_seq/publicData/expression")
geneExpMtx = 'BCells.genes.tpm.long.matrix'
# With the count matrix, cts, and the sample information, colData
cts <- data.frame(read.csv(geneExpMtx, sep="\t", row.names=NULL))
#filtered <- apply(cts, 1, function(x) length(x[x>1])>=round(sampleSize/2))
#cts <- cts[filtered,]

plotGeneCount <- function(geneId, geneName, type) {
  pdf(paste(geneName, ".", type ,".pdf", sep=""), paper="a4r")
  subSelect <- cts[cts$gene_id==geneId,]
  row.names(subSelect) = subSelect$sample
  #subSelect <- plotCounts(ddsTC, geneId, 
  #                   intgroup = c("cell","age"), transform=FALSE, returnData = TRUE)
  subSelect$cell_type <- factor(subSelect$cell,levels = c("earlyB", "proB", "preB", "immatureB"))
  print(ggplot(subSelect,
    aes(x = cell_type, y = expression, color = age, group = age)) + 
    geom_point(size=2) + stat_summary(fun.y=mean, geom="line", size = 2) +
    labs(y = paste("Gene Expression: ", type), x = "Cell Types"))
  garbage <- dev.off()
}

geneList <- list(list("ENSG00000138336.9", "TET1"), list("ENSG00000168769.13", "TET2"), list("ENSG00000187605.15", "TET3"), 
    list("ENSG00000145388.15", "METTL14"), list("ENSG00000198492.16", "YTHDF2"))

for (i in 1:length(geneList)) {
  geneId <- geneList[i][1][[1]][[1]]
  geneName <- geneList[i][1][[1]][[2]]
  plotGeneCount(geneId, geneName, "TPM")
}

geneExpMtx = 'BCells.genes.fpkm.long.matrix'
cts <- data.frame(read.csv(geneExpMtx, sep="\t", row.names=NULL))
for (i in 1:length(geneList)) {
  geneId <- geneList[i][1][[1]][[1]]
  geneName <- geneList[i][1][[1]][[2]]
  plotGeneCount(geneId, geneName, "FPKM")
}
