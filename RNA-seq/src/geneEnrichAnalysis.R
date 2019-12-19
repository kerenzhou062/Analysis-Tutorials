#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",     "h",  0,   "logical",      "Show help information",
    "method",   "me", 1,   "character",    "GO|KEGG",
    "dbName",   "db", 1,   "character",    "org.Hs.eg.db|org.Mm.eg.db",
    "item",     "it", 1,   "character",    "ENSEMBL|SYMBOL|ENTREZID",
    "pval",     "p",  1,   "numeric",      "pval cutoff, default is 0.05",
    "qval",     "q",  1,   "numeric",      "qval cutoff default is 0.2",
    "input",    "i",  1,   "character",    "GeneList input file",
    "organism", "or", 1,   "character",    "hsa|mmu(kegg)",
    "pdf",      "o",  1,   "character",    "Output pdf file",
    "type",     "t",  1,   "character",    "BP|MF|CC",
    "text" ,    "te", 1,   "character",    "Output the items in txt file",
    "max" ,     "m",  1 ,  "numeric",      "Maximum items in pdf"
), byrow=TRUE, ncol=5)

args <- getopt(command)

if(!is.null(args$help)){
    cat(paste(getopt(command, usage = T),"\n"))
    q()
}

# isRightType = sum(apply(cbind(c("all","up","down"),rep(args$type, length(c("all","up","down")))), 1, function(x){if(x[1]==x[2]){return(1)}else{return(0)}}))
# if(!isRightType){
#     cat(paste(getopt(command, usage = T),"\n"))
#     print("type has to be TPM,RPM,FPKM")
#     q()
# }


suppressMessages(library('ggplot2'))
suppressMessages(library('clusterProfiler'))
suppressMessages(library(args$dbName, character.only = TRUE))

gene <- read.table(args$input,sep="\t",header=F)

# if(args$type == "up"){
#     genelist = as.character(gene[gene$V2>1,]$V1)
# }else if(args$type == "down"){
#     genelist = as.character(gene[gene$V2<1,]$V1)
# }else{
#     genelist = as.character(gene$V1)
# }

genelist <- as.character(gene$V1)
geneID <- bitr(genelist, fromType = args$item, toType = "ENTREZID", OrgDb = args$dbName)

if (args$method == "GO") {
  geneEnrich <- enrichGO(
      geneID$ENTREZID,
      OrgDb = args$dbName,
      ont = args$type,
      keyType = "ENTREZID",
      pvalueCutoff = args$pval,
      pAdjustMethod = "BH",
      qvalueCutoff = args$qval,
      readable = TRUE
    )
}else if (args$method == "KEGG") {
  geneEnrich <- enrichKEGG(
      geneID$ENTREZID,
      organism = args$organism,
      pvalueCutoff = args$pval,
      pAdjustMethod = "BH",
      qvalueCutoff = args$qval
    )
}

# exl = gsub("\\.pdf",".txt",args$out)
#h <- 5 + args$max/10
pdf(args$pdf,width=10, heigh=5)
    dotplot(geneEnrich, showCategory = args$max)
dev.off()

write.table(geneEnrich, args$text, sep = "\t")

#print("done")

