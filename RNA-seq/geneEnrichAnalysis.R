#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",      "h",  0,   "logical",      "Show help information",
    "category",  "c",  1,   "character",    "hallmark|c2-7 (https://www.gsea-msigdb.org/gsea)",
    "dbName",    "d",  1,   "character",    "Used database: org.Hs.eg.db|org.Mm.eg.db",
    "geneset",   "s",  1,   "character",    "geneset (https://www.gsea-msigdb.org/gsea)",
    "gotype",    "y",  1,   "character",    "GO type: ALL|BP|MF|CC",
    "input",     "i",  1,   "character",    "2-column input file (NO HEADER, 1stCol:gene_id, 2ndCol:value)",
    "item",      "e",  1,   "character",    "Type of item: ENSEMBL|SYMBOL|ENTREZID",
    "max",       "a",  1 ,  "numeric",      "Maximum items in pdf [30]",
    "minGSSize", "n",  1 ,  "numeric",      "Minimal size of genes annotated for testing [10]",
    "maxGSSize", "z",  1 ,  "numeric",      "Maximal size of genes annotated for testing [500]",
    "nperm",     "r",  1 ,  "numeric",      "Number of permutations (GSEA) [1000]",
    "organism",  "g",  1,   "character",    "Organism: human|mouse",
    "output",    "o",  1,   "character",    "Output directory [./]",
    "prefix",    "f",  1,   "character",    "Prefix of output file [geneEnrich]",
    "pval",      "p",  1,   "numeric",      "pval cutoff [0.05]",
    "qval",      "q",  1,   "numeric",      "qval cutoff [0.2]",
    "type",      "t",  1,   "character",    "Type of analysis: GO|KEGG|GSEA"
), byrow=TRUE, ncol=5)

# function
ShowHelp <- function(object, param, reverse=FALSE, bool=FALSE) {
  if (reverse) {
    if (bool) {
      judge <- !isTRUE(object)
    }else{
      judge <- !is.null(object)
    }
  }else{
    if (bool) {
      judge <- isTRUE(object)
    }else{
      judge <- is.null(object)
    }
  }
  if (judge) {
    if (param != 'none') {
      cat(paste("None valid ", param, "!\n"))
    }
    cat(paste(getopt(command, usage = T),"\n"))
    q()
  }
}

LoadPacakge <- function(name) {
  suppressMessages(library(name, character.only = TRUE))
  cat(paste("Load package: ", name, ".\n"))
}

WriteText <- function(file, enrichResult, type=NULL, reverse=FALSE, bool=FALSE) {
  # output text file
  output.file <- file(file, "wb")
  write.table(as.data.frame(enrichResult), sep="\t", eol = "\n", 
              quote = FALSE, row.names=TRUE, file=output.file)
  close(output.file)
  if (! is.null(type)) {
    if (type == 'GSEA') {
      scan(pipe(paste("sed -i '1 s/ID/ID1\tID2/' ", file, sep = "")))
    }
  }
}

args <- getopt(command)

# argument check and assignment
ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$input, '-i|--input')
ShowHelp(args$item, '-e|--item')
ShowHelp(args$type, '-t|--type')
ShowHelp(args$organism, '-g|--organism')

# check variables
itemVector <- c('ENSEMBL', 'SYMBOL', 'ENTREZID')
if ( isFALSE(args$item %in% itemVector) ) {
  ShowHelp(args$item, '-e|--item')
}

methodVector <- c('GO', 'KEGG', 'GSEA')
if ( isFALSE(args$type %in% methodVector) ) {
  ShowHelp(args$type, '-t|--type')
}

organismVector <- c('human', 'mouse')
if ( isFALSE(args$organism %in% organismVector) ) {
  ShowHelp(args$organism, '-g|--organism')
}

categoryVector <- c('hallmark', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7')
if (args$type == 'GSEA') {
  ShowHelp(args$category, '-c|--category')
  if ( isFALSE(args$category %in% categoryVector) ) {
    ShowHelp(args$category, '-c|--category')
  }
}

gotypeVector <- c('ALL', 'BP', 'MF', 'CC')
if (args$type == 'GO') {
  ShowHelp(args$gotype, '-y|--gotype')
  if ( isFALSE(args$gotype %in% gotypeVector) ) {
    ShowHelp(args$gotype, '-y|--gotype')
  }
}

# default arguments
if ( is.null(args$pval) ) {
  args$pval <- 0.05
}

if ( is.null(args$qval) ) {
  args$qval <- 0.2
}

if ( is.null(args$max) ) {
  args$max <- 30
}

if ( is.null(args$minGSSize) ) {
  args$minGSSize <- 10
}

if ( is.null(args$maxGSSize) ) {
  args$maxGSSize <- 500
}

if ( is.null(args$nperm) ) {
  args$nperm <- 1000
}

if ( is.null(args$prefix) ) {
  args$prefix <- 'geneEnrich'
}

if ( is.null(args$output) ) {
  args$output <- './'
}

if (args$organism == 'human'){
  dbName <- 'org.Hs.eg.db'
}else if (args$organism == 'mouse') {
  dbName <- 'org.Mm.eg.db'
}

if (args$type == 'KEGG') {
  if (args$organism == 'human'){
    organism <- 'hsa'
  }else if (args$organism == 'mouse') {
    organism <- 'mmu'
  }
}

# create output directory
if (args$type == "GO" || args$type == "KEGG") {
  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
}

## main program
LoadPacakge('ggplot2')
LoadPacakge('clusterProfiler')
LoadPacakge(dbName)

geneMtx <- read.table(args$input, sep="\t", header=F)
names(geneMtx) <- c("id", "value")
## assume 1st column is ID
## 2nd column is FC

if (args$item == 'ENSEMBL') {
  ## removing gene version
  ## ENSG00000262662.1 -> ENSG00000262662
  geneMtx$id <- gsub("\\..+$", "", geneMtx$id)
}

if (args$item != 'ENTREZID' ) {
  # convert item to ENTREZID and filter
  geneBitr <- bitr(geneMtx$id, fromType = args$item, toType = "ENTREZID", OrgDb = dbName)
  names(geneBitr)[1] <- paste("oldItem")
  # merge geneBitr and geneMtx
  combineData <- merge(geneBitr, geneMtx, by.x='oldItem', by.y="id")
  geneBitr <- bitr(combineData$ENTREZID, fromType = 'ENTREZID', toType = "SYMBOL", OrgDb = dbName)
  ## becomes: ENTREZID SYMBOL oldItem value
  combineData <- merge(geneBitr, combineData, by.x="ENTREZID", by.y="ENTREZID")
}else{
  ## filter non-valided ids
  geneBitr <- bitr(geneMtx$id, fromType = args$item, toType = "SYMBOL", OrgDb = dbName)
  combineData <- merge(geneBitr, geneMtx, by.x="ENTREZID", by.y="id")
  combineData <- combineData[c("ENTREZID", "SYMBOL", "value")]
}

# feature 1: numeric vector
geneList <- combineData$value

# feature 2: named vector
if (args$type == "GO") {
  names(geneList) <- as.character(combineData$ENTREZID)
}else if (args$type == "KEGG") {
  names(geneList) <- as.character(combineData$ENTREZID)
}else{
  names(geneList) <- as.character(combineData$SYMBOL)
}

# feature 3: decreasing orde
geneList <-sort(geneList, decreasing = TRUE)

# get gene list
gene <- names(geneList)

# pdf name
dotplotPdf = file.path(args$output, paste(args$prefix, ".dotplot.pdf", sep=""))
#text file
textFile = file.path(args$output, paste(args$prefix, ".geneEnrich.txt", sep=""))

# enrichment analysis
if (args$type == "GO") {
  ## run enrichGO
  geneEnrich <- enrichGO(
    gene,
    OrgDb = dbName,
    ont = args$gotype,
    keyType = "ENTREZID",
    minGSSize = args$minGSSize,
    maxGSSize = args$maxGSSize,
    pvalueCutoff = args$pval,
    pAdjustMethod = "BH",
    qvalueCutoff = args$qval,
    readable = FALSE
  )
  ## print to pdf
  pdf(dotplotPdf, paper='a4r', height=0)
  dotplot(geneEnrich, showCategory = args$max)
  garbage <- dev.off()
  WriteText(textFile, geneEnrich)
}else if (args$type == "KEGG") {
  ## run enrichKEGG
  geneEnrich <- enrichKEGG(
    gene,
    keyType = 'kegg',
    minGSSize = args$minGSSize,
    maxGSSize = args$maxGSSize,
    organism = organism,
    pvalueCutoff = args$pval,
    pAdjustMethod = "BH",
    qvalueCutoff = args$qval
  )
  ## print to pdf
  pdf(dotplotPdf, paper='a4r', height=0)
  dotplot(geneEnrich, showCategory = args$max)
  garbage <- dev.off()
  WriteText(textFile, geneEnrich)
}else if (args$type == "GSEA") {
  LoadPacakge('dplyr')
  LoadPacakge('msigdf')
  LoadPacakge('enrichplot')
  ## get MSigDB dataset from msigdf
  if (args$organism == 'human' ) {
    gseaMSigDB <- msigdf.human %>% filter(category_code == args$category) %>% select(geneset, symbol) %>% as.data.frame
  }else if (args$organism == 'mouse') {
    gseaMSigDB <- msigdf.mouse %>% filter(category_code == args$category) %>% select(geneset, symbol) %>% as.data.frame
  }
  if (! is.null(args$geneset)) {
    gseaMSigDB <- gseaMSigDB[gseaMSigDB$geneset == args$geneset, ]
  }
  ## run GSEA
  geneEnrich <- GSEA(
    geneList,
    nPerm = args$nperm,
    minGSSize = args$minGSSize,
    maxGSSize = args$maxGSSize,
    TERM2GENE = gseaMSigDB,
    pvalueCutoff = args$pval,
    pAdjustMethod = "BH",
    seed = TRUE,
    by = "DOSE"
  )
  ## print to pdf
  genesetVector <- geneEnrich$ID
  if (length(genesetVector) > 0) {
    dir.create(args$output, showWarnings = FALSE, recursive = TRUE)
    for (i in 1:length(genesetVector)) {
      geneset <- genesetVector[i]
      pvalue <- formatC(geneEnrich$pvalue[i], format = "e", digits = 2)
      padjust <- formatC(geneEnrich$p.adjust[i], format = "e", digits = 2)
      title <- paste("Statistical significance ", '(pvalue:', pvalue, ', p.adjust:', padjust, ')', sep="")
      ## print to pdf
      pdfName <- paste(args$prefix, '.', geneset, ".gsea.pdf", sep="")
      pdf(file.path(args$output, pdfName), paper='a4r', height=0)
      print(gseaplot2(geneEnrich, geneSetID=i, title = title, base_size=10))
      dev.off()
    }
    WriteText(textFile, geneEnrich, type='GSEA')
  }
}
