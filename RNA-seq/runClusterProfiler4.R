#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",      "h",  0,   "logical",      "Show help information",
    "by",        "b",  1,   "character",    "by method in GSEA: fgsea|DOSE",
    "category",  "c",  1,   "character",    "hallmark|C2-C7 (https://www.gsea-msigdb.org/gsea)",
    "gotype",    "y",  1,   "character",    "GO type: ALL|BP|MF|CC",
    "geneset",   "l",  1,   "character",    "user input geneset matrix file (2 columns, gs_name,entrez_gene)",
    "input",     "i",  1,   "character",    "2-column input file (NO HEADER+TAB, 1stCol:gene_id, 2ndCol:value)",
    "item",      "e",  1,   "character",    "Type of item: ENSEMBL|SYMBOL|ENTREZID",
    "max",       "a",  1 ,  "numeric",      "Maximum items in pdf [20]",
    "minGSSize", "n",  1 ,  "numeric",      "Minimal size of genes annotated for testing [10]",
    "maxGSSize", "x",  1 ,  "numeric",      "Maximal size of genes annotated for testing [500]",
    "noCnet",    "k",  0,   "logical",      "skip GO cnet plots",
    "noHeatmap", "d",  0,   "logical",      "skip GO heatmap plots",
    "nperm",     "r",  1 ,  "numeric",      "Number of permutations (GSEA) [1000]",
    "organism",  "g",  1,   "character",    "Organism: human|mouse",
    "output",    "o",  1,   "character",    "Output directory [./]",
    "prefix",    "f",  1,   "character",    "Prefix of output file [geneEnrich]",
    "padjust",   "j",  1,   "character",   "The method to use for adjusting p-values, see ?p.adjust",
    "pval",      "p",  1,   "numeric",      "pval cutoff [0.05]",
    "qval",      "q",  1,   "numeric",      "qval cutoff [0.2]",
    "simplify",  "s",  0,   "logical",      "using simplify() function on GO enrichment",
    "type",      "t",  1,   "character",    "Type of analysis: GO|KEGG|GSEA|wiki|reactome|disease"
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
  cat(paste("Load package: ", name, 'version ', packageVersion(name), ".\n"))
}

WriteText <- function(file, enrichResult, bool=FALSE, type=NULL ) {
  # output text file
  if(isTRUE(bool)) {
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
}

GetEnrichCount <- function(geneEnrich){
  c <- geneEnrich
  count <- nrow(c[c$pvalue < geneEnrich@pvalueCutoff & c$p.adjust < geneEnrich@pvalueCutoff & c$qvalue < geneEnrich@qvalueCutoff,])
  return(count)
}

TryUpsetPlot <- function(geneEnrich, maxCat) {
  tryCatch(
    expr = {
      p <- upsetplot(geneEnrich, n = maxCat)
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )
}

PlotSave <- function(filename, object, h=9, w=9) {
  tryCatch(
    expr = {
      ggplot2::ggsave(filename,
         height = h,
         width = w,
         plot = object,
         dpi=300
      )
      TRUE
    },
    error = function(e) {
      file.remove(filename)
      FALSE
    }
  )
}

PlotDotPlotWriteText <- function(geneEnrich, maxCat, pdfName, txtFile) {
  if (! is.null(geneEnrich)) {
    p <- enrichplot::dotplot(geneEnrich, showCategory = maxCat, font.size = 10)
    plotBool <- PlotSave(pdfName, p, 9, 12)
    WriteText(txtFile, geneEnrich, plotBool)
  }
}

# parse arguments
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

byVector <- c('fgsea', 'DOSE')
if ( isFALSE(args$by %in% byVector) ) {
  ShowHelp(args$by, '-b|--by')
}

methodVector <- c('GO', 'KEGG', 'GSEA', 'wiki', 'reactome', 'disease')
if ( isFALSE(args$type %in% methodVector) ) {
  ShowHelp(args$type, '-t|--type')
}

organismVector <- c('human', 'mouse')
if ( isFALSE(args$organism %in% organismVector) ) {
  ShowHelp(args$organism, '-g|--organism')
}

categoryVector <- c('H', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8')
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
  args$max <- 20
}

if ( is.null(args$minGSSize) ) {
  args$minGSSize <- 10
}

if ( is.null(args$maxGSSize) ) {
  args$maxGSSize <- 500
}

if ( is.null(args$padjust) ) { args$padjust = "BH" }

if ( is.null(args$by) ) {
  args$by <- "fgsea"
}

if ( is.null(args$noCnet) ) {
  args$noCnet <- FALSE
}

if ( is.null(args$noHeatmap) ) {
  args$noHeatmap <- FALSE
}

if ( is.null(args$simplify) ) {
  args$simplify <- FALSE
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

if (args$type == 'wiki' || args$type == 'GSEA') {
  if (args$organism == "human"){
    organism = "Homo sapiens"
  }else{
    organism = "Mus musculus"
  }
}

# create output directory
dir.create(args$output, showWarnings = FALSE, recursive = TRUE)

## print running information
if (args$type == 'GO') {
  cat(paste("\n\nRunning gene enrichment", args$prefix, "GO", args$gotype, ".\n", sep=", "))
}else if (args$type == 'GSEA') {
  cat(paste("\n\nRunning gene enrichment", args$prefix, "GSEA", args$category, ".\n", sep=", "))
}else{
  cat(paste("\n\nRunning gene enrichment", args$prefix, "KEGG", ".\n", sep=", "))
}

## main program
LoadPacakge('ggplot2')
LoadPacakge('clusterProfiler')
LoadPacakge('ChIPseeker')
LoadPacakge(dbName)
LoadPacakge("enrichplot")
LoadPacakge("GOSemSim")
LoadPacakge("DOSE")

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
  geneBitr <- bitr(geneMtx$id, fromType = args$item, toType = "ENTREZID", OrgDb = get(dbName))
  names(geneBitr)[1] <- paste("oldItem")
  # merge geneBitr and geneMtx
  combineData <- merge(geneBitr, geneMtx, by.x='oldItem', by.y="id")
  geneBitr <- bitr(combineData$ENTREZID, fromType = 'ENTREZID', toType = "SYMBOL", OrgDb = get(dbName))
  ## becomes: ENTREZID SYMBOL oldItem value
  combineData <- merge(geneBitr, combineData, by.x="ENTREZID", by.y="ENTREZID")
}else{
  ## filter non-valided ids
  geneBitr <- bitr(geneMtx$id, fromType = args$item, toType = "SYMBOL", OrgDb = get(dbName))
  combineData <- merge(geneBitr, geneMtx, by.x="ENTREZID", by.y="id")
  combineData <- combineData[c("ENTREZID", "SYMBOL", "value")]
}

# feature 1: numeric vector
geneList <- combineData$value
# feature 2: named vector
names(geneList) <- as.character(combineData$ENTREZID)
# feature 3: decreasing orde
geneList <-sort(geneList, decreasing = TRUE)
# get gene list
gene <- names(geneList)

# pdf name
dotplotPdf = file.path(args$output, paste(args$prefix, ".dotplot.pdf", sep=""))
upsetPdf = file.path(args$output, paste(args$prefix, ".upset.pdf", sep=""))
ridgeplotPdf = file.path(args$output, paste(args$prefix, ".ridgeplot.pdf", sep=""))
cnetplotPdf = file.path(args$output, paste(args$prefix, ".cnetplot.pdf", sep=""))
heatplotPdf = file.path(args$output, paste(args$prefix, ".heatplot.pdf", sep=""))
head(geneList)
#text file
txtFile = file.path(args$output, paste(args$prefix, ".geneEnrich.txt", sep=""))
# enrichment analysis
if (args$type == "GO") {
  ## run enrichGO
  geneEnrich <- enrichGO(
    gene          = gene,
    OrgDb         = get(dbName),
    ont           = args$gotype,
    keyType       = "ENTREZID",
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval,
    readable      = TRUE
  )
  ## print to pdf
  # remove redundent GO terms
  if (isTRUE(args$simplify)) {
    geneEnrich <- simplify(geneEnrich)
  }
  ## dotplot plot
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
  ## cnet plot
  if (isFALSE(args$noCnet)) {
    p <- cnetplot(geneEnrich, foldChange=geneList, categorySize="pvalue", showCategory=args$max)
    plotBool <- PlotSave(cnetplotPdf, p, 25, 25)
  }
  ## upset plot
  plotBool <- TryUpsetPlot(geneEnrich, args$max)
  if (isTRUE(plotBool)) {
    p <- upsetplot(geneEnrich, n = args$max)
    plotBool <- PlotSave(upsetPdf, p, 11, 22)
  }
  ### emapplot
  enrichItemCount <- GetEnrichCount(geneEnrich)
  geneEnrich <- pairwise_termsim(geneEnrich, showCategory = enrichItemCount)
  ## for group
  emapplotPdf = file.path(args$output, paste(args$prefix, ".emapplot.group.pdf", sep=""))
  p <- emapplot(geneEnrich, showCategory=enrichItemCount, node_label="group", clusterFunction = stats::kmeans, cex_line=0.6)
  plotBool <- PlotSave(emapplotPdf, p, 12, 12)
  ## for category
  emapplotPdf = file.path(args$output, paste(args$prefix, ".emapplot.category.pdf", sep=""))
  p <- emapplot(geneEnrich, showCategory=enrichItemCount, node_label="category", clusterFunction = stats::kmeans, cex_line=0.6, cex_label_category=0.6)
  plotBool <- PlotSave(emapplotPdf, p, 12, 12)
  ### heatmap
  if (isFALSE(args$noHeatmap)) {
    p <- heatplot(geneEnrich, foldChange=geneList, showCategory=args$max)
    plotBool <- PlotSave(heatplotPdf, p, 25, 25)
  }
}else if (args$type == "KEGG") {
  ## run enrichKEGG
  geneEnrich <- enrichKEGG(
    gene          = gene,
    keyType       = 'kegg',
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    organism      = organism,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval
  )
  geneEnrich <- DOSE::setReadable(geneEnrich, get(dbName), keyType = "ENTREZID")
  ## print to pdf
  dotplotPdf = file.path(args$output, paste(args$prefix, ".dotplot.KEGG.pdf", sep=""))
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
  ## upset plot
  plotBool <- TryUpsetPlot(geneEnrich, args$max)
  if (isTRUE(plotBool)) {
    upsetPdf = file.path(args$output, paste(args$prefix, ".upset.KEGG.pdf", sep=""))
    p <- upsetplot(geneEnrich, n = args$max)
    plotBool <- PlotSave(upsetPdf, p, 11, 22)
  }
  ## run enrichMKEGG
  geneEnrich <- enrichMKEGG(
    gene          = gene,
    keyType       = 'kegg',
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    organism      = organism,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval
  )
  geneEnrich <- DOSE::setReadable(geneEnrich, get(dbName), keyType = "ENTREZID")
  ## dotplot
  dotplotPdf = file.path(args$output, paste(args$prefix, ".dotplot.MKEGG.pdf", sep=""))
  txtFile = file.path(args$output, paste(args$prefix, ".geneEnrich.MKEGG.txt", sep=""))
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
  ## upset plot
  plotBool <- TryUpsetPlot(geneEnrich, args$max)
  if (isTRUE(plotBool)) {
    upsetPdf = file.path(args$output, paste(args$prefix, ".upset.MKEGG.pdf", sep=""))
    p <- upsetplot(geneEnrich, n = args$max)
    plotBool <- PlotSave(upsetPdf, p, 11, 22)
  }
}else if(args$type == "wiki") {
  geneEnrich <- enrichWP(
    gene          = gene,
    organism      = organism,
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval
  )
  geneEnrich <- DOSE::setReadable(geneEnrich, get(dbName), keyType = "ENTREZID")
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
}else if(args$type == "reactome") {
  LoadPacakge("ReactomePA")
  geneEnrich <- enrichPathway(
    gene          = gene,
    organism      = args$organism,
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval,
    readable      = TRUE
  )
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
}else if(args$type == "disease") {
  ## Disease Ontology
  geneEnrich <- enrichDO(
    gene          = gene,
    ont           = "DO",
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval,
    readable      = TRUE
  )
  dotplotPdf = file.path(args$output, paste(args$prefix, ".dotplot.DO.pdf", sep=""))
  txtFile = file.path(args$output, paste(args$prefix, ".enrichedDO.txt", sep=""))
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
  ## network of cancer gene
  geneEnrich <- enrichNCG(
    gene          = gene,
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval,
    readable      = TRUE
  )
  dotplotPdf = file.path(args$output, paste(args$prefix, ".dotplot.NCG.pdf", sep=""))
  txtFile = file.path(args$output, paste(args$prefix, ".enrichedNCG.txt", sep=""))
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
  ## disease gene network
  geneEnrich <- enrichDGN(
    gene          = gene,
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    qvalueCutoff  = args$qval,
    readable      = TRUE
  )
  dotplotPdf = file.path(args$output, paste(args$prefix, ".dotplot.DGN.pdf", sep=""))
  txtFile = file.path(args$output, paste(args$prefix, ".enrichedDGN.txt", sep=""))
  PlotDotPlotWriteText(geneEnrich, args$max, dotplotPdf, txtFile)
}else if (args$type == "GSEA") {
  LoadPacakge('dplyr')
  LoadPacakge('enrichplot')
  ## get MSigDB dataset from msigdbr
  if (is.null(args$geneset)) {
    LoadPacakge('msigdbr')
    gseaMSigDB <- msigdbr(species = organism, category = args$category) %>% dplyr::select(gs_name, entrez_gene)
  }else{
    gseaMSigDB <- read.table(args$geneset, sep="\t", header=F)
    names(geneMtx) <- c("gs_name", "entrez_gene")
  }
  ## run GSEA
  geneEnrich <- GSEA(
    geneList      = geneList,
    nPerm         = args$nperm,
    eps           = 0.0,
    minGSSize     = args$minGSSize,
    maxGSSize     = args$maxGSSize,
    TERM2GENE     = gseaMSigDB,
    pvalueCutoff  = args$pval,
    pAdjustMethod = args$padjust,
    seed          = TRUE,
    by            = args$by
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
      gseaplotPdf <- file.path(args$output, pdfName)
      p <- gseaplot2(geneEnrich, geneSetID=i, title = title, base_size=10)
      plotBool <- PlotSave(gseaplotPdf, p)
    }
    ## print to txt file
    txtFile = file.path(args$output, paste(args$prefix, ".", args$category, ".geneEnrich.txt", sep=""))
    WriteText(txtFile, geneEnrich, TRUE, type='GSEA')
  }
}

# delete empty Rplots.pdf generated by ggsave()
fn <- "Rplots.pdf"
setwd(getwd())
if (file.exists(fn)) {
  file.remove(fn)
}
