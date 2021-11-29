#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information (https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html)",
    "autoBatch",    "a",   0,  "logical",     "Use 'design' column to remove hidden batch effect (work with empirical)",
    "adjp",         "q",   2,  "numeric",     "adjp cutoff (0.1)",
    "batchMethod",  "b",   1,  "character",   "Remove hidden batch effect (none|empirical|spikeins|upper|median|full)",
    "control",      "c",   1,  "character",   "Name for control design in colData",
    "counts",       "g",   1,  "character",   "Gene counts matrix (gene_id|tx_id in 1st column)",
    "design",       "d",   1,  "character",   "Design for construction of DESeqDataSet (colname in colData)",
    "fitType",      "y",   1,  "character",   "the fitType in DESeq() (parametric[default]localmean|glmGamPoi)",
    "formula",      "k",   1,  "character",   "The design formula in DESeqDataSetFromMatrix() function (eg. genotype + treatment + genotype:treatment), will overwrite --control, --treat, --design",
    "filter",       "f",   2,  "integer",     "Filter out genes less than # counts across 1/2 samples",
    "keepSpike",    "S",   0,  "logical",     "Keep spikeins reads when passing to DESeq()",
    "name",         "n",   1,  "character",   "Name for extract the results formula in results() (eg. genotypewt.treatmentheat for wt (heat vs normal) vs ko (heat vs normal) )",
    "norcounts",    "N",   0,  "logical",     "output normalized Counts",
    "output" ,      "o",   1,  "character",   "Output directory",
    "padjust",      "j",   1,  "character",   "The method to use for adjusting p-values, see ?p.adjust",
    "prefix",       "e",   1,  "character",   "Prefix for output",
    "pval",         "p",   2,  "numeric",     "pval cutoff (0.05)",
    "reduced",      "i",   1,  "character",   "The reduced formula for DESeq(ddds, test=\"LRT\", reduced=)",
    "relevel",      "w",   1,  "character",   "Relevel the sample data for comparison (eg. 'genotype:wt,treatment:normal')",
    "ruvgCount",    "u",   2,  "numeric",     "Counts cutoff for filtering count matrix with --batchMethod empirical (5)",
    "ruvgLogfc",    "l",   2,  "numeric",     "pvalue cutoff for filtering count matrix with --batchMethod empirical (1)",
    "ruvgPval",     "v",   2,  "numeric",     "log2CF cutoff for filtering count matrix with --batchMethod empirical (0.3)",
    "sampleMtx",    "m",   1,  "character",   "Sample relationships matrix",
    "shrink",       "s",   1,  "character",   "Shrinkage method for DE results (none|normal|apeglm[default]|ashr)",
    "spiRegex",     "r",   1,  "character",   "Name pattern of spikeins in gene_id (ERCC-)",
    "treat",        "t",   1,  "character",   "Name for treat design in colData",
    "test",         "T",   1,  "character",   "The test method for p-value (Wald|LRT)"
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

# parsing arguments
args <- getopt(command)

ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$counts, '-g|--counts')
ShowHelp(args$sampleMtx, '-s|--sampleMtx')

if (is.null(args$formula)) {
  ShowHelp(args$design, '-d|--design')
  ShowHelp(args$control, '-c|--control')
  ShowHelp(args$treat, '-t|--treat')
}

if ( is.null(args$test) ) {
  args$test = 'Wald'
}else{
  testVetor <- c('Wald', 'LRT')
  bool <- isFALSE(args$test %in% testVetor)
  ShowHelp(bool, '-T|--test', FALSE, TRUE)
}

if ( is.null(args$batchMethod) ) {
  args$batchMethod = 'none'
}else{
  bmVetor <- c('none', 'empirical', 'spikeins', 'upper')
  bool <- isFALSE(args$batchMethod %in% bmVetor)
  ShowHelp(bool, '-b|--batchMethod', FALSE, TRUE)
}

if ( is.null(args$shrink) ) {
  args$shrink = 'apeglm'
}else{
  shrinkVetor <- c('none', 'normal', 'apeglm', 'ashr')
  bool <- isFALSE(args$shrink %in% shrinkVetor)
  ShowHelp(bool, '-s|--shrink', FALSE, TRUE)
}

if ( is.null(args$fitType) ) {
  args$fitType = 'parametric'
}else{
  fitTypeVetor <- c("parametric", "local", "mean")
  bool <- isFALSE(args$fitType %in% fitTypeVetor)
  ShowHelp(bool, '-y|--fitType', FALSE, TRUE)
}

# default values
if ( is.null(args$filter) ) { args$filter = 1 }
if ( is.null(args$pval) ) { args$pval = 0.05 }
if ( is.null(args$adjp) ) { args$adjp = 0.1 }
if ( is.null(args$spiRegex) ) { args$spiRegex = 'ERCC-' }
if ( is.null(args$keepSpike) ) { args$keepSpike = FALSE }
if ( is.null(args$norcounts) ) { args$norcounts = FALSE }
if ( is.null(args$padjust) ) { args$padjust = "BH" }
if ( is.null(args$prefix) ) { args$prefix = 'result' }
if ( is.null(args$output) ) { args$output = './' }
if ( is.null(args$autoBatch) ) { args$autoBatch = FALSE }
if ( is.null(args$ruvgCount) ) { args$ruvgCount = 5 }
if ( is.null(args$ruvgLogfc) ) { args$ruvgLogfc = 1 }
if ( is.null(args$ruvgPval) ) { args$ruvgPval = 0.3 }
# load DESeq2 and perform Differential Analysis
LoadPacakge('DESeq2')
LoadPacakge('ggplot2')
LoadPacakge("RColorBrewer")
## load arguments
batchMethod <- args$batchMethod
geneCountMtx <- args$counts
sampleMtx <- args$sampleMtx
filter <- args$filter
fitType <- args$fitType
design <- args$design
control <- args$control
treat <- args$treat
test <- args$test
pvalCutoff <- args$pval
padjCuotff <- args$adjp
shrink <- args$shrink
prefix <- args$prefix
output <- args$output
spiRegex <- args$spiRegex
keepSpike <- args$keepSpike
autoBatch <- args$autoBatch

# With the count matrix, cts, and the sample information, colData
cts <- as.matrix(read.csv(geneCountMtx, sep="\t", row.names=1))
cts <-round(cts, 0)
colData <- read.csv(sampleMtx, row.names=1, sep="\t")

if (! is.null(design)) {
  # judge if --control and --treat in colData
  contrast <- c(design, treat, control)
  designs <- unique(as.character(colData[[design]]))
  if ( !(control %in% designs) ) {
    bool <- TRUE
    ShowHelp(bool, '-c|--control', FALSE, TRUE)
  }
  
  if ( !(treat %in% designs) ) {
    bool <- TRUE
    ShowHelp(bool, '-t|--treat', FALSE, TRUE)
  }
}

if (test == "LRT") {
  if (is.null(args$reduced)) {
    bool <- TRUE
    ShowHelp(bool, '-i|--reduced', FALSE, TRUE)
  }
}

# check all sample rownames in geneCountMtx colNames
all(rownames(colData) %in% colnames(cts))
cts <- cts[, rownames(colData)]

# filter out gSpikein_phiX174, ENCODE
genes <- rownames(cts)[grep('_phiX174', rownames(cts), invert=TRUE)]
cts <- cts[genes,]

sampleSize <- nrow(colData)
# pre-filtering, counts > args$filter in at least half of the samples
filtered <- apply(cts, 1, function(x) length(x[x>filter])>=round(sampleSize/2))
cts <- cts[filtered,]

if (!is.null(args$design)) {
  name <- paste( design, treat, 'vs', control, sep="_")
}

if (!is.null(args$name)) {
  name <- args$name
}

# if used spikein, use "RUVSeq" to Estimating the factors of unwanted variation using control genes
if( args$batchMethod == "spikeins" | args$batchMethod == "upper" ){
  ## Removing hidden batch effects using spike-in controls by RUVg
  cat(paste('Removing hidden batch effects with ', args$batchMethod,' (RUVg)!\n', sep=""))
  LoadPacakge('RUVSeq')
  ## seperate to genes and spikes
  genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
  spikes <- rownames(cts)[grep(spiRegex, rownames(cts))]
  set <- newSeqExpressionSet(as.matrix(cts), phenoData = colData)
  if (args$batchMethod == "upper" || args$batchMethod == "median" || args$batchMethod == "full") {
    ## use the betweenLaneNormalization function of EDASeq to normalize the
    ## data using upper-quartile (UQ) normalization
    normalizeSet <- betweenLaneNormalization(set, which=args$batchMethod, offset=TRUE)
  }else{
    ## RUVg: Estimating the factors of unwanted variation using control genes
    normalizeSet <- RUVg(set, spikes, k=1)
  }
  ## re-construct cts, filter out spike-ins if --keepSpike not set
  if(! keepSpike){
    cts <- cts[genes,]
  }
  colData <- pData(normalizeSet)
  if (!is.null(design)) {
    if (batchMethod == "spikeins") {
      designFormula <- as.formula(paste("~ W_1 +", design, sep=" "))
    }else{
      designFormula <- as.formula(paste("~ ", design, sep=" "))
    }
  }
}else{
  ## filter out spike-ins if --keepSpike set
  if(! keepSpike){
    genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
    cts <- cts[genes,]
  }
  if (!is.null(design)) {
    designFormula <- as.formula(paste("~", design, sep=""))
  }
}

if (!is.null(args$formula)) {
  if( args$batchMethod == "spikeins"){
    designFormula <- as.formula(paste("~ W_1 +", args$formula, sep=""))
  }else{
    designFormula <- as.formula(paste("~", args$formula, sep=""))
  }
}

# construct a DESeqDataSet object
if (args$batchMethod == "upper" || args$batchMethod == "median" || args$batchMethod == "full") {
  dds <- DESeqDataSetFromMatrix(countData = counts(normalizeSet),
                                colData = pData(normalizeSet),
                                design = designFormula)
  normFactors <- exp(-1 * offst(normalizeSet))
  normFactors <- normFactors / exp(rowMeans(log(normFactors)))
  normalizationFactors(dds) <- normFactors
}else if (args$batchMethod == "spikeins") {
  dds <- DESeqDataSetFromMatrix(countData = counts(normalizeSet),
                                colData = colData,
                                design = designFormula)
}else{
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = colData,
                                design = designFormula)
}

## relevel the dds
if (! is.null(args$relevel)) {
  ## "phenotype:sgNS,fraction:nuc"
  relevelVector <- unlist(strsplit(args$relevel, ","))
  ##[1] "phenotype:sgNS" "fraction:nuc"
  for (eachRelevel in relevelVector) {
    eachRelevelVector <- unlist(strsplit(eachRelevel, ":"))
    relevelName <- eachRelevelVector[1]
    relevelVal <- eachRelevelVector[2]
    dds[[relevelName]] <- relevel(dds[[relevelName]], relevelVal)
  }
}

if (args$batchMethod == "upper" || args$batchMethod == "median" || args$batchMethod == "full") {
  nothing <- 0
}else{
  featureData <- data.frame(gene=rownames(cts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
}

if (! is.null(args$reduced) && test != "Wald") {
  reduced <- as.formula(paste("~", args$reduced, sep=""))
}

# runing DESeq
if (test == "LRT") {
  dds <- DESeq(dds, fitType=fitType, test=test, reduced = reduced)
}else{
  dds <- DESeq(dds, fitType=fitType, test=test)
}

# Removing hidden batch effects using RUVg, --batchMethod: empirical
if (args$batchMethod == 'empirical') {
  cat('Removing hidden batch effects using empirical control genes.\n')
  cat('Using empirical control genes by looking at the genes that do not have a small p-value\n')
  LoadPacakge('RUVSeq')
  res <- results(dds, pAdjustMethod=args$padjust, contrast=c(design, treat, control))
  ## removing hidden batch effect
  set <- newSeqExpressionSet(counts(dds))
  idx <- rowSums(counts(set) > args$ruvgCount) == sampleSize
  set <- set[idx, ]
  set <- betweenLaneNormalization(set, which="upper")
  if (autoBatch) {
    empirical <- rownames(cts)
    tDesigns <- designs[designs != control]
    for (tDesign in tDesigns) {
      tmpRes <- results(dds, pAdjustMethod=args$padjust, contrast=c(design, tDesign, control))
      tmpSet <- newSeqExpressionSet(counts(dds))
      tmpIdx <- rowSums(counts(tmpSet) > args$ruvgCount) == sampleSize
      tmpSet <- tmpSet[tmpIdx, ]
      tmpSet <- betweenLaneNormalization(tmpSet, which="upper")
      tmpNotSig <- rownames(tmpRes)[which(tmpRes$pvalue > args$ruvgPval & abs(tmpRes$log2FoldChange) < args$ruvgLogfc)]
      tmpEmpirical <- rownames(tmpSet)[ rownames(tmpSet) %in% tmpNotSig ]
      empirical <- intersect(empirical, tmpEmpirical)
    }
  }else{
    notSig <- rownames(tmpRes)[which(tmpRes$pvalue > args$ruvgPval & abs(tmpRes$log2FoldChange) < args$ruvgLogfc)]
    empirical <- rownames(set)[ rownames(set) %in% notSig ]
  }
  cat('empirical genes:')
  print(length(empirical))
  set <- RUVg(set, empirical, k=2)
  ## assign W_1 and W_2 to dd
  dds$W1 <- set$W_1
  dds$W2 <- set$W_2
  ## re-design factors
  if (!is.null(args$design)) {
    design(dds) <- as.formula(paste("~ W1 + W2 +", design, sep=" "))
  }
  if (!is.null(args$formula)) {
    design(dds) <- as.formula(paste("~ W1 + W2 +", args$formula, sep=" "))
  }
  if (test == "LRT") {
    dds <- DESeq(dds, fitType=fitType, test=test, reduced = reduced)
  }else{
    dds <- DESeq(dds, fitType=fitType, test=test)
  }
}

resultsNameVector <- resultsNames(dds)
resultsNameVector
if ( name %in% resultsNameVector ) {
  name
}else{
  contrast
}

if ( shrink != 'none' ) {
  if ( name %in% resultsNameVector ) {
    res <- lfcShrink(dds, coef=name, type=shrink)
  }else{
    if (shrink == 'apeglm') {
      print('change shrink method to "ashr"!')
      shrink <- 'ashr'
    }
    res <- lfcShrink(dds, contrast=contrast, type=shrink)
  }
}else{
  if ( name %in% resultsNameVector ) {
    res <- results( dds, pAdjustMethod=args$padjust, name=name )
  }else{
    res <- results( dds, pAdjustMethod=args$padjust, contrast=contrast )
  }
}

# to avoid NA adjusted p-values
res$padj <- p.adjust(res$pvalue, method="BH")

# plot MA plot
maPlotPdf <- file.path(output, paste(prefix, ".MA.pdf", sep=""))
pdf(maPlotPdf, paper='a4r', height=0)
DESeq2::plotMA(res, ylim=c(-6,6), alpha=padjCuotff, cex=0.6)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
garbage <- dev.off()

# plot the histogram of the p values
histoPlotPdf <- file.path(output, paste(prefix, ".pvalue.histogram.pdf", sep=""))
pdf(histoPlotPdf, paper='a4r', height=0)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", 
     xlab = 'p-value', main = 'Histogram of p-value (baseMean > 1)')
garbage <- dev.off()

# The ratio of small p values for genes binned by mean normalized count.
if (test == "Wald") {
  resLFC1 <- results(dds, pAdjustMethod=args$padjust, lfcThreshold=1)
  qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
  bins <- cut(resLFC1$baseMean, qs)
  levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
  fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
                            mean(p < .05, na.rm = TRUE))
  barPlotPdf <- file.path(output, paste(prefix, ".pvalueNorCounts.bar.pdf", sep=""))
  pdf(barPlotPdf, paper='a4r', height=0)
  barplot(fractionSig, xlab = "Mean normalized count",
                       ylab = "Fraction of small p values")
  garbage <- dev.off()
}

# output result
res <- cbind(uniqId = rownames(res), res)
resultFile <- file.path(output, paste(prefix, ".DESeq2.txt", sep=""))
output.file <- file(resultFile, "wb")
write.table(as.data.frame(res), sep="\t", eol = "\n", 
            quote = FALSE, row.names=FALSE, file=output.file)
close(output.file)

# significant output result
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < padjCuotff & pvalue < pvalCutoff)
resultFile <- file.path(output, paste(prefix, ".sig.DESeq2.txt", sep=""))
output.file <- file(resultFile, "wb")
write.table(as.data.frame(resSig), sep="\t", eol = "\n", 
            quote = FALSE, row.names=FALSE, file=output.file)
close(output.file)

# normalzed counts
if (isTRUE(args$norcounts)) {
  resultFile <- file.path(output, paste(prefix, ".normalized.counts.txt", sep=""))
  output.file <- file(resultFile, "wb")
  normalzedCounts <- counts(dds, normalized=TRUE)
  normalzedCounts <- cbind(uniqId = rownames(normalzedCounts), normalzedCounts)
  write.table(as.data.frame(normalzedCounts), sep="\t", eol = "\n", 
              quote = FALSE, row.names=FALSE, file=output.file)
  close(output.file)
}
