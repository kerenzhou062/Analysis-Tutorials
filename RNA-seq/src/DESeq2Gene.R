#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "adjp",         "q",   2,  "numeric",     "adjp cutoff (0.1)",
    "batchMethod",  "b",   1,  "character",   "Remove hidden batch effect (none|RUVg|spikeins)",
    "control",      "c",   1,  "character",   "Name for control design in colData",
    "counts",       "g",   1,  "character",   "Gene counts matrix",
    "design",       "d",   1,  "character",   "Design for construction of DESeqDataSet (colname in colData)",
    "filter",       "f",   2,  "integer",     "Filter out genes less than # counts across all samples",
    "keepSpike",    "S",   0,  "logical",     "Keep spikeins reads when passing to DESeq()",
    "output" ,      "o",   1,  "character",   "Output directory",
    "prefix",       "e",   1,  "character",   "Prefix for output",
    "pval",         "p",   2,  "numeric",     "pval cutoff (0.05)",
    "ruvgCount",    "u",   2,  "integer",     "Counts cutoff for filtering count matrix with --batchMethod RUVg (5)",
    "selectCol",    "C",   2,  "character",   "Only 'selectCol' in samplemtx was selected as contrast in results()",
    "selectRow",    "R",   2,  "character",   "Only 'selectRow' from 'selectCol' was selected as contrast in results()",
    "sampleMtx",    "m",   1,  "character",   "Sample relationships matrix",
    "shrink",       "s",   1,  "character",   "Shrinkage method for DE results (none|normal|apeglm[default]|ashr)",
    "spiRegex",     "r",   1,  "character",   "Name pattern of spikeins in gene_id (Spikein-ERCC-|ERCC-)",
    "treat",        "t",   1,  "character",   "Name for treat design in colData",
    "test",         "T",   1,  "character",   "test method for p-value (Wald|LRT)"
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
  cat(paste("Load package: ", name, ".\n"))
}

# parsing arguments
args <- getopt(command)

ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$counts, '-g|--counts')
ShowHelp(args$sampleMtx, '-s|--sampleMtx')
ShowHelp(args$design, '-d|--design')
ShowHelp(args$control, '-c|--control')
ShowHelp(args$treat, '-t|--treat')

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
  bmVetor <- c('none', 'RUVg', 'spikeins')
  bool <- isFALSE(args$batchMethod %in% bmVetor)
  ShowHelp(bool, '-b|--batchMethod', FALSE, TRUE)
}

if ( is.null(args$shrink) ) {
  args$shrink = 'apeglm'
}else{
  shrinkVetor <- c('none', 'normal', 'apeglm', 'ashr')
  bool <- isFALSE(args$shrink %in% shrinkVetor))
  ShowHelp(bool, '-s|--shrink', FALSE, TRUE)
}

if (!is.null(args$selectCol)) {
  bool <- is.null(args$selectRow)
  ShowHelp(bool, '-C|--selectCol & -R|--selectRow', FALSE, TRUE)
}

if (!is.null(args$selectRow)) {
  bool <- is.null(args$selectCol)
  ShowHelp(bool, '-C|--selectCol & -R|--selectRow', FALSE, TRUE)
}

# default values
if ( is.null(args$pval) ) { args$pval = 0.05 }
if ( is.null(args$adjp) ) { args$adjp = 0.1 }
if ( is.null(args$spiRegex) ) { args$spiRegex = 'ERCC-' }
if ( is.null(args$keepSpike) ) { args$keepSpike = FALSE }
if ( is.null(args$prefix) ) { args$prefix = 'result' }
if ( is.null(args$output) ) { args$output = './' }
if ( is.null(args$ruvgCount) ) { args$ruvgCount = 5 }

# load DESeq2 and perform Differential Analysis
suppressMessages(library('DESeq2'))
suppressMessages(library('ggplot2'))
LoadPacakge('DESeq2')
LoadPacakge('ggplot2')
## load arguments
batchMethod <- args$batchMethod
geneCountMtx <- args$counts
sampleMtx <- args$sampleMtx
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

# With the count matrix, cts, and the sample information, colData
cts <- as.matrix(read.csv(geneCountMtx, sep="\t", row.names="gene_id"))
cts <-round(cts,0)
colData <- read.csv(sampleMtx, row.names=1, sep="\t")

# check all sample rownames in geneCountMtx colNames
all(rownames(colData) %in% colnames(cts))
cts <- cts[, rownames(colData)]

# filter out gSpikein_phiX174, ENCODE
genes <- rownames(cts)[grep('_phiX174', rownames(cts), invert=TRUE)]
cts <- cts[genes,]

sampleSize <- nrow(colData)
# pre-filtering, counts > args$filter in at least half of the samples
if(!is.null(args$filter)){
  filter <- apply(cts, 1, function(x) length(x[x>args$filter])>=round(sampleSize/2))
  cts <- cts[filter,]
}

if (!is.null(args$selectCol)) {
  name <- paste( design,control, '.', selectCol, selectRow, sep="")
}else{
  name <- paste( design, treat, 'vs', control, sep="_")
}

# if used spikein, use "RUVSeq" to Estimating the factors of unwanted variation using control genes
if( args$batchMethod == "spikeins" ){
  ## Removing hidden batch effects using spike-in controls by RUVg
  cat('Removing hidden batch effects with spike-ins (RUVg)!\n')
  suppressMessages(library('RUVSeq'))
  LoadPacakge('RUVSeq')
  ## seperate to genes and spikes
  genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
  spikes <- rownames(cts)[grep(spiRegex, rownames(cts))]
  set <- newSeqExpressionSet(as.matrix(cts), phenoData = colData)
  ## use the betweenLaneNormalization function of EDASeq to normalize the
  ## data using upper-quartile (UQ) normalization
  set <- betweenLaneNormalization(set, which="upper")
  ## RUVg: Estimating the factors of unwanted variation using control genes
  spikeNorSet <- RUVg(set, spikes, k=1)
  ## pass spikeNorSet to DESeq2
  ## re-construct cts, filter out spike-ins if --keepSpike set
  if(! keepSpike){
    cts <- cts[genes,]
  }
  colData <- pData(spikeNorSet)

  if (!is.null(args$selectCol)) {
    designFormula <- as.formula(paste("~ W1 +", design, '+', design, ':', selectCol, sep=" "))
  }else{
    designFormula <- as.formula(paste("~ W1 +", design, sep=" "))
  }
}else{
  ## filter out spike-ins if --keepSpike set
  if(! keepSpike){
    genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
    cts <- cts[genes,]
  }
  designFormula <- as.formula(paste("~", design, sep=" "))
}

# construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = designFormula)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
# Note on factor levels
dds[[design]] <- factor(dds[[design]], levels = c(control, treat))
dds <- DESeq(dds, test=test)

# Removing hidden batch effects using RUVg, --batchMethod: RUVg
if (args$batchMethod == 'RUVg') {
  cat('Removing hidden batch effects using RUVg.\n')
  cat('Using empirical control genes by looking at the genes that do not have a small p-value\n')
  suppressMessages(library('RUVSeq'))
  LoadPacakge('RUVSeq')

  featureData <- data.frame(gene=rownames(cts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  # perform DE analysis before passing to RUVg
  dds[[design]] <- factor(dds[[design]], levels = c(control, treat))
  dds <- DESeq(dds, test=test)
  res <- results(dds, contrast=c(design, treat, control))
  # removing hidden batch effect
  set <- newSeqExpressionSet(counts(dds), phenoData = colData)
  idx <- rowSums(counts(set) > args$ruvgCount) >= round(sampleSize/2)
  set <- set[idx, ]
  set <- betweenLaneNormalization(set, which="upper")
  notSig <- rownames(res)[which(res$pvalue > .1)]
  empirical <- rownames(set)[ rownames(set) %in% notSig ]
  set <- RUVg(set, empirical, k=2)
  ## assign W_1 and W_2 to dd
  dds$W1 <- set$W_1
  dds$W2 <- set$W_2
  ## re-design factors
  if (!is.null(args$selectCol)) {
    design(dds) <- as.formula(paste("~ W1 + W2 +", design, '+', design, ':', selectCol, sep=" "))
  }else{
    design(dds) <- as.formula(paste("~ W1 + W2 +", design, sep=" "))
  }
  dds <- DESeq(dds, test=test)
}

if ( shrink != 'none' ) {
  res <- lfcShrink(dds, coef=name, type=shrink)
}else{
  res <- results( dds, name=name )
}

# plot MA plot
maPlotPdf <- file.path(output, paste(prefix, ".MA.pdf", sep=""))
pdf(maPlotPdf, paper='a4r', height=0)
plotMA(res, ylim=c(-6,6), alpha=padjCuotff, cex=0.8)
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
resLFC1 <- results(dds, lfcThreshold=1)
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

# output result
resultFile <- file.path(output, paste(prefix, ".DE.txt", sep=""))
output.file <- file(resultFile, "wb")
write.table(as.data.frame(res), sep="\t", eol = "\n", 
            quote = FALSE, row.names=TRUE, file=output.file)
close(output.file)

# significant output result
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < padjCuotff)
resSig <- subset(resSig, pvalue < pvalCutoff)
resultFile <- file.path(output, paste(prefix, ".DE.sig.txt", sep=""))
output.file <- file(resultFile, "wb")
write.table(as.data.frame(resSig), sep="\t", eol = "\n", 
            quote = FALSE, row.names=TRUE, file=output.file)
close(output.file)
