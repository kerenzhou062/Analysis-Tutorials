#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "adjp",         "q",   2,  "numeric",     "adjp cutoff (0.1)",
    "control",      "cn",  1,  "character",   "Name for control design in colData",
    "counts",       "ct",  1,  "character",   "Gene counts matrix",
    "filter",       "f",   2,  "integer",     "Filter out genes less than # counts across all samples",
    "design",       "de",  1,  "character",   "Design for construction of DESeqDataSet (colname in colData)",
    "pval",         "p",   2,  "numeric",     "pval cutoff (0.05)",
    "keepCol",      "kc",  2,  "character",   "Only 'excol' in samplemtx keeped before passing to DESeq()",
    "keepRow",      "kr",  2,  "character",   "Only 'exrow' from 'excol' in samplemtx keeped before passing to DESeq()",
    "removeSp",     "rs",  1,  "logical",     "Remove spikeins before reads passing to DESeq()",
    "sampleMtx",    "sm",  1,  "character",   "Sample relationships matrix",
    "shrink",       "sr",  1,  "character",   "Shrinkage method for DE results (none|normal|apeglm[default]|ashr)",
    "spikein",      "sp",  0,  "logical",     "EstimateSizeFactors with spikeins",
    "spiRegex",     "pa",  1,  "character",   "Name pattern of spikeins in gene_id (Spikein-ERCC-|ERCC-)",
    "treat",        "tr",  1,  "character",   "Name for treat design in colData",
    "test",         "t",   1,  "character",   "test method for p-value (Wald|LRT)",
    "prefix",       "pr",  1,  "character",   "Prefix for output",
    "output" ,      "o",   1,  "character",   "Output directory"
), byrow=TRUE, ncol=5)

## parsing arguments
args <- getopt(command)

if(!is.null(args$help)){
  cat(paste(getopt(command, usage = T),"\n"))
  q()
}

if(is.null(args$design) || is.null(args$control) || is.null(args$treat)){
  cat(paste(getopt(command, usage = T),"\n"))
  q()
}

if(is.null(args$test)){
  cat(paste(getopt(command, usage = T),"\n"))
  q()
}

if ( is.null(args$test) ) {
  args$test = 'Wald'
}else{
  testVetor <- c('Wald', 'LRT')
  if (isFALSE(args$test %in% testVetor)) {
    cat("None valid -t|--test!\n")
    cat(paste(getopt(command, usage = T),"\n"))
    q()
  }
}

if ( is.null(args$pval) ) { args$pval = 0.05 }
if ( is.null(args$adjp) ) { args$adjp = 0.1 }
if ( is.null(args$spiregex) ) { args$spiregex = 'ERCC-' }

if ( is.null(args$shrink) ) {
  args$shrink = 'apeglm'
}else{
  shrinkVetor <- c('none', 'normal', 'apeglm', 'ashr')
  if (isFALSE(args$shrink %in% shrinkVetor)) {
    cat("None valid -sr|--shrink!\n")
    cat(paste(getopt(command, usage = T),"\n"))
    q()
  }
}

## load DESeq2 and perform Differential Analysis
suppressMessages(library('DESeq2'))
suppressMessages(library('ggplot2'))
## load arguments
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

# With the count matrix, cts, and the sample information, colData
cts <- as.matrix(read.csv(geneCountMtx, sep="\t", row.names="gene_id"))
cts <-round(cts,0)
colData <- read.csv(sampleMtx, row.names=1, sep="\t")
#colData <- colData[,c("condition","type")]

## check all sample rownames in geneCountMtx colNames
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

# if used spikein, use "RUVSeq" to Estimating the factors of unwanted variation using control genes
if( !is.null(args$spikein) ){
  suppressMessages(library('RUVSeq'))
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
  ## re-construct cts, filter out spike-ins if --removeSp set
  if(!is.null(args$removeSp)){
    cts <- cts[genes,]
  }
  colData <- pData(spikeNorSet)
  designFormula <- as.formula(paste("~", "W_1", "+", design, sep=" "))
}else{
  ## filter out spike-ins if --removeSp set
  if(!is.null(args$removeSp)){
    genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
    cts <- cts[genes,]
  }
  designFormula <- as.formula(paste("~", design, sep=" "))
}

## keep colData controled by --keepRow and --keepCol
if (!is.null(args$keepCol) & !is.null(args$keepRow)) {
  colData <- colData[colData[[args$keepCol]] == args$keepRow, ]
}

# reconstruct cts dataframe, remove unwanted samples
all(rownames(colData) %in% colnames(cts))
cts <- cts[, rownames(colData)]

# construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = designFormula)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# Note on factor levels
dds[[design]] <- factor(dds[[design]], levels = c(control, treat))
# dds$condition <- relevel(dds$condition, ref = "untreated")
dds <- DESeq(dds, test=test)
res <- results(dds, contrast=c(design, treat, control))

if ( shrink != 'none' ) {
  coefName = paste(design, treat, 'vs', control, sep="_")
  res <- lfcShrink(dds, coef=coefName, type=shrink)
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
