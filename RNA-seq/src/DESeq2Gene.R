#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "show help information",
    "filter",       "f",   0,  "logical",     "filter out genes with 0 counts across all samples",
    "spikein",      "sp",  0,  "logical",     "estimateSizeFactors with spikeins",
    "spiregex",     "pa",  1,  "character",   "name pattern of spikeins in gene_id (Spikein-ERCC-|ERCC-)",
    "counts",       "ct",  1,  "character",   "Gene counts matrix",
    "sampleMtx",    "mt",  1,  "character",   "Sample relationships matrix",
    "design",       "de",  1,  "character",   "design for construction of DESeqDataSet (colname in colData)",
    "pval",         "p",   2,  "numeric",     "pval cutoff (0.05)",
    "adjp",         "q",   2,  "numeric",     "adjp cutoff (0.1)",
    "shrink",       "sh",  0,  "character",   "shrinkage method for DE results (none|normal|apeglm[default]|ashr)",
    "control",      "c",   1,  "character",   "name for control design in colData",
    "treat",        "t",   1,  "character",   "name for treat design in colData",
    "prefix",       "pr",  1,  "character",   "prefix for output",
    "output" ,      "o",   1,  "character",   "output directory"
), byrow=TRUE, ncol=5)

## parsing arguments
args <- getopt(command)

if(!is.null(args$help)){
  cat(paste(getopt(command, usage = T),"\n"))
  q()
}

if ( is.null(args$pval) ) { args$pval = 0.05 }
if ( is.null(args$adjp) ) { args$adjp = 0.1 }
if ( is.null(args$spiregex) ) { args$spiregex = 'ERCC-' }

if ( is.null(args$shrink) ) {
  args$shrink = 'apeglm'
}else{
  shrinkVetor <- c('none', 'normal', 'apeglm', 'ashr')
  if (isFALSE(args$shrink %in% shrinkVetor)) {
    cat("None valid -sh|--shrink!\n")
    cat(paste(getopt(command, usage = T),"\n"))
    q()
  }
}

## load DESeq2 and perform Differential Analysis
suppressMessages(library('DESeq2'))
## load arguments
geneCountMtx <- args$counts
sampleMtx <- args$sampleMtx
design <- args$design
control <- args$control
treat <- args$treat
pvalCutoff <- args$pval
padjCuotff <- args$adjp
prefix <- args$prefix
output <- args$output
spiRegex <- args$spiregex

# With the count matrix, cts, and the sample information, colData
cts <- as.matrix(read.csv(geneCountMtx, sep="\t", row.names="gene_id"))
colData <- read.csv(sampleMtx, row.names=1, sep="\t")
#colData <- colData[,c("condition","type")]

## check all sample rownames in geneCountMtx colNames
all(rownames(colData) %in% colnames(cts))
cts <- cts[, rownames(colData)]

# if used spikein, use "RUVSeq" to Estimating the factors of unwanted variation using control genes
if( !is.null(args$spikein) ){
  suppressMessages(library('RUVSeq'))
  library('RColorBrewer')
  colors <- brewer.pal(3, "Set2")

  ## seperate to genes and spikes
  genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
  spikes <- rownames(cts)[grep(spiRegex, rownames(cts))]
  spif <- as.factor(colData[[design]])
  set <- newSeqExpressionSet(as.matrix(cts),
         phenoData = data.frame(spif, row.names=colnames(cts)))
  ## use the betweenLaneNormalization function of EDASeq to normalize the
  ## data using upper-quartile (UQ) normalization
  set <- betweenLaneNormalization(set, which="upper")
  ## plot RLE and PCA before RUVg
  rlePlotPdf <- file.path(output, paste(prefix, ".beforeSpikein.RLE.pdf", sep="."))
  pdf(rlePlotPdf, paper='a4r', height=0)
  plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[spif])
  garbage <- dev.off()
  pcaPlotPdf <- file.path(output, paste(prefix, ".beforeSpikein.PCA.pdf", sep="."))
  pdf(pcaPlotPdf, paper='a4r', height=0)
  plotPCA(set, col=colors[spif], cex=1.2)
  garbage <- dev.off()
  ## RUVg: Estimating the factors of unwanted variation using control genes
  spikeNorSet <- RUVg(set, spikes, k=1)
  ## plot RLE and PCA after RUVg
  rlePlotPdf <- file.path(output, paste(prefix, ".afterSpikein.RLE.pdf", sep="."))
  pdf(rlePlotPdf, paper='a4r', height=0)
  plotRLE(spikeNorSet, outline=FALSE, ylim=c(-4, 4), col=colors[spif])
  garbage <- dev.off()
  pcaPlotPdf <- file.path(output, paste(prefix, ".afterSpikein.PCA.pdf", sep="."))
  pdf(pcaPlotPdf, paper='a4r', height=0)
  plotPCA(spikeNorSet, col=colors[spif], cex=1.2)
  garbage <- dev.off()
  ## pass spikeNorSet to DESeq2
  ### rename "spif" -> design
  spikeColData <- pData(spikeNorSet)
  colnames(spikeColData)[1] <- design
  ## re-construct cts and colData
  cts <- cts[genes,]
  colData <- spikeColData
  designFormula <- as.formula(paste("~", "W_1", "+", design, sep=" "))
}else{
  genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
  cts <- cts[genes,]
  designFormula <- as.formula(paste("~", design, sep=" "))
}

# construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = designFormula)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# pre-filtering, counts > 0
if(!is.null(args$filter)){
  keep <- rowSums(counts(dds)) > 0
  dds <- dds[keep,]
}

# Note on factor levels
dds[[design]] <- factor(dds[[design]], levels = c(control, treat))
# dds$condition <- relevel(dds$condition, ref = "untreated")
dds <- DESeq(dds)
res <- results(dds)
coefName = paste(design, control, 'vs', treat, sep="_")
res <- results(dds, name=coefName)
res <- results(dds, contrast=c(design, treat, control))

if ( args$shrink != 'none' ) {
  res <- lfcShrink(dds, coef=coefName, type=args$shrink)
}

# plot MA plot
maPlotPdf <- file.path(output, paste(prefix, "MA.pdf", sep="."))
pdf(maPlotPdf, paper='a4r', height=0)
plotMA(res, ylim=c(-6,6), alpha=padjCuotff, cex=0.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
garbage <- dev.off()

# output result
resultFile <- file.path(output, paste(prefix, ".DE.txt", sep="."))
output.file <- file(resultFile, "wb")
write.table(as.data.frame(res), sep="\t", eol = "\n", 
            quote = FALSE, row.names=TRUE, file=output.file)
close(output.file)

# significant output result
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < padjCuotff, pvalue < pvalCutoff)
resultFile <- file.path(output, paste(prefix, ".DE.sig.txt", sep="."))
output.file <- file(resultFile, "wb")
write.table(as.data.frame(res), sep="\t", eol = "\n", 
            quote = FALSE, row.names=TRUE, file=output.file)
close(output.file)
