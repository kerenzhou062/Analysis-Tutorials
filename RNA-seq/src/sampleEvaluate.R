#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",   "Show help information",
    "filter",       "f",   2,  "integer",   "Filter out genes less than # counts across all samples",
    "spikein",      "sp",  0,  "logical",   "EstimateSizeFactors with spikeins",
    "spiregex",     "pa",  1,  "character", "Name pattern of spikeins in gene_id (Spikein-ERCC|ERCC)",
    "counts",       "ct",  1,  "character", "Gene counts matrix",
    "samplemtx",    "sm",  1,  "character", "Sample relationships matrix",
    "design1",      "d1",  1,  "character", "Design for construction of DESeqDataSet (colname in colData)",
    "design2",      "d2",  1,  "character", "Design for construction of DESeqDataSet (colname in colData)",
    "poiheatmap",   "po",  1,  "character", "Sample distance plot using Poisson Distance",
    "glmpca",       "gp",  0,  "logical",   "PCA plot using Generalized PCA",
    "normalize",    "nr",  1,  "character", "Normalize method for raw counts (auto|vst|rlog)",
    "prefix",       "pr",  1,  "character", "Prefix for output",
    "output" ,      "o",   1,  "character", "Output directory"
), byrow=TRUE, ncol=5)

## parsing arguments
args <- getopt(command)

if(!is.null(args$help)){
  cat(paste(getopt(command, usage = T),"\n"))
  q()
}

if ( is.null(args$normalize) ) {
  args$normalize = 'auto'
}else{
  norMethodVetor <- c('auto', 'vst', 'rlog')
  if (isFALSE(args$normalize %in% norMethodVetor)) {
    cat("None valid -nr|--normalize!\n")
    cat(paste(getopt(command, usage = T),"\n"))
    q()
  }
}

## functions
VarDefined <- function(x) {
  tryCatch ( 
    if( class(x) == 'logical' ) {
      TRUE
    } else{ 
      TRUE
    },
    error=function(e) {
      FALSE
    }
  )
}

NormalizeData <- function(method, sampleSize) {
  ## normalize data with vst or rlog
  if (method == 'vst') {
    normalizeData <- vst(dds, blind = TRUE)
  }else if (method == 'rlog') {
    normalizeData <- rlog(dds, blind = TRUE)
  }else {
    ## if sampleSize <= 30, rlog; else vst
    if (sampleSize <= 30) {
      normalizeData <- rlog(dds, blind = TRUE)
    }else{
      normalizeData <- vst(dds, blind = TRUE)
    }
  }
  return(normalizeData)
}

## load DESeq2 and perform Differential Analysis
suppressMessages(library('DESeq2'))
## load arguments
geneCountMtx <- args$counts
sampleMtx <- args$samplemtx
design1 <- args$design1
design2 <- args$design2
prefix <- args$prefix
output <- args$output
spiRegex <- args$spiregex
normalize <- args$normalize

# With the count matrix, cts, and the sample information, colData
cts <- as.matrix(read.csv(geneCountMtx, sep="\t", row.names="gene_id"))
colData <- read.csv(sampleMtx, row.names=1, sep="\t")
## check all sample rownames in geneCountMtx colNames
all(rownames(colData) %in% colnames(cts))
cts <- cts[, rownames(colData)]

# if used spikein, use "RUVSeq" to Estimating the factors of unwanted variation using control genes
if( !is.null(args$spikein) ){
  suppressMessages(library('RUVSeq'))
  suppressMessages(library('RColorBrewer'))
  colors <- brewer.pal(3, "Set2")

  ## seperate to genes and spikes
  genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
  spikes <- rownames(cts)[grep(spiRegex, rownames(cts))]
  set <- newSeqExpressionSet(as.matrix(cts), phenoData = colData)
  ## use the betweenLaneNormalization function of EDASeq to normalize the
  ## data using upper-quartile (UQ) normalization
  set <- betweenLaneNormalization(set, which="upper")
  ## plot RLE and PCA before RUVg
  rlePlotPdf <- file.path(output, paste(prefix, ".beforeSpikein.RLE.pdf", sep=""))
  pdf(rlePlotPdf, paper='a4r', height=0)
  plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[spif])
  garbage <- dev.off()
  pcaPlotPdf <- file.path(output, paste(prefix, ".beforeSpikein.PCA.pdf", sep=""))
  pdf(pcaPlotPdf, paper='a4r', height=0)
  plotPCA(set, col=colors[spif], cex=1.2)
  garbage <- dev.off()
  ## RUVg: Estimating the factors of unwanted variation using control genes
  spikeNorSet <- RUVg(set, spikes, k=1)
  ## plot RLE and PCA after RUVg
  rlePlotPdf <- file.path(output, paste(prefix, ".afterSpikein.RLE.pdf", sep=""))
  pdf(rlePlotPdf, paper='a4r', height=0)
  plotRLE(spikeNorSet, outline=FALSE, ylim=c(-4, 4), col=colors[spif])
  garbage <- dev.off()
  pcaPlotPdf <- file.path(output, paste(prefix, ".afterSpikein.PCA.pdf", sep=""))
  pdf(pcaPlotPdf, paper='a4r', height=0)
  plotPCA(spikeNorSet, col=colors[spif], cex=1.2)
  garbage <- dev.off()
  ## pass spikeNorSet to DESeq2
  ## re-construct cts, filter out spike-ins
  cts <- cts[genes,]
  designFormula <- as.formula(paste("~", "W_1", "+", design, sep=" "))
}else{
  ## filter out spike-ins 
  genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
  cts <- cts[genes,]
  designFormula <- as.formula(paste("~", design, sep=" "))
}

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

# use the Poisson Distance to calculate sample distance 
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))

sampleSize <- nrow(colData)
sampleDisPdf <- file.path(output, paste(prefix, ".sd.heatmap.pdf", sep=""))
if(!is.null(args$poiheatmap)) {
  suppressMessages(library("PoiClaClu"))
  poisd <- PoissonDistance(t(counts(dds)))
  sampleDistMatrix <- as.matrix( poisd$dd )
  rownames(sampleDistMatrix) <- paste( dds[[design1]], dds[[design2]], sep=" - " )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors)
  garbage <- dev.off()
}else{
  normalizeData <- NormalizeData(normalize, sampleSize)
  sampleDists <- dist(t(assay(normalizeData)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( normalizeData[[design1]], normalizeData[[design2]], sep = " - " )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  garbage <- dev.off()
}


# plot PCA

if(!is.null(args$glmpca)) {
  suppressMessages(library("glmpca"))
  gpca <- glmpca(counts(dds), L=2)
  gpca.dat <- gpca$factors
  gpca.dat[[design1]] <- dds[[design1]]
  gpca.dat[[design2]] <- dds[[design2]]
  pcaPdf <- file.path(output, paste(prefix, ".GLM.PCA.pdf", sep=""))
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell)) +
    geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
  garbage <- dev.off()
}else{
  if (! VarDefined(normalizeData)) {
    normalizeData <- NormalizeData(normalize, sampleSize)
  }
  
  pcaData <- plotPCA(normalizeData, intgroup = c( design1, design2), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaPdf <- file.path(output, paste(prefix, ".PCA.pdf", sep=""))
  ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle(paste("PCA plot with normalized data", " (", normalize, ")", sep=""))
  garbage <- dev.off()
}

# plot MDS
mdsPdf <- file.path(output, paste(prefix, ".MDS.pdf", sep=""))
if(!is.null(args$poiheatmap)) {
  mdsPois <- as.data.frame(colData(dds)) %>%
     cbind(cmdscale(sampleDistMatrix))
  ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
    geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
  garbage <- dev.off()
}else{
  mds <- as.data.frame(colData(normalizeData))  %>%
           cbind(cmdscale(sampleDistMatrix))
  ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
    geom_point(size = 3) + coord_fixed() +
    ggtitle(paste("MDS plot with normalized data", " (", normalize, ")", sep=""))
  garbage <- dev.off()
}
