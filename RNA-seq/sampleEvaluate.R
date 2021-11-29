#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",   "Show help information",
    "batchMethod",  "b",   1,  "character",   "Remove hidden batch effect (none|spikeins|upper)",
    "counts",       "g",   1,  "character", "Gene counts matrix (gene_id|tx_id in 1st column)",
    "design1",      "d",   1,  "character", "Design for construction of DESeqDataSet (colname in colData)",
    "design2",      "D",   1,  "character", "Design for construction of DESeqDataSet (colname in colData)",
    "filter",       "f",   2,  "integer",   "Filter out genes less than # counts across all samples",
    "glmPca",       "l",   0,  "logical",   "PCA plot using Generalized PCA",
    "keepSpike",    "k",   1,  "logical",   "Keep spikeins reads when passing to DESeq()",
    "normalize",    "n",   1,  "character", "Normalize method for raw counts (auto|vst|rlog)",
    "output" ,      "o",   1,  "character", "Output directory",
    "poiHeatmap",   "H",   0,  "logical",   "Sample distance plot using Poisson Distance",
    "prefix",       "e",   1,  "character", "Prefix for output",
    "sampleMtx",    "m",   1,  "character", "Sample relationships matrix",
    "spiRegex",     "r",   1,  "character", "Name pattern of spikeins in gene_id (ERCC)"
), byrow=TRUE, ncol=5)

# functions
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

NormalizeData <- function(dds, method, sampleSize) {
  ## normalize data with vst or rlog
  if (method == 'vst') {
    normalizeData <- vst(dds, blind = TRUE)
  }else if (method == 'rlog') {
    normalizeData <- rlog(dds, blind = TRUE)
  }else if (method == 'none') {
    normalizeData <- dds
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

RatioGgplot <- function(data1, data2) {
  span1 <- max(data1) - min(data1)
  span2 <- max(data2) - min(data2)
  maxVal = max(c(span1, span2))
  minVal = min(c(span1, span2))
  ratio <- round(maxVal/minVal, digits = 1)
  return(ratio)
}

AxisLim <- function(data1, data2) {
  maxVal = max(c(max(data1), max(data2)))
  minVal = min(c(min(data1), min(data2)))
  return (c(minVal, maxVal))
}

# parsing arguments
args <- getopt(command)

ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$counts, '-g|--counts', FALSE)
ShowHelp(args$design1, '-d|--design1', FALSE)
ShowHelp(args$design2, '-D|--design2', FALSE)

if ( is.null(args$normalize) ) {
  args$normalize = 'auto'
}else{
  norMethodVetor <- c('auto', 'vst', 'rlog')
  bool <- isFALSE(args$normalize %in% norMethodVetor)
  ShowHelp(bool, '-n|--normalize', FALSE, TRUE)
}

if ( is.null(args$batchMethod) ) {
  args$batchMethod = 'none'
}else{
  bmVetor <- c('none', 'spikeins', 'upper')
  bool <- isFALSE(args$batchMethod %in% bmVetor)
  ShowHelp(bool, '-b|--batchMethod', FALSE, TRUE)
}

if ( is.null(args$test) ) {
  args$test = 'Wald'
}else{
  testVetor <- c('Wald', 'LRT')
  bool <- isFALSE(args$test %in% testVetor)
  ShowHelp(bool, '-T|--test', FALSE, TRUE)
}

# default values
if ( is.null(args$filter) ) { args$filter = 1 }
if ( is.null(args$spiRegex) ) { args$spiRegex = 'ERCC-' }
if ( is.null(args$glmPca) ) { args$glmPca = FALSE }
if ( is.null(args$poiHeatmap) ) { args$poiHeatmap = FALSE }
if ( is.null(args$prefix) ) { args$prefix = 'result' }
if ( is.null(args$output) ) { args$output = './' }
if ( is.null(args$keepSpike) ) { args$keepSpike = FALSE }

# load arguments
geneCountMtx <- args$counts
sampleMtx <- args$sampleMtx
filter <- args$filter
test <- args$test
design1 <- args$design1
design2 <- args$design2
prefix <- args$prefix
output <- args$output
spiRegex <- args$spiRegex
normalize <- args$normalize
keepSpike <- args$keepSpike
glmPca <- args$glmPca
poiHeatmap <- args$poiHeatmap

# load DESeq2
LoadPacakge('DESeq2')
LoadPacakge('ggplot2')
LoadPacakge('dplyr')
# With the count matrix, cts, and the sample information, colData
cts <- as.matrix(read.csv(geneCountMtx, sep="\t", row.names=1))
cts <-round(cts, 0)
colData <- read.csv(sampleMtx, row.names=1, sep="\t")
## check all sample rownames in geneCountMtx colNames
all(rownames(colData) %in% colnames(cts))
cts <- cts[, rownames(colData)]

# filter out gSpikein_phiX174, ENCODE
genes <- rownames(cts)[grep('_phiX174', rownames(cts), invert=TRUE)]
cts <- cts[genes,]

if( args$batchMethod == "spikeins" | args$batchMethod == "upper" ){
  ## Removing hidden batch effects using spike-in controls by RUVg
  cat(paste('Removing hidden batch effects with ', args$batchMethod,' (RUVg)!\n', sep=""))
  LoadPacakge('RUVSeq')
  LoadPacakge('dplyr')
  LoadPacakge("RColorBrewer")
  set <- newSeqExpressionSet(as.matrix(cts), phenoData = colData)
  if (args$batchMethod == "upper") {
    ## use the betweenLaneNormalization function of EDASeq to normalize the
    ## data using upper-quartile (UQ) normalization
    norSet <- betweenLaneNormalization(set, which="upper")
  }else{
    ## RUVg: Estimating the factors of unwanted variation using control genes
    spikes <- rownames(cts)[grep(spiRegex, rownames(cts))]
    norSet <- RUVg(set, spikes, k=1)
  }
  # plot the PCA plot after normalization
  #colors <- colorRampPalette(brewer.pal(8, "Set2"))(nrow(colData))
  #pcaPlotPdf <- file.path(output, paste(prefix, ".ruvseq_pca.pdf", sep=""))
  #pdf(pcaPlotPdf, paper='a4r', height=20, width=20)
  #plotPCA(norSet, cex=1.2)
  cts <- norSet@assayData[["normalizedCounts"]]
  cts <-round(cts, 0)
}

sampleSize <- nrow(colData)
# pre-filtering, counts > args$filter in at least half of the samples
filtered <- apply(cts, 1, function(x) length(x[x>filter])>=round(sampleSize/2))
cts <- cts[filtered,]

if(! keepSpike){
  genes <- rownames(cts)[grep(spiRegex, rownames(cts), invert=TRUE)]
  cts <- cts[genes,]
}

designFormula <- as.formula(paste("~ ", design1, sep=""))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = designFormula)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# runing DESeq
dds <- DESeq(dds, test=test)

# use the Poisson Distance to calculate sample distance 
LoadPacakge('pheatmap')
LoadPacakge('RColorBrewer')

sampleDisPdf <- file.path(output, paste(prefix, ".sd.heatmap.pdf", sep=""))
if(poiHeatmap) {
  suppressMessages(library("PoiClaClu"))
  poisd <- PoissonDistance(t(counts(dds)))
  sampleDistMatrix <- as.matrix( poisd$dd )
  rownames(sampleDistMatrix) <- paste( dds[[design1]], dds[[design2]], sep=" - " )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(sampleDisPdf, paper='a4r', height=0)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = poisd$dd,
           clustering_distance_cols = poisd$dd,
           col = colors)
  garbage <- dev.off()
}else{
  normalizeData <- NormalizeData(dds, normalize, sampleSize)
  sampleDists <- dist(t(assay(normalizeData)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( normalizeData[[design1]], normalizeData[[design2]], sep = " - " )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf(sampleDisPdf, paper='a4r', height=0)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  garbage <- dev.off()
}

# plot PCA
## scale_shape_manual(values=1:nlevels(factor(gpca.dat[[design2]]))) +
pcaPdf <- file.path(output, paste(prefix, ".PCA.pdf", sep=""))
pdf(pcaPdf, paper='a4r', height=0)
if(glmPca) {
  suppressMessages(library("glmpca"))
  gpca <- glmpca(counts(dds), L=2)
  gpca.dat <- gpca$factors
  gpca.dat[[design1]] <- dds[[design1]]
  gpca.dat[[design2]] <- dds[[design2]]
  #ratio <- RatioGgplot(gpca.dat$dim1, gpca.dat$dim2)
  axisLimVector <- AxisLim(gpca.dat$dim1, gpca.dat$dim2)
  title <- "glmpca - Generalized PCA"
  print(ggplot(gpca.dat, aes_string(x = 'dim1', y = 'dim2', color = design1, shape = design2)) +
    geom_point(size =3) + 
    coord_fixed(xlim = axisLimVector, ylim = axisLimVector) + 
    ggtitle(title))
  garbage <- dev.off()
}else{
  if (! VarDefined(normalizeData)) {
    normalizeData <- NormalizeData(dds, normalize, sampleSize)
  }
  pcaData <- plotPCA(normalizeData, intgroup = c( design1, design2), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  #ratio <- RatioGgplot(pcaData$PC1, pcaData$PC2)
  axisLimVector <- AxisLim(pcaData$PC1, pcaData$PC2)
  title <- paste("PCA plot with normalized data", " (", normalize, ")", sep="")
  print(ggplot(pcaData, aes_string(x = 'PC1', y = 'PC2', color = design1, shape = design2)) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed(xlim = axisLimVector, ylim = axisLimVector) +
    ggtitle(title))
  garbage <- dev.off()
}

# plot MDS
if(!is.null(args$poiHeatmap)) {
  mds <- as.data.frame(colData(dds)) %>% cbind(cmdscale(sampleDistMatrix))
  title <- "MDS with PoissonDistances"
}else{
  mds <- as.data.frame(colData(normalizeData))  %>% cbind(cmdscale(sampleDistMatrix))
  title <- paste("MDS plot with normalized data", " (", normalize, ")", sep="")
}

#ratio <- RatioGgplot(mds$`1`, mds$`2`)
axisLimVector <- AxisLim(mds$`1`, mds$`2`)

mdsPdf <- file.path(output, paste(prefix, ".MDS.pdf", sep=""))
pdf(mdsPdf, paper='a4r', height=0)
ggplot(mds, aes_string(x = '`1`', y = '`2`', color = design1, shape = design2)) +
  geom_point(size = 3) + 
  coord_fixed(xlim = axisLimVector, ylim = axisLimVector) + 
  ggtitle(title)
garbage <- dev.off()
