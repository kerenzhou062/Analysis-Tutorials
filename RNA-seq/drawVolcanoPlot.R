#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",       "Show help information",
    "fc",           "f",   1,  "numeric",       "log2FoldChange cutoff [1.0]",
    "fcCol",        "c",   1,  "character",     "log2FoldChange column name [log2FoldChange]",
    "input",        "i",   1,  "character",     "Input matrix (with header, <tab> separated)",
    "idCol",        "d",   1,  "character",     "Column name used for row.names",
    "idName",       "n",   1,  "character",     "IDs used to highlight (',' separated)",
    "geneName",     "g",   1,  "character",     "Gene names used to highlight (',' separated) (if set, corresponds to --idName)",
    "pointSize",    "e",   1,  "numeric",       "pointSize value for EnhancedVolcano [2.0]",
    "sig",          "s",   1,  "numeric",       "Significance cutoff",
    "sigCol",       "k",   1,  "character",     "Significance column name (pvalue, padj or others) [padj]",
    "sigType",      "p",   1,  "character",     "Type of the significance (pvalue|padj) [padj]",
    "title",        "t",   1,  "character",     "Title of output",
    "xlab",         "x",   1,  "character",     "Name of xlab",
    "xlim",         "a",   1,  "numeric",       "limit of xlab",
    "ylab",         "y",   1,  "character",     "Name of ylab",
    "ylim",         "b",   1,  "numeric",       "limit of ylab",
    "output" ,      "o",   1,  "character",     "Output pdf file (e.g. ./results.pdf)"
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

# parsing arguments
args <- getopt(command)

ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$input, '-i|--input')
ShowHelp(args$idCol, '-d|--idCol')
ShowHelp(args$output, '-o|--output')

# default values
if ( is.null(args$fc) ) { args$fc = 1.0 }
if ( is.null(args$sig) ) { args$sig = 0.1 }
if ( is.null(args$fcCol) ) { args$fcCol = 'log2FoldChange' }
if ( is.null(args$sigCol) ) { args$sigCol = 'padj' }
if ( is.null(args$title) ) { args$title = 'Volcano plot' }
if ( is.null(args$pointSize) ) { args$pointSize = 2.0 }
if ( is.null(args$sigType) ) { args$sigType = 'padj' }
if ( is.null(args$xlab) ) { args$xlab = 'Log2 Fold Change' }
if ( is.null(args$ylab) ) { args$ylab = '-log(padj)' }

LoadPacakge('EnhancedVolcano')

## load arguments
input <- args$input
idCol <- args$idCol
fc <- args$fc
sig <- args$sig
fcCol <- args$fcCol
pointSize <- args$pointSize
sigCol <- args$sigCol
sigType <- args$sigType
title <- args$title
xlab <- args$xlab
ylab <- args$ylab
output <- args$output

if ( sigType != 'padj' ) {
  sigType = 'pvalue'
}

if ( sigType == 'padj' ) {
  legendLabels <- c('Not Sig.', 'Log2(FC)', 'padj', 'padj & Log2(FC)')
}else{
  legendLabels <- c('Not Sig.', 'Log2(FC)', 'p-value', 'p-value & Log2(FC)')
}

# import data matrix
resData <- read.csv(input, header=TRUE, sep="\t", quote = "", row.names=c(idCol))

# get xlim
maxVal = max( resData[c(fcCol)] )
minVal = min( resData[c(fcCol)] )

if ( is.null(args$xlim) ) {
  xlimMax = ceiling (max( c(abs(maxVal), abs(minVal)) ))
}else{
  xlimMax = args$xlim
}

if ( is.null(args$ylim) ) {
  ylimMax = -log10( min( resData[c(sigCol)] ) )
}else{
  ylimMax = args$ylim
}

# try to replace geneid by geneName
if ( ! is.null(args$idName)) {
  selectLabIdVec <- unlist(strsplit(args$idName, ","))
  if ( ! is.null(args$geneName) ) {
    selectLabNameVec <- unlist(strsplit(args$geneName, ","))
    if ( length(selectLabIdVec) != length(selectLabNameVec) ) {
      cat(paste("--idName not equal to --geneName\n"))
      q()
    }else{
      for(i in 1:length(selectLabIdVec)){
        id = selectLabIdVec[i]
        gene = selectLabNameVec[i]
        rownames(resData)[rownames(resData) == id] <- gene
      }
    }
    selectLab = selectLabNameVec
  }else{
    selectLab = selectLabIdVec
  }
}else{
  selectLab = c('')
}

# plot the volcano plot
p <- EnhancedVolcano(resData,
    lab = rownames(resData),
    x = fcCol,
    y = sigCol,
    xlim = c(-xlimMax, xlimMax),
    ylim = c(0, ylimMax),
    selectLab = selectLab,
    xlab = xlab,
    ylab = ylab,
    title = title,
    pCutoff = sig,
    FCcutoff = fc,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.8,
    pointSize = pointSize,
    labSize = 3.0,
    colAlpha = 1,
    legendLabels=legendLabels,
    legendPosition = 'top',
    legendLabSize = 10,
    legendIconSize = 4.0)

pdf(output, paper='a4r', height=0)

print(p)

garbage <- dev.off()
