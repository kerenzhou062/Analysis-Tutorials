#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "gzip",         "z",   0,  "logical",     "the --input matrix is gzip",
    "input",        "i",   1,  "character",   "input data matrix (long)",
    "method",       "m",   1,  "character",   "statistic method for comparing means [wilcox.test|t.test]",
    "output",       "o",   1,  "character",   "output directory",
    "prefix",       "p",   1,  "character",   "{target}.{prefix}.pdf",
    "refgroup",     "r",   1,  "character",   "column name of reference group",
    "targetCol",    "c",   1,  "character",   "column name of target names (eg.: gene_name)",
    "target",       "t",   1,  "character",   "',' separated target names (eg. METTL3,METTL14)",
    "xaxis",        "x",   1,  "character",   "column name used as values on x-axis (eg.:cancerType)",
    "yaxis",        "y",   1,  "character",   "column name used as values on y-axis (eg.:expression)"
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

plotGeneCount <- function(target, pdfName) {
  pdf(pdfName, paper="a4r")
  selectData <- inputData[inputData[args$targetCol] == target,]
  pd = position_dodge(width = 0.5)
  plot <- ggplot(selectData, aes_string(x=args$xaxis, y=args$yaxis, fill=args$xaxis)) +
    stat_boxplot(geom="errorbar", position=pd, width=0.2) +
    geom_boxplot(outlier.shape=NA) + 
    geom_jitter(position=position_jitter(width=0.1, height=0.1)) + 
    theme(axis.text.x = element_text(angle=80, hjust=1, vjust=0.98)) +
    labs(y = paste("Expression Level: ", args$prefix), x = "")
  ## add significance level
  if (args$refgroup != 'none') {
    plot <- plot + 
      stat_compare_means(method = args$method, ref.group = args$refgroup, 
          label = "p.signif" )
  }
  print(plot)
  garbage <- dev.off()
}

# parsing arguments
args <- getopt(command)

ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$input, '-i|--input')
ShowHelp(args$xaxis, '-x|--xaxis')
ShowHelp(args$yaxis, '-y|--yaxis')
ShowHelp(args$targetCol, '-c|--targetCol')
ShowHelp(args$target, '-t|--target')
ShowHelp(args$output, '-o|--output')

# default values
if ( is.null(args$gzip) ) { args$gzip = FALSE }
if ( is.null(args$method) ) { args$method = "wilcox.test" }
if ( is.null(args$prefix) ) { args$prefix = 'boxplot' }
if ( is.null(args$refgroup) ) { args$refgroup = 'none' }
if ( is.null(args$output) ) { args$output = './' }

# load libraries
LoadPacakge('ggplot2')
LoadPacakge('ggpubr')

# split target name into vector
targetList <- unlist(strsplit(args$target, ","))

if (args$gzip) {
  dataMtx <- gzfile(args$input,'rt')
}else{
  dataMtx <- args$input
}

# read inputData from input matrix
inputData <- read.csv(dataMtx, sep="\t", row.names=NULL)

# draw boxplots
for(i in 1:length(targetList)) {
  target <- targetList[i]
  pdfName <- file.path(args$output, paste(target, ".", args$prefix, ".pdf", sep=""))
  plotGeneCount(target, pdfName)
}
