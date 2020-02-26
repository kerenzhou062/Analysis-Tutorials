#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "gzip",         "z",   0,  "logical",     "the --input matrix is gzip",
    "input",        "i",   1,  "character",   "input expression matrix (long)",
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
  cat(paste("Load package: ", name, ".\n", sep="\t"))
}

plotGeneCount <- function(geneName, pdfName) {
  pdf(pdfName, paper="a4r")
  selectData <- expData[expData[args$targetCol] == geneName,]
  pd = position_dodge(width = 0.5)
  plot <- ggplot(selectData, aes_string(x=args$xaxis, y=args$yaxis, fill=args$xaxis)) +
    stat_boxplot(geom="errorbar", position=pd, width=0.2) +
    geom_boxplot() + 
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
if ( is.null(args$prefix) ) { args$prefix = 'FPKM' }
if ( is.null(args$refgroup) ) { args$refgroup = 'none' }
if ( is.null(args$output) ) { args$output = './' }

# load libraries
suppressMessages(library('ggplot2'))
suppressMessages(library('ggpubr'))

LoadPacakge('ggplot2')
LoadPacakge('ggpubr')

# split target name into vector
geneList <- unlist(strsplit(args$target, ","))

if (args$gzip) {
  expMtx <- gzfile(args$input,'rt')
}else{
  expMtx <- args$input
}

# read expData from input matrix
expData <- read.csv(expMtx, sep="\t", row.names=NULL)

# draw boxplots
for (i in 1:length(geneList)) {
  geneName <- geneList[i]
  pdfName <- file.path(args$output, paste(geneName, ".", args$prefix, ".pdf", sep=""))
  plotGeneCount(geneName, pdfName)
}
