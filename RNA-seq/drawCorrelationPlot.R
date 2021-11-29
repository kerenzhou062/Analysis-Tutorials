#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "drop",         "d",   1,  "character",   "drop columns of dataframe (start,end)",
    "input",        "i",   1,  "character",   "input expression matrix used to construct dataframe (row:gene, column:sample, row.names = 1st column)",
    "method",       "m",   1,  "character",   "method for evaluate the correlations [pearson|spearman|kendall]",
    "output",       "o",   1,  "character",   "output directory",
    "prefix",       "p",   1,  "character",   "prefix of output pdf or png",
    "size",         "s",   1,  "int",         "point size in the plot",
    "type",         "t",   1,  "character",   "output image type [pdf|png]",
    "xgene",        "x",   1,  "character",   "expression data of gene_id used for x axis",
    "ygene",        "y",   1,  "character",   "expression data of gene_id used for y axis"
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

PlotCorr <- function(data, args) {
  if (args$type == "png") {
    outputFile <- file.path(args$output, paste(args$prefix, ".png", sep=""))
    png(file = outputFile, width = 4, height = 4, units = 'in', res = 300, bg = "transparent")
  }else{
    outputFile <- file.path(args$output, paste(args$prefix, ".pdf", sep=""))
    pdf(outputFile, paper="a4r")
  }
  plot <- ggscatter(data, x = args$xgene, y = args$ygene,
          color = "black", shape = 20, size = args$size,
          add = "reg.line", # Add regressin line
          conf.int = TRUE,
          add.params = list(color = "#1E90FF", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE,
          cor.method = args$method,
          cor.coeff.args = list(method = args$method, label.sep = "\n"),
          xlab = args$xgene, ylab = args$ygene)
  print(plot)
  garbage <- dev.off()
}

# parsing arguments
args <- getopt(command)

ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$input, '-i|--input')
ShowHelp(args$xgene, '-x|--xgene')
ShowHelp(args$ygene, '-y|--ygene')

if ( is.null(args$method) ) {
  args$method = 'pearson'
}else{
  methodVetor <- c('pearson', 'spearman', 'kendall')
  bool <- isFALSE(args$method %in% methodVetor)
  ShowHelp(bool, '-m|--method', FALSE, TRUE)
}

if ( is.null(args$type) ) {
  args$type = 'png'
}else{
  typeVetor <- c('png', 'pdf')
  bool <- isFALSE(args$type %in% typeVetor)
  ShowHelp(bool, '-t|--type', FALSE, TRUE)
}

# default values
if ( is.null(args$size) ) { args$size = 3 }
if ( is.null(args$prefix) ) { args$prefix = 'coexpression' }
if ( is.null(args$output) ) { args$output = './' }

# load libraries
LoadPacakge('ggplot2')
LoadPacakge("ggpubr")

myData <- read.csv(args$input, sep="\t", header = TRUE, row.names = 1)

if (! is.null(args$drop) ){
  dropList = unlist(strsplit(args$drop, ","))
  dropStart = strtoi(dropList[1])
  dropEnd = strtoi(dropList[2])
  myData <- myData[, -c(dropStart:dropEnd)]
}

myData <- as.data.frame(t(myData))

if( (args$xgene %in% names(myData)) & (args$ygene %in% names(myData)) ) {
  PlotCorr(myData, args)
}
