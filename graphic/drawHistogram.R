#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "fcol",         "f",   1,  "character",   "name of column that to be filtered",
    "fval",         "v",   1,  "character",   "kept value in --fcol",
    "gzip",         "z",   0,  "logical",     "the --input matrix is gzip",
    "input",        "i",   1,  "character",   "input data matrix (wide)",
    "keycol",       "k",   1,  "character",   "',' separated key column names kept with tidyr.[1st column]",
    "legend",       "l",   1,  "character",   "legend name [condition]",
    "output",       "o",   1,  "character",   "output pdf name",
    "position",     "p",   1,  "character",   "position of geom_bar (dodge, stack) [stack]",
    "xcol",         "x",   1,  "character",   "name of column showed on x-axis [1st column]",
    "ylab",         "y",   1,  "character",   "name of y-lab"
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
ShowHelp(args$ylab, '-y|--ylab')

# default values
if ( is.null(args$gzip) ) { args$gzip = FALSE }
if ( is.null(args$color) ) { args$color = 'Dark2' }
if ( is.null(args$keycol) ) { args$keycol = "1st"}
if ( is.null(args$legend) ) { args$legend = "condition"}
if ( is.null(args$position) ) { args$position = "stack"}
if ( is.null(args$xcol) ) { args$xcol = "1st"}
if ( is.null(args$output) ) { args$output = 'histogram.pdf' }

# load libraries
LoadPacakge('ggplot2')
LoadPacakge('viridis')
LoadPacakge('tidyr')

if (args$gzip) {
  dataMtx <- gzfile(args$input,'rt')
}else{
  dataMtx <- args$input
}

# read expData from input matrix
inputData <- read.table(dataMtx, sep="\t", header=TRUE)
inputData[,1] <- factor(inputData[,1])

## split target name into vector
colNameVector <- colnames(inputData)
if (args$keycol == '1st') {
  keycol = colNameVector[1]
}else{
  keycol = args$keycol
}

if (args$xcol == '1st') {
  xcol = colNameVector[1]
}else{
  xcol = args$xcol
}

# select wrap data
## split target name into vector
setColNameVector <- unlist(strsplit(keycol, ","))
targetCols <- colNameVector[!colNameVector %in% setColNameVector]
tidyrData <- gather(inputData, key = condition, value = value, targetCols)
## rename condition to --legend
names(tidyrData)[names(tidyrData)=="condition"] <- args$legend

## keep data with --fcol
if ((! is.null(args$fcol)) & (! is.null(args$fval))) {
  names(tidyrData)[names(tidyrData)==args$fcol] <- "newFilterCol"
  tidyrData <- tidyrData[tidyrData$newFilterCol == args$fval,]
}

plot <- ggplot(tidyrData, aes_string(fill=args$legend, y="value", x=xcol)) + 
    geom_bar(position=args$position, stat="identity") +
    scale_fill_viridis(discrete = T) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      text=element_text(size=12),
      axis.text.x = element_text(angle=50, hjust=1, vjust=1)) +
    xlab("") + ylab(args$ylab)

pdf(args$output, paper = "a4r", width = 0, height = 0)
print(plot)
gabage <- dev.off()
