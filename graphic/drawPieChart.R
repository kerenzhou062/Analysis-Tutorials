#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "input",        "i",   1,  "character",   "STAR mapping stats matrix from getStarMapStats.py",
    "label",        "l",   1,  "character",   "name of label column",
    "output",       "o",   1,  "character",   "output pdf",
    "title",        "t",   1,  "character",   "title for pie chart",
    "value",        "v",   1,  "character",   "name of value column (count)"
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
ShowHelp(args$label, '-l|--label')
ShowHelp(args$value, '-v|--value')
ShowHelp(args$output, '-o|--output')

# default values
if ( is.null(args$title) ) { args$title = "Pie Chart" }
if ( is.null(args$output) ) { args$output = './piechart.pdf' }

# load libraries
LoadPacakge('plotly')

## read countData from input matrix
countData <- read.table(args$input, row.names=NULL, header = TRUE, sep="\t")
#
## draw piechart

p <- plot_ly(countData,
  labels = as.formula(paste('~', args$label, sep="")), 
  values = as.formula(paste('~', args$value, sep="")),
  type = 'pie') %>%
  layout(title = args$title,
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

options(warn=-1)

# webshot::install_phantomjs()
export(p, file = args$output)
