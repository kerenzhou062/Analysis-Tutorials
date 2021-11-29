#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",   "Show help information",
    "input" ,       "i",   1,  "character", "Input directory"
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
  cat(paste("Load package: ", name, ".\n"))
}

# parsing arguments
args <- getopt(command)

ShowHelp(args$help, 'none', TRUE)
ShowHelp(args$input, '-i|--input', FALSE)

# default values
if ( is.null(args$input) ) { args$input = './' }

# load arguments
input <- args$input

setwd(input)

# load DESeq2
suppressMessages(library('readDepth'))

# create a readDepth object, then fill it by
# reading in the params, setting up the environment,
# creating the model, and choosing optimal bin size
rdo = new("rdObject")

# calculate depth of coverage in each bin
rdo = readDepth(rdo)

# correct the reads for mapability. This example uses a conservative
# threshold of 0.75. In other words, if a bin is less than 75% mapable,
# it's depth is set to NA. This prevents overcorrection.
rdo.mapCor = rd.mapCorrect(rdo, minMapability=0.75)

# do LOESS-based GC correction.
rdo.mapCor.gcCor = rd.gcCorrect(rdo.mapCor)

# segment the data using CBS. If you notice artifacts in the output, such
# as regions of gain that span centromeres, you might try adding the
# "rmGaps=FALSE" parameter. If you're using data with very high coverage
# (say, greater than 10x), consider adding "minWidth=3" (maybe even 4 or 5)
# to reduce the number of false positives (at the expense of sensitivity)
segs = rd.cnSegments(rdo.mapCor.gcCor)

# write all the segments out to the output directory
writeSegs(segs)

# If you want just the alterations, you can write those out too
writeAlts(segs,rdo)

#write the window size and CN gain/loss thresholds to the outdir
writeThresholds(rdo)

# (optional) save an image of your R session so that you can come
# back and rerun parts of the analysis without redoing the
# whole thing.
