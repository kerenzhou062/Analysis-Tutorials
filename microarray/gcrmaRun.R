#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",   "Show help information",
    "avereps" ,     "a",   0,  "logical",   "use avereps() to unique genes",
    "gse" ,         "g",   1,  "character", "GSE accession number",
    "gcrma" ,       "r",   0,  "logical",   "Use gcrma() instead of rma()",
    "input" ,       "i",   1,  "character", "input directory of CEL files (.CELL)",
    "cdf" ,         "c",   1,  "character", "cdf file for probe set (eg: hgu133acdf)",
    "output" ,      "o",   1,  "character", "output result file",
    "probe" ,       "p",   1,  "character", "probe-gene matrix from probeToGene.R"
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
ShowHelp(args$probe, '-p|--probe', FALSE)
ShowHelp(args$output, '-o|--output', FALSE)

# default values
if ( is.null(args$avereps) ) { args$avereps = FALSE }
if ( is.null(args$gcrma) ) { args$gcrma = FALSE }
if ( is.null(args$output) ) { args$output = './probe_to_gene.txt' }
if ( is.null(args$gse) ) { args$gse = NULL }
if ( is.null(args$cdf) ) { args$cdf = NULL }

# load arguments
avereps <- args$avereps
gcrma <- args$gcrma
gse <- args$gse
input <- args$input
cdf <- args$cdf
probe <- args$probe
output <- args$output

# load libraries
suppressMessages(library("R.utils"))
suppressMessages(library('affy'))

LoadPacakge("R.utils")
LoadPacakge("affy")

if (gcrma) {
  suppressMessages(library('gcrma'))
  LoadPacakge("gcrma")
}

suppressMessages(library('gcrma'))

if (avereps) {
  suppressMessages(library('limma'))
  LoadPacakge("limma")
}

# ref: https://www.biostars.org/p/53870/

if (! is.null(cdf)) {
  suppressMessages(library(cdf, character.only = TRUE))
  LoadPacakge(cdf)
}

if (! is.null(gse)) {
  #Download the CEL file package if --gse
  getGEOSuppFiles(gse)
  gsePath <- file.path(input, gse)
  setwd(gsePath)
  ## Unpack the CEL files
  untar(paste(gse, "_RAW.tar", sep=""), exdir="data")
  ## setwd to input/gse/data
  dataPath <- file.path(gsePath, "data")
  setwd(dataPath)
  ## decompress CEL.gz
  gzipCels <- list.files("./", pattern = "gz")
  sapply(paste("./", gzipCels, sep="/"), gunzip)
  cels <- list.files("./", pattern = "CEL")
}else{
  setwd(input)
  cels <- list.files("./", pattern = "CEL")
}

rawData <- ReadAffy(verbose=TRUE, filenames=cels, cdfname=cdf)

#perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
if (gcrma) {
  dataRmaNorm <- gcrma(rawData, affinity.info=rawData)
}else{
  dataRmaNorm <- rma(rawData)
}

#Get the important stuff out of the data - the expression estimates for each array
rma = exprs(dataRmaNorm)

#Format values to 5 decimal places
rma=format(rma, digits=5)

#Map probe sets to gene symbols

probeMtx <- read.table(probe, sep="\t", row.names = 1, header = TRUE, quote = "")

commonProbes <- intersect(rownames(rma), rownames(probeMtx))
rma <- rma[commonProbes,]
probeMtx <- probeMtx[commonProbes,]

probeId <- row.names(rma)
geneId <- as.character(probeMtx$ensembl_gene_id)
geneSymbol <- as.character(probeMtx$external_gene_name)

rma=cbind(probeId, geneId, geneSymbol, rma)

if (avereps) {
  rownames(rma) <- as.vector(rma[,"geneSymbol"])
  rma = avereps(rma)
}

# remove .CEL from column names
colnames(rma) <- gsub(pattern = ".CEL|.cel",replacement = "", x  = colnames(rma))

resultFile <- file.path(output)
output.file <- file(resultFile, "wb")
write.table(rma, sep="\t", eol = "\n", 
            quote = FALSE, row.names=FALSE, file=output.file)
close(output.file)

if (! is.null(gse)) {
  ## remove all downloaded files
  gsePath <- file.path(input, gse)
  system(paste0("rm -rf ", gsePath))
}
