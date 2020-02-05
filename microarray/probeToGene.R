#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",   "Show help information",
    "dataset" ,     "d",   1,  "character", "dataset name in biomart [hsapiens_gene_ensembl]",
    "gse" ,         "g",   1,  "character", "GSE accession number",
    "grep" ,        "r",   1,  "character", "strings for filtering attr(gset, 'names') (eg:GPL)",
    "output" ,      "o",   1,  "character", "output file",
    "probe" ,       "p",   1,  "character", "probe set id (eg:affy_hg_u133a)"
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
ShowHelp(args$gse, '-g|--gse', FALSE)
ShowHelp(args$grep, '-r|--grep', FALSE)
ShowHelp(args$output, '-o|--output', FALSE)

# default values
if ( is.null(args$output) ) { args$output = './probe_to_gene.txt' }

# load arguments
gse <- args$gse
dataset <- args$dataset
grep <- args$grep
probe <- args$probe
output <- args$output

# load libraries
suppressMessages(library("readr", lib.loc = .libPaths()[2]))
suppressMessages(require('GEOquery'))
suppressMessages(require('Biobase'))
suppressMessages(require('biomaRt'))
#https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html

# get data matrix from gse
gset <- getGEO(gse, GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) {
  idx <- grep(grep, attr(gset, "names"))
}else{
  idx <- 1
}

gset <- gset[[idx]]

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset(dataset, mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    probe,
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = probe,
  values = rownames(exprs(gset)), uniqueRows=TRUE)

names(annotLookup)[1] <- "probe_id"

# keep orginal order
annotLookup <- annotLookup[match(rownames(gset), annotLookup$probe_id),]
annotLookup <- annotLookup[!is.na(annotLookup$probe_id),]
# remove duplicate probe_id, just keep the first one
annotLookup <- annotLookup[!duplicated(annotLookup$probe_id),]

resultFile <- file.path(output)
output.file <- file(resultFile, "wb")
write.table(annotLookup, sep="\t", eol = "\n", 
            quote = FALSE, row.names=FALSE, file=output.file)
close(output.file)
