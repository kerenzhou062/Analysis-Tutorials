#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "input",        "i",   1,  "character",   "input matrix",
    "prefix",       "p",   1,  "character",   "output prefix",
    "output",       "o",   1,  "character",   "output directory"
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
  cat(paste("Load package: ", name, ".\n"))
}

# parsing arguments
args <- getopt(command)

ShowHelp(args$input, '-i|--input')

if ( is.null(args$prefix) ) { args$prefix = 'survival' }
if ( is.null(args$output) ) { args$output = './' }

input <- args$input
prefix <- args$prefix
output <- args$output

suppressMessages(library('survival'))
suppressMessages(library("survminer"))

survivalData <- read.csv(input, sep="\t", row.names="Sample", header=TRUE)
sfit <- survfit(Surv(SurvivalDays, Status)~Category, data=survivalData)
maxDay <- max(survivalData$SurvivalDays)
surPlotPdf <- file.path(output, paste(prefix, ".survival.pdf", sep=""))
pdf(surPlotPdf, paper='a4r', height=0, onefile=FALSE)

#names(sfit$strata) <- gsub("Category=", "", names(sfit$strata))

ggsurvplot(sfit, linetype = "strata", 
            conf.int = FALSE, pval = TRUE,
            pval.method=TRUE,
            legend.title="Expression",
            surv.median.line="hv", pval.coord=c(maxDay * 0.8, 0.9),
            pval.method.coord =c(maxDay * 0.8, 1),
            palette = "Dark2",
            title=paste("Kaplan-Meier Curve in BRCA: ", prefix, sep=""))

garbage <- dev.off()

# output p-value
pval <- unlist(strsplit(surv_pvalue(sfit)$pval.txt, " "))[3]
text = paste(prefix, pval, sep="\t")

pvalFile <- file.path(output, paste(prefix, ".pvalue.txt", sep=""))
output.file <- file(pvalFile, "wb")
write.table(text, file = output.file, col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
close(output.file)
