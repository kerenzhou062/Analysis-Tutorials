#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",         "h",   0,  "logical",     "Show help information",
    "feature",      "f",   1,  "character",   "feature (coding, exon, intron, full) [coding]",
    "input",        "i",   1,  "character",   "input bin matrix derived from metagene.py",
    "output",       "o",   1,  "character",   "Output directory (./)",
    "prefix",       "p",   1,  "character",   "Output prefix (metagene)",
    "smooth",       "s",   0,  "logical",     "Smooth the curve",
    "ylab",         "y",   1,  "character",   "The name of ylab"
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

# default values
if ( is.null(args$smooth) ) { args$smooth = FALSE }
if ( is.null(args$prefix) ) { args$prefix = 'metagene' }
if ( is.null(args$output) ) { args$output = './' }
if ( is.null(args$feature) ) { args$feature = 'coding' }
if ( is.null(args$ylab) ) { args$ylab = 'Peaks Coverage (%)' }

LoadPacakge("ggplot2")
LoadPacakge("ggpubr")
LoadPacakge("reshape")
LoadPacakge("RColorBrewer")
LoadPacakge("grid")

data = read.table(file=args$input, sep="\t", header=TRUE, row.names=NULL)
valueData <- data[ , -which(names(data) %in% c("feature"))]
melData <- melt(valueData,id=c("bin"))
colourCount <- length(data) - 2
getPalette <- colorRampPalette(brewer.pal(6, "Dark2"))

annoCoord = (max(melData$value)-min(melData$value)) / 20

## library("forcats")
## geom_area(aes(fill= fct_reorder(variable, value, .desc = TRUE)), position = 'identity')
metaGenePlot <- ggplot(melData, aes(x=bin, y=value, group=variable)) +
    geom_area(aes(fill= variable), position = 'identity') +
    geom_line(aes(colour= variable), size=0.5) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      text=element_text(size=12),
      legend.position = 'right',
      legend.title = element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      plot.margin = unit(c(1,1,2,1), "lines")) +
    scale_color_manual(values = getPalette(colourCount)) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    xlab("") + ylab(args$ylab)

if (args$feature != 'coding') {
  featureXlabCoord = as.integer((max(data$bin) + min(data$bin)) / 2)
  featureText <- textGrob(args$feature, gp=gpar(fontface="bold"))
  metaGenePlot <- metaGenePlot +
    annotation_custom(featureText,xmin=featureXlabCoord,xmax=featureXlabCoord,ymin=-annoCoord,ymax=-annoCoord)
}else{
  colourCount <- length(data) - 2
  getPalette <- colorRampPalette(brewer.pal(5, "Dark2"))
  utr5XlabMin <- min(data[data$feature == 'utr5',]$bin)
  utr5XlabMax <- max(data[data$feature == 'utr5',]$bin)
  cdsXlabMin <- min(data[data$feature == 'cds',]$bin)
  cdsXlabMax <- max(data[data$feature == 'cds',]$bin)
  utr3XlabMin <- min(data[data$feature == 'utr3',]$bin)
  utr3XlabMax <- max(data[data$feature == 'utr3',]$bin)
  ## feature text position in x lab
  utr5XlabCoord <- as.integer((utr5XlabMin + utr5XlabMax) / 2)
  cdsXlabCoord <- as.integer((cdsXlabMin + cdsXlabMax) / 2)
  utr3XlabCoord <- as.integer((utr3XlabMin + utr3XlabMax) / 2)
  ## feature text 
  utr5Text <- textGrob("5' UTR", gp=gpar(fontface="bold"))
  cdsText <- textGrob("CDS", gp=gpar(fontface="bold"))
  utr3Text <- textGrob("3' UTR", gp=gpar(fontface="bold"))
  ## metagene
  metaGenePlot <- metaGenePlot +
    geom_vline(aes(xintercept=utr5XlabMax), colour="#4F3833", linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=cdsXlabMax), colour="#4F3833", linetype="dashed", size=0.5) +
    annotation_custom(utr5Text,xmin=utr5XlabCoord,xmax=utr5XlabCoord,ymin=-annoCoord,ymax=-annoCoord) +
    annotation_custom(cdsText,xmin=cdsXlabCoord,xmax=cdsXlabCoord,ymin=-annoCoord,ymax=-annoCoord) +
    annotation_custom(utr3Text,xmin=utr3XlabCoord,xmax=utr3XlabCoord,ymin=-annoCoord,ymax=-annoCoord)
}

metaGenePlot <- metaGenePlot + coord_cartesian(clip = "off")

## print plot to pdf
metagenePlotPdf <- file.path(args$output, paste(args$prefix, ".metagene.pdf", sep=""))

pdf(metagenePlotPdf, width = 8, height = 5)
print(metaGenePlot)
gabage <- dev.off()
