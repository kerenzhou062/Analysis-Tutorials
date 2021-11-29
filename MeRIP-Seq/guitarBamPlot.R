#!/usr/bin/env Rscript

suppressMessages(library('getopt'))

command =  matrix(c(
    "help",     "h", 0, "logical", "print helps",
    "sqlite",   "s", 1, "character", "txdb_file(e.g. gencode.v24lift37.annotation.sqlite)",
    "name",     "n", 1, "character", "name matched input files (',' seperated)",
    "input",    "i", 1, "character", "input files (',' seperated)",
    "neighbor", "N", 0, "logical", "awake NeighborDNA in GuitarPlot()",
    "height",   "H", 1, "integer", "PDF height",
    "width",    "W", 1, "integer", "PDF width",
    "output",   "o", 1, "character", "output file prefix (include path)"
), byrow=TRUE, ncol=5)

args <- getopt(command)

if(!is.null(args$help)){
  cat(paste(getopt(command, usage = T),"\n"))
  q()
}

suppressMessages(library('Guitar'))

featureList = list()
inputNameList <- strsplit(args$name, ",")[[1]]
inputFileList <- strsplit(args$input, ",")[[1]]
length <- length(inputNameList)
for (i in c(1:length)){
  featureList[inputNameList[i]] = readGAlignments(inputFileList[i])
}

txdbFile <- args$sqlite
TxdbOb <- loadDb(txdbFile)
gcTxdb <- makeGuitarCoordsFromTxDb(TxdbOb, noBins=100, minimalComponentLength=50)

if(!is.null(args$neighbor)){
  GuitarPlot(gfeatures = featureList,
    GuitarCoordsFromTxDb = gcTxdb,
    includeNeighborDNA =TRUE,
    rescaleComponent=FALSE,
    PDFHeight=args$height,
    PDFWidth=args$width,
    saveToPDFprefix = args$output)
}else{
  GuitarPlot(gfeatures = featureList,
    GuitarCoordsFromTxDb = gcTxdb,
    includeNeighborDNA =FALSE,
    rescaleComponent=FALSE,
    PDFHeight=args$height,
    PDFWidth=args$width,
    saveToPDFprefix = args$output)
}

