#!/usr/bin/env Rscript

usage = "\
Rscript chipseq-venn-diagramr [OPTIONS] PEAK-FILE(s)
"
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-o","--outputdir"),default="./", help="output directory."),
  make_option(c("-f","--filename"),default="venn-diagram.pdf", help="output file name"),
  make_option(c("-g","--genome"),default="./", help="annotation genome")
)
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) < 2 || length(inputs) > 5) {
  write("warning: number of peak files supported is 2 to 5", stderr()); quit(save='no')
}

suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(tools))

files = list()
for (i in 1:length(inputs)){
  # check file size and skip input files that are too short (will cause error in readPeakFile)
  fileconn = file(inputs[i], open = "r")
  input_length <- length(readLines(fileconn))
  if (input_length > 10) {
    files <- append(files, list(inputs[i]))
    names(files)[length(files)] <- basename(dirname(inputs[i]))
  }
}

if(!dir.exists(opt$outputdir)){dir.create(opt$outputdir)}

peak.Sets <- lapply(files, readPeakFile)


tryCatch({
  pdf(file = paste0(opt$outputdir,"/",opt$filename), pointsize=6 )
  vennplot(peak.Sets)
  title(main=opt$filename)
  dev.off()
 },
 error=function(e){
  dev.off()
  write(conditionMessage(e),stdout())
 }
)
