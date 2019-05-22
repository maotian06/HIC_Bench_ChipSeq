#!/usr/bin/env Rscript

usage = "\
Rscript chipseq-peakanno.r [OPTIONS] PEAK-FILE
"
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-o","--outputdir"), default = "./", help = "output directory."),
  make_option(c("-d","--promoter"), default = "./", help = "extending promoter upstream and downstream by nt."),
  make_option(c("-g","--genome"), default = "./", help = "annotation genome")
)
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs)!=1) { write(paste0("USAGE: ", usage), stderr()); quit(save = 'no') }

# exit if input file is too short (will cause error in readPeakFile)
fileconn = file(inputs, open = "r")
input_length <- length(readLines(fileconn))
if (input_length < 10) {
  write("too few peaks so skipping", stderr())
  quit(save = 'no')
}

if(grepl("hg19|hg38",opt$genome)){
  TXDB_NAME = paste0("TxDb.Hsapiens.UCSC.",opt$genome,".knownGene");
  annoDb_name = paste0("org.Hs.eg.db");
}else if(grepl("mm9|mm10",opt$genome)){
  TXDB_NAME = paste0("TxDb.Mmusculus.UCSC.",opt$genome,".knownGene");
  annoDb_name = paste0("org.Mm.eg.db");
}else {
  cat("ERROR: Unsupported Genome\n");
  quit(save='no')
}
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(TXDB_NAME,character.only=TRUE))

txdb <- get(TXDB_NAME)

peak <- readPeakFile(inputs)

if(!dir.exists(opt$outputdir)){dir.create(opt$outputdir)}

pdf(file = paste0(opt$outputdir,"/peaks-coverage.pdf"))
covplot(peak, weightCol="V5")
dev.off()

peakAnno <- annotatePeak(peak, tssRegion=c(-as.numeric(opt$promoter), as.numeric(opt$promoter)), TxDb=txdb, annoDb=annoDb_name)

pdf(file = paste0(opt$outputdir,"/anno-piechart.pdf"))
plotAnnoPie(peakAnno)
dev.off()

pdf(file = paste0(opt$outputdir,"/upsetplot.pdf"),width = 9,height=4.5,onefile = F)
upsetplot(peakAnno,vennpie=T)
dev.off()

write.table(peakAnno,file=paste0(opt$outputdir,"/peak_anno.tsv"),quote=F,sep="\t",row.names =F)

sessionInfo()
