#!/usr/bin/env Rscript

usage = "\
Rscript chipseq-summ-peakanno.r [OPTIONS] PEAK-FILE(s)
"
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-o","--outputdir"),default="./", help="output directory."),
  make_option(c("-d","--promoter"),default="./", help="extending promoter upstream and downstream by nt."),
  make_option(c("-g","--genome"),default="./", help="annotation genome")
)
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs)<1) { write(paste("USAGE: ",usage,sep=''),stderr()); quit(save='no') }


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
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(TXDB_NAME,character.only=TRUE))
txdb <- get(TXDB_NAME)


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

if (length(files) == 0) {
  write("all files have too few peaks so skipping", stderr())
  quit(save = 'no')
}

if (!dir.exists(opt$outputdir)) { dir.create(opt$outputdir) }


#promoter <- getPromoters(TxDb=txdb, upstream=as.numeric(opt$promoter), downstream=as.numeric(opt$promoter))
#tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

#pdf(file = paste0(opt$outputdir,"/average-profile.pdf"))
#plotAvgProf(tagMatrixList, xlim=c(-as.numeric(opt$promoter), as.numeric(opt$promoter)), conf=0.95,resample=500, facet="row")
#dev.off()

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-as.numeric(opt$promoter), as.numeric(opt$promoter)), verbose=FALSE)

pdf(file = paste0(opt$outputdir,"/annotation-barplot.pdf"))
plotAnnoBar(peakAnnoList)
dev.off()

pdf(file = paste0(opt$outputdir,"/distribution-to-TSS.pdf"))
plotDistToTSS(peakAnnoList,title="Distribution of peaks loci\nrelative to TSS")
dev.off()

tryCatch({
pdf(file = paste0(opt$outputdir,"/KEGG-pathway.pdf"))
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
plot(compKEGG, showCategory = 15, font.size = 6, title = "KEGG Pathway Enrichment Analysis")
dev.off()
 },
 error=function(e){
  dev.off()
  write(conditionMessage(e),stdout())
 }
)

if(file.size(paste0(opt$outputdir,"/KEGG-pathway.pdf"))=="3611"){system(paste0("rm -rf ",opt$outputdir,"/KEGG-pathway.pdf"))}

sessionInfo()
