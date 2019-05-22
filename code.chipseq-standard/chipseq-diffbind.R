#!/usr/bin/env Rscript


# output width
options(width = 120)

# process command-line arguments (only arguments after --args)
args = commandArgs(trailingOnly = TRUE)
outDir = args[1]
sampleSheetCsv = args[2]
normMethod = args[3]          # e.g. DESEQ2
genome = args[4]
genomeDir = args[5]
blockFactor = args[6]

# load libraries
library(DiffBind, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(ChIPpeakAnno, quietly = TRUE)

# setup normalization 
if (normMethod=='DESEQ2') {
  print(normMethod)
  normMethodObj = DBA_DESEQ2
  normMethodBlock = DBA_DESEQ2_BLOCK
} else if (normMethod=='EDGER') {
  print(normMethod)
  normMethodObj = DBA_EDGER
  normMethodBlock = DBA_EDGER_BLOCK
} else {
  write("Error: unknown normalization method!", stderr())
  quit(save='no')
}


# extract differentially bound peaks and annotate them
generateDiffBindReport = function(dba, contrast, th = 0.05, method = DBA_DESEQ2, reps = FALSE, tss, mart.df, out.dir = ".") {

  library(DiffBind, quietly = TRUE)
  library(ChIPpeakAnno, quietly = TRUE)

  # extract contrast group names
  contrasts = as.matrix(dba.show(dba, bContrasts = TRUE))
  group1 = contrasts[as.character(contrast), "Group1"]
  group2 = contrasts[as.character(contrast), "Group2"]

  # generate report (GRanges object)
  message("[generateDiffBindReport] generate diffbind report")
  message("[generateDiffBindReport] contrast num: ", contrast)
  message("[generateDiffBindReport] group1: ", group1)
  message("[generateDiffBindReport] group2: ", group2)
  message("[generateDiffBindReport] threshold: ", th)
  db.gr = dba.report(dba, contrast = contrast, method = method, th = th, bCounts = reps, DataType = DBA_DATA_GRANGES)

  message("[generateDiffBindReport] num sig peaks: ", length(db.gr))

  # annotate if any significant results were found
  if (length(db.gr)) {

    message("[generateDiffBindReport] annotate peaks")
    db.ann.gr = annotatePeakInBatch(db.gr, AnnotationData = tss, PeakLocForDistance = "middle",
                                    FeatureLocForDistance = "TSS", output = "shortestDistance", select = "first")

    # add gene symbols
    message("[generateDiffBindReport] add gene symbols")
    db.ann.df = merge(as.data.frame(db.ann.gr), mart.df, by.x = c("feature"), by.y = c("ensembl_gene_id"), all.x = TRUE)

    # keep just the relevant columns
    db.ann.df = db.ann.df[c(
      "seqnames", "start", "end", "feature", "external_gene_name", "gene_biotype",
      "start_position", "end_position",
      "insideFeature", "distancetoFeature", "shortestDistance", "fromOverlappingOrNearest")]

    # merge bed and annotations
    ann.merged = merge(as.data.frame(db.gr), db.ann.df,
                       by.x = c("seqnames", "start", "end"),
                       by.y = c("seqnames", "start", "end"),
                       all = TRUE)

    # sort
    ann.merged$seqnames = as.character(ann.merged$seqnames)
    ann.merged = ann.merged[with(ann.merged, order(seqnames, start)), ]

    # generate file name
    th = format(th, nsmall = 2)
    contrast.name = paste0(group1, "-vs-", group2)
    if("Block1Val" %in% colnames(contrasts))
    {
      contrast.name = paste0(contrast.name, ".blocking")
    }
    filename_base = paste0(out.dir, "/diff_bind.", contrast.name, ".q", gsub("\\.", "", th))
    filename_csv = paste0(filename_base, ".csv")
    filename_bed = paste0(filename_base, ".bed")
    filename_up_bed = paste0(filename_base, ".up.bed")
    filename_dn_bed = paste0(filename_base, ".dn.bed")

    # save full table
    message("[generateDiffBindReport] save as: ", filename_csv)
    write.csv(ann.merged, file = filename_csv, row.names = FALSE)
    Sys.sleep(1)

    # save peaks in BED format
    bed_cols = c("seqnames", "start", "end")
    message("[generateDiffBindReport] save as: ", filename_bed)
    write.table(ann.merged[, bed_cols], file = filename_bed, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    Sys.sleep(1)
    message("[generateDiffBindReport] save as: ", filename_up_bed)
    write.table(subset(ann.merged, Fold > 0)[, bed_cols], file = filename_up_bed, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    Sys.sleep(1)
    message("[generateDiffBindReport] save as: ", filename_dn_bed)
    write.table(subset(ann.merged, Fold < 0)[, bed_cols], file = filename_dn_bed, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    Sys.sleep(1)
  }
}

message(" ========== load data ========== ")

# load data
sampleSheet = read.csv(sampleSheetCsv)
print(sampleSheet)
# multi-threading possible, but fails at dba.analyze() with larger datasets (>100,000 peaks)
# db_config = data.frame(AnalysisMethod = normMethodObj, th = 0.05, cores = 4)
db_config = data.frame(AnalysisMethod = normMethodObj, th = 0.05, RunParallel = FALSE)
db = dba(sampleSheet = sampleSheet, bCorPlot = FALSE, config = db_config)
print(db)

message(" ========== calculate binding matrix ========== ")

# calculate a binding matrix with scores based on read counts for every sample (affinity scores),
# rather than confidence scores for only those peaks called in a specific sample (occupancy scores)
db = dba.count(db, score = DBA_SCORE_RPKM)

# automatic contrasts
db = dba.contrast(db, categories = DBA_CONDITION, minMembers = 2)

# save
dbContrastRData = paste0(outDir, "/db.contrast.RData")
save(db, file = dbContrastRData)

message(" ========== generate plots ========== ")

# colors (combine multiple to prevent running out of colors)
colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))
colors = unique(colors)

# heatmap
heatmapPDF = paste0(outDir, "/plot.heatmap.fpkm.pdf")
pdf(heatmapPDF, family = "Palatino", pointsize = 10)
dba.plotHeatmap(db, score = DBA_SCORE_RPKM, colScheme = "Blues",
                colSideCols = colors, rowSideCols = colSideCols, cexRow = 0.8, cexCol = 0.8)
dev.off()

# pca based on conditions
pcaPDF = paste0(outDir, "/plot.pca.fpkm.pdf")
pdf(pcaPDF, family = "Palatino", pointsize = 10)
dba.plotPCA(db, DBA_CONDITION, score = DBA_SCORE_RPKM, label = DBA_ID, vColors = colors, labelSize = 0.6)
dev.off()

message(" ========== perform differential binding analysis ========== ")

# perform differential binding affinity analysis
db = dba.analyze(db, bFullLibrarySize = TRUE, bCorPlot = FALSE)

# save
dbAnalyzeRData = paste0(outDir, "/db.analyze.RData")
save(db, file = dbAnalyzeRData)

# show contrasts
contrasts = dba.show(db, bContrasts = TRUE)
contrasts = as.data.frame(contrasts)
print(contrasts)

message(" ========== retrieve biomart annotations ========== ")

# retrieve annotations
#if (genome == "hg19") {
#  martEns = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", verbose=F)
#}
#if (genome == "mm10") {
#  martEns = useMart(host="sep2017.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", verbose=F)
#}
# martEnsDF = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), mart=martEns)
# martEnsTSS = getAnnotation(mart=martEns, featureType="TSS")
martEnsDF = readRDS(paste(genomeDir,"martEnsDF.rds",sep='/'))
martEnsTSS = readRDS(paste(genomeDir,"martEnsTSS.rds",sep='/'))

message(" ========== generate reports ========== ")

# generate report for each possible contrast with different cutoffs
for (i in 1:length(row.names(contrasts))) {
  print(contrasts[i,])
  # using try to prevent "execution halted" that kills script if there are no significant results
  try(generateDiffBindReport(dba=db, contrast=i, th=1.00, method=normMethodObj, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
  try(generateDiffBindReport(dba=db, contrast=i, th=0.20, method=normMethodObj, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
  try(generateDiffBindReport(dba=db, contrast=i, th=0.05, method=normMethodObj, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
}

# repeat with blocking factor if blocking factor parameter was passed
if (!is.na(blockFactor)) {
  message(" ========== re-analyze with blocking factor ========== ")

  # automatic contrasts with blocking factor
  db = dba.contrast(db, categories = DBA_CONDITION, block = DBA_REPLICATE)
  db = dba.analyze(db, bFullLibrarySize = TRUE, bCorPlot = FALSE)
  dbAnalyzeBlockingRData = paste0(outDir, "/db.analyze.blocking.RData")
  save(db, file = dbAnalyzeBlockingRData)

  # show contrasts
  contrasts = dba.show(db, bContrasts = TRUE)
  contrasts = as.data.frame(contrasts)
  print(contrasts)

  message(" ========== generate reports with blocking factor ========== ")

  # generate report for each possible contrast
  for (i in 1:length(row.names(contrasts))) {
    print(contrasts[i,])
    # using try to prevent "execution halted" that kills script if there are no significant results
    try(generateDiffBindReport(dba=db, contrast=i, th=1.00, method=normMethodBlock, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
    try(generateDiffBindReport(dba=db, contrast=i, th=0.20, method=normMethodBlock, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
    try(generateDiffBindReport(dba=db, contrast=i, th=0.05, method=normMethodBlock, tss=martEnsTSS, mart.df=martEnsDF, reps=T, out.dir=outDir))
  }
}



# end
