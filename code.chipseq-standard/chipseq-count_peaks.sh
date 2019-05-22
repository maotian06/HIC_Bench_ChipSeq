#!/bin/bash

## USAGE: chipseq-count_peaks.sh /path/to/outdir /path/to/dir_with_BED_files 
## DESCRIPTION: This script will find all "peaks.bed" files in a parent directory
## and output a single tab-separated table with the number of peaks per sample


#~~~~~~ GET SCRIPT ARGS ~~~~~~~~~~~~#
outdir="$1"
# make sure outdir exists
mkdir -p "$outdir"

peaks_dir="$2" # the pipeline branch


#~~~~~~ SET OUTPUTS ~~~~~~~~~~~~#
output_table="${outdir}/peaks_stats.tsv"
echo -e 'Peaks\tSample' > "$output_table"

#~~~~~~ FIND PEAKS ~~~~~~~~~~~~#
peaks_files=$(find $peaks_dir -type f -name peaks.bed | LC_ALL=C sort)

#~~~~~~ COUNT PEAKS ~~~~~~~~~~~~#

echo "Counting pipeline peaks"
echo -e "Input directory is:\n${peaks_dir}\n"
echo -e "Output Directory is:\n${outdir}\n"
echo -e "Output table is:\n${output_table}\n"


for i in $peaks_files; do
echo -e "Found peak file:\n${i}\n"

# get the sample ID from the dirname
tmp_sampleID=$(basename $(dirname "$i"))
echo "Sample is $tmp_sampleID"

# get the number of peaks 
num_peaks=$(cat $i | wc -l)
echo "num_peaks is $num_peaks"

# print it to the table
echo -e "${num_peaks}\t${tmp_sampleID}" >> "$output_table"

done


#~~~~~~ PLOT PEAKS ~~~~~~~~~~~~#
# pass the table file to R for plotting
  Rscript --slave --no-save --no-restore - "$output_table" "$peaks_dir" "$outdir" <<EOF
  #!/usr/bin/env Rscript
  # R code for plotting the peaks table
  
  # ~~~~~ GET SCRIPT ARGS ~~~~~~~ #
  args <- commandArgs(TRUE); cat("Script args are:\n"); args
  
  # path to the peaks table file
  peaks_table_file <- args[1]
  
  # get the sample branch
  # peaks_branch <- basename(dirname(peaks_table_file))
  peaks_branch <- args[2]
  
  outdir <- args[3]
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
  
  # load the file into a dataframe
  peaks_table_df<-read.table(peaks_table_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  # convert the Sample column entries into rownames
  rownames(peaks_table_df) <- peaks_table_df[["Sample"]]
  
  
  # plot layout setup
  plot_layout_matrix<-structure(c(1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 
                                2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 
                                2L, 2L, 2L), .Dim = c(8L, 4L), .Dimnames = list(NULL, c("V1", 
                                                                                        "V2", "V3", "V4")))
  pdf(file = file.path(outdir, "peaks_count_barplots.pdf"),width = 8,height = 8)
  # setup the panel layout
  layout(plot_layout_matrix)
  # need to set this for some reason
  # plot margins # c(bottom, left, top, right) # default is c(5, 4, 4, 2) + 0.1
  par(mar=c(0,0,5,0))
  # call blank plot to fill the first panel
  plot(1,type='n',axes=FALSE,xlab="",ylab="",main = paste0("Peaks\n",peaks_branch),cex.main=1.5) 
  # set up the Legend in the first panel
  # legend("bottom",legend=colnames(overlap_df),bty = "n",cex=1.0) # fill=BARPLOT_COLORS,,ncol=length(BARPLOT_COLORS)
  # set some plot margin parameters to fit the names
  par(mar=c(5,14,0,2)+ 0.1) 
  barplot(t(peaks_table_df),
        # main=peaks_branch,
        cex.names = 0.7,
        horiz = T,
        # col=BARPLOT_COLORS,
        border=NA,
        las=1,
        # cex.names=Names_scale,
        xlab="Number of peaks",
        space=0.6
  ) 
  dev.off()
EOF
