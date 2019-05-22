#!/bin/tcsh
set echo

source ./code/code.main/custom-tcshrc         # shell settings (must be included in all scripts)

##
## USAGE: chipseq-homer_tagdir.tcsh OUTPUT-DIR PARAMETER-SCRIPT ALIGNMENT-BRANCH OBJECT(S)
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params
# scripts-send2err "-- Parameters: "
# scripts-send2err "- macs = $macs_params"
# scripts-send2err "- annotation = $annot_params"

# create path
scripts-create-path $outdir/

# set parameters
source $params
if (! $?NSLOTS) then
  set threads = 1
else
  set threads = $NSLOTS
endif



# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------
echo "outdir is $outdir"
echo "params is $params"
echo "branch is $branch"
echo "objects is $objects"
echo "genome is $genome"
echo "genome_dir is $genome_dir"
echo "threads is $threads"
set var = `pwd`   

echo "pwd is $var"

# echo "${branch}/${objects}/alignments.bam"

# http://homer.salk.edu/homer/ngs/annotation.html
# examples:
# makeTagDirectory $outdir/ "${branch}/${objects}/alignments.bam"
# findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
# findMotifsGenome.pl ${branch}/${objects}/peaks.bed $genome $outdir/ -p $threads -preparsedDir $outdir/preparsed $motif_options
# annotatePeaks.pl <peak/BED file> <genome>   > <output file>

annotatePeaks.pl "${branch}/${objects}/peaks.bed" "$genome" -annStats $outdir/annotation_stats.txt -go $outdir/gene_ontolgy -genomeOntology $outdir/genome_ontolgy > $outdir/annotated_peaks.txt


# convert to bed format
pos2bed.pl $outdir/annotated_peaks.txt > $outdir/annotated_peaks.bed

# without gene genome ontology:
# annotatePeaks.pl "${branch}/${objects}/peaks.bed" "$genome" -annStats $outdir/annotation_stats.txt > $outdir/annotated_peaks.txt



# get the number of types of peaks
set Peak_Type_Stats = $outdir/peak_type_stats.txt
tail -n +2 $outdir/annotated_peaks.txt | cut -f8 | cut -d '(' -f1 | sort | uniq -c | sed -e 's/ *//' -e 's/ /\t/' -e "s/'//" -e 's/NA/Other/' -e 's/ //' > $Peak_Type_Stats


# plot the number of peaks per type
module unload r; module load r/3.3.0

Rscript --vanilla - $outdir $Peak_Type_Stats <<HERE
  ## R code
  cat("\nR loaded\n")
  # get R script args
  args <- commandArgs(TRUE)
  cat("Script args are:\n")
  args

  OutDir<-args[1]
  PeakTypesTable_file<-args[2]
  
  #quit()
  
  # read in the table from the file
  PeakTypesTable<-read.table(file = PeakTypesTable_file,header = FALSE,sep = "\t")
  PeakTypesTable

  # format the table
  rownames(PeakTypesTable)<-PeakTypesTable[["V2"]]
  PeakTypesTable<-PeakTypesTable["V1"]
  colnames(PeakTypesTable)<-"NumberPeaks"
  PeakTypesTable
  
  # get the total peaks
  TotalPeaks<-sum(PeakTypesTable[["NumberPeaks"]])
  
  # Calculate the percent 
  PeakTypesTable[["PercentPeaks"]]<-signif((PeakTypesTable[["NumberPeaks"]]/TotalPeaks)*100,digits=3)
  
  # save the values to be plotted into a transposed matrix, since thats what the barplot() likes
  Peaks_Raw_Matrix<-t(as.matrix(PeakTypesTable["NumberPeaks"]))
  Peaks_Pcnt_Matrix<-t(as.matrix(PeakTypesTable["PercentPeaks"])) # don't actually need this..
  
  
  # make a barplot
  # # default par: par() mar # [1] 5.1 4.1 4.1 2.1
  pdf(file = paste0(OutDir,"/peaks_type_barplots.pdf"),width = 8,height = 8)
  par(mar=c(5.1, 7.1, 4.1, 2.1))
  barplot(Peaks_Raw_Matrix,main="Number of peaks",horiz = TRUE,las=1,cex.names=1)
  dev.off()
  
  # write a CSV of the final table
  # # peel off the rownames into a separate vector
  Types<-row.names(PeakTypesTable)
  write.csv(x = cbind(Types,PeakTypesTable), file = paste0(OutDir,"/peaks_type_table.csv"),quote = FALSE,row.names = FALSE)
  save.image(file = paste0(OutDir,"/R_session.Rdata") )


HERE

# https://github.com/stevekm/Bioinformatics/blob/88ab3a029cd99978b33676bcae3c2302574053bd/HOMER_motif_analysis/peak_types_plot.R


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



