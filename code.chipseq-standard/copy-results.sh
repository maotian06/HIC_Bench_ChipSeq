#!/bin/bash

##
## USAGE: scripts-copy-chipseq.sh /path/to/outdir /path/to/project_dir 
## this script copies the results of a chip-seq analysis over to a Results dir for a client. 
## 

#~~~ Shell Options ~~~~#
# remember whether extglob was originally set, so we know whether to unset it
shopt -q extglob; extglob_set=$?
# set extglob if it wasn't originally set.
((extglob_set)) && shopt -s extglob
# Note, 0 (true) from shopt -q is "false" in a math context.

shopt -q nullglob; nullglob_set=$?
((nullglob_set)) && shopt -s nullglob

shopt -q globstar; globstar_set=$?
((globstar_set)) && shopt -s globstar

#~~~~~~~~~~~~~~~~~~#

# make sure that the correct number of script arguments were used
if [ $# != 2 ] # if not enough args provided
then
  grep '^##' $0 # print out lines from the script that start with '##' and exit
  exit
fi



OUT_DIR="$(readlink -m $1)" # get full path, through symlinks
PROJ_DIR="$(readlink -m $2)"
# Algn_stats_dir="$(readlink -m $3)"
# FastQC_dir="$(readlink -m $3)"
Pipeline_dir="$PROJ_DIR"/pipeline


# make sure the OUT_DIR exists
mkdir -p $OUT_DIR


#
##
### COPY THE BIG WIGS
##
# check if the $PROJ_DIR has peaks results dir
if [[ -d $Pipeline_dir/align/results ]]; then
  Pipeline_align_dir="$Pipeline_dir/align/results"

  # make outdir dir
  OUT_DIR_align="$OUT_DIR/bigwig"
  mkdir -p $OUT_DIR_align

  # glob the bigwig files
  BIG_WIGS=($Pipeline_align_dir/**/*.bw)

  for i in "${BIG_WIGS[@]}"; do
    # get the sample name from its dir name
    q=$(basename $(dirname "$i" ) )

    # strip out things we don't need in the name
    z=${q/*./}

    # copy the file with new name (don't overwrite if already present!)
    echo "${z}.bw"
    rsync -rt "$i" "$OUT_DIR_align"/"${z}.bw"
  done
else
  echo "BIG WIG RESULTS NOT FOUND"
fi
###
##
#

#
##
### COPY THE BAMS
# check if the $PROJ_DIR has peaks results dir
if [[ -d $Pipeline_dir/align/results ]]; then
  Pipeline_align_dir="$Pipeline_dir/align/results"

  # make outdir dir
  OUT_DIR_align="$OUT_DIR/bam"
  mkdir -p $OUT_DIR_align

  # glob the bam files
  BAM_list=($Pipeline_align_dir/**/*.bam)

  for i in "${BAM_list[@]}"; do
    # get the sample name from its dir name
    q=$(basename $(dirname "$i" ) )
    # strip out things we don't need in the name
    z=${q/*./}
    echo "${z}.bam"
    rsync -rt "$i" "$OUT_DIR_align"/"${z}.bam"
    rsync -rt "${i}.bai" "$OUT_DIR_align"/"${z}.bam.bai"
  done
else
  echo "BAM RESULTS NOT FOUND"
fi
###
##
#

#
##
### COPY THE PEAKS
if [[ -d $Pipeline_dir/peaks/results ]]; then
  Pipeline_peaks_dir="$Pipeline_dir/peaks/results"

  # make the outdir
  OUT_DIR_peaks=$OUT_DIR/peaks
  mkdir -p $OUT_DIR_peaks

  # glob the peaks files we want
  FILES=($Pipeline_peaks_dir/**/@(macs_peaks.xls|peaks.bed|peak-scores.bed|macs_peaks.narrowPeak|macs_summits.bed|peaks.tsv) )

  for i in "${FILES[@]}"; do
    # flatten the filepath to make the filename # consider changing this later ?
    TMP_NAME=$(echo ${i#*results/} | tr / _ )
    echo "$OUT_DIR_peaks/$TMP_NAME"
    rsync -rt $i $OUT_DIR_peaks/$TMP_NAME
  done
fi
#


#
##
### copy the PCA
if [ -d $Pipeline_dir/pca/results ]; then
  Pipeline_pca_dir="$Pipeline_dir/pca/results"
  # make the outdir
  OUT_DIR_pca="$OUT_DIR/pca"
  mkdir -p $OUT_DIR/pca

  # glob the PCA files we want
  FILES=($Pipeline_pca_dir/**/@(report.raw.pdf|report.qnorm.pdf|report.mnorm.pdf|labels.tsv|report.mnorm.pdf) )

  for i in "${FILES[@]}"; do
    # flatten the filepath to make the filename # consider changing this later ?
    TMP_NAME=$(echo ${i#*pca.standard/} | tr / _ )
    echo "$OUT_DIR_pca/$TMP_NAME"
    rsync -rt $i $OUT_DIR_pca/$TMP_NAME
  done
fi


# Copy heatmaps
if [ -d $Pipeline_dir/heatmaps/results ]; then
  Pipeline_heatmap_dir="$Pipeline_dir/heatmaps/results"
  # make the outdir
  OUT_DIR_heatmap="$OUT_DIR/heatmaps"
  mkdir -p $OUT_DIR_heatmap

  # glob the PCA files we want
  FILES=($Pipeline_heatmap_dir/**/@(clustering.tif) )
  for i in "${FILES[@]}"; do
    # flatten the filepath to make the filename # consider changing this later ?
    # strip this : align.by_sample.bowtie2/all-samples/
    TMP_NAME=$(echo ${i//align.by_sample.bowtie2\/})
    TMP_NAME=$(echo ${TMP_NAME//all-samples\/})
    TMP_NAME=$(echo $TMP_NAME | tr / _ )
    echo "$OUT_DIR_heatmap/$TMP_NAME"
    rsync -rt $i $OUT_DIR_heatmap/$TMP_NAME
  done
fi




# Copy qc fingerprints; chip-fingerprint.pdf
if [ -d $Pipeline_dir/qc/results ]; then
  Pipeline_qc_dir="$Pipeline_dir/qc/results"
  # make the outdir
  OUT_DIR_qc="$OUT_DIR/qc"
  mkdir -p $OUT_DIR_qc

  # glob the PCA files we want
  FILES=($Pipeline_qc_dir/**/@(chip-fingerprint.pdf) )

  for i in "${FILES[@]}"; do
    # flatten the filepath to make the filename # consider changing this later ?
    TMP_NAME=$(echo ${i#*align.by_sample.bowtie2\/} | tr / _ )
    # strip this : align.by_sample.bowtie2/all-samples/
    # TMP_NAME=$(echo ${i//align.by_sample.bowtie2\/} | tr / _ )
    # TMP_NAME=$(echo ${TMP_NAME//all-samples\/})
    # TMP_NAME=$(echo $TMP_NAME | tr / _ )
    rsync -rt $i $OUT_DIR_qc/$TMP_NAME
  done
fi


#
##
# Copy the fastqc results
if [ -d $Pipeline_dir/fastqc/results ]; then
  Pipeline_fastqc_dir="$Pipeline_dir/fastqc/results"
  # make the outdir
  OUT_DIR_fastqc="$OUT_DIR/fastqc"
  mkdir -p $OUT_DIR_fastqc
  rsync -rt --exclude "job.*" --exclude ".db*" "$Pipeline_fastqc_dir/" "$OUT_DIR_fastqc"
fi
# if [ ! -d $OUT_DIR/$(basename $FastQC_dir) ]; then
#   cp -avn "$FastQC_dir" $OUT_DIR/$(basename $FastQC_dir)
# fi
# Copy qc fiongerprints; chip-fingerprint.pdf




#  Skip this for now.. update this code later too!
# # copy the Peaks table
# if [ -d $Pipeline_dir/peaktable/results ]; then
#   Pipeline_peakstabl_dir=$Pipeline_dir/peaktable/results
#   
#   # make the outdir
#   $OUT_DIR_peaktable="$OUT_DIR/peaktable"
#   mkdir -p $OUT_DIR_peaktable
#   
#   FILES=$(find $Pipeline_peakstabl_dir -type f \( ! -name "*.win" ! -name "job*" ! -name "obj.*" ! -name "*.sh" ! -name "*.branch" \) )
#   for i in $FILES; do
#     # use the filepath as the dirname;
#     TMP_NAME=$(echo ${i#*peaktable.standard/} | tr / _ )
#     # echo $TMP_NAME
#     cp -avn $i $Pipeline_peakstabl_dir/$TMP_NAME
#   done
# 
#   # rsync -a --exclude "job*" --exclude "ref*" $PROJ_DIR/peaktable/results/peaktable.standard/ $OUT_DIR/peaktable
# fi
# ###
# ##
# #


# Copy the alignment summary stats, if it doesn't exist
# if [ ! -d $OUT_DIR/$(basename $Algn_stats_dir) ]; then
#   cp -avn "$Algn_stats_dir" $OUT_DIR/$(basename $Algn_stats_dir)
# fi
if [ -d "$Pipeline_dir/align-stats/results" ]; then
  Pipeline_align_stats_dir="$Pipeline_dir/align-stats/results"
  # make the outdir
  OUT_DIR_algnstats="$OUT_DIR/align-stats"
  mkdir -p "$OUT_DIR_algnstats"

  rsync -rt --exclude "job.*" "$Pipeline_dir/align-stats/results/" "$OUT_DIR_algnstats"
fi


# copy the diffbind results diffbind
if [ -d $Pipeline_dir/diffbind/results ]; then
  if [ ! -d $OUT_DIR/diffbind ]; then
    Pipeline_diffbind_dir="$Pipeline_dir/diffbind/results"
    # make the outdir
    OUT_DIR_diffbind="$OUT_DIR/diffbind"
    mkdir -p $OUT_DIR_diffbind
    rsync -rt --exclude "job.*" --exclude ".db*" "$Pipeline_diffbind_dir/" "$OUT_DIR_diffbind"
  fi
fi


# Copy the sample sheet
rsync -rt "$Pipeline_dir/inputs/sample-sheet.tsv" $OUT_DIR/sample-sheet.tsv


# copy the report
if [ ! -f $OUT_DIR/chipseq_report.pdf ]; then
  rsync -rt "$PROJ_DIR/report/chipseq_report.pdf" $OUT_DIR/chipseq_report.pdf
fi

# unset globs if it wasn't originally set
((extglob_set)) && shopt -u extglob
((nullglob_set)) && shopt -u nullglob
((globstar_set)) && shopt -u globstar


# set outdir permissions to allow read/write access to group members
# chmod g+r -R $OUT_DIR*

