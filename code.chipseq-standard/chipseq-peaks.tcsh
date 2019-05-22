#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # shell settings (must be included in all scripts)

##
## USAGE: chipseq-peaks.tcsh OUTPUT-DIR PARAMETER-SCRIPT ALIGNMENT-BRANCH OBJECT(S)
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
set caller = `echo $params | cut -d '.' -f2 | cut -d '_' -f1`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- caller = $caller"
scripts-send2err "- caller_params = $caller_params"
scripts-send2err "- annotation = $annot_params"

if (! $?NSLOTS) then
  set threads = 8
else
  set threads = $NSLOTS
endif

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# determine input files
set treatment_aln = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`

if ("$caller" == 'sicer') then
  set caller = "SICER-rb.sh"
endif

mkdir $outdir/treatment_preprocess
if ( $#treatment_aln == 1 ) then
  scripts-send2err "Linking bam files, command:"
  scripts-send2err "ln -sf `readlink -f $treatment_aln` $outdir/treatment_preprocess/merged.bam"
  ln -sf `readlink -f $treatment_aln` $outdir/treatment_preprocess/merged.bam
else
  scripts-send2err "Merging bam files, command:"
  scripts-send2err "samtools merge -@ $threads $outdir/treatment_preprocess/merged.bam $treatment_aln"
  samtools merge -@ $threads $outdir/treatment_preprocess/merged.bam $treatment_aln
endif
if ("$caller" == 'SICER-rb.sh') then
  scripts-send2err "Converting bam files to BED, command:"
  scripts-send2err "bamToBed -i $outdir/treatment_preprocess/merged.bam \> $outdir/treatment_preprocess/merged.bed"
  bamToBed -i $outdir/treatment_preprocess/merged.bam > $outdir/treatment_preprocess/merged.bed
  set treatment_aln_ps = $outdir/treatment_preprocess/merged.bed
  set treatment_bed = merged.bed
else if ("$caller" == 'homer') then
  scripts-send2err "Creating TagDirectory, command:"
  scripts-send2err "makeTagDirectory $outdir/treatment_preprocess/ $outdir/treatment_preprocess/merged.bam"
  makeTagDirectory $outdir/treatment_preprocess/ $outdir/treatment_preprocess/merged.bam
  set treatment_aln_ps = $outdir/treatment_preprocess
else
  set treatment_aln_ps = $outdir/treatment_preprocess/merged.bam
endif

set control_aln = 
set control_aln_ps = 
set control_bed = 
if ($use_input != 'false') then
  set control_samples = `./code/read-sample-sheet.tcsh $sheet "$objects" control | sort -u | grep -v '^NA$'`
  if ("$control_samples" != "") then
    set control_aln = `echo $control_samples | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`
    mkdir $outdir/control_preprocess
    if ( $#control_aln == 1 ) then
      scripts-send2err "Linking bam files, command:"
      scripts-send2err "ln -sf `readlink -f $control_aln` $outdir/control_preprocess/merged.bam"
      ln -sf `readlink -f $control_aln` $outdir/control_preprocess/merged.bam
    else
      scripts-send2err "Merging bam files, command:"
      scripts-send2err "samtools merge -@ $threads $outdir/control_preprocess/merged.bam $control_aln"
      samtools merge -@ $threads $outdir/control_preprocess/merged.bam $control_aln
    endif
    if ("$caller" == 'SICER-rb.sh') then
      scripts-send2err "Converting bam files to BED, command:"
      scripts-send2err "bamToBed -i $outdir/control_preprocess/merged.bam \> $outdir/control_preprocess/merged.bed"
      bamToBed -i $outdir/control_preprocess/merged.bam > $outdir/control_preprocess/merged.bed
      ln -sf ../control_preprocess/merged.bed $outdir/treatment_preprocess/control_merged.bed
      set control_aln_ps = $outdir/treatment_preprocess/control_merged.bed
      set control_bed = control_merged.bed
      set caller = "SICER.sh"
    else if ("$caller" == 'homer') then
      scripts-send2err "Creating TagDirectory, command:"
      scripts-send2err "makeTagDirectory $outdir/control_preprocess/ $outdir/control_preprocess/merged.bam"
      makeTagDirectory $outdir/control_preprocess/ $outdir/control_preprocess/merged.bam
      set control_aln_ps = "-i $outdir/control_preprocess"
     else if ("$caller" == 'macs' ) then
      set control_aln_ps = "-c $outdir/control_preprocess/merged.bam"
     else
      set control_aln_ps = "$outdir/control_preprocess/merged.bam"
    endif
  endif
endif

if ("$caller" == 'macs') then
# run macs2
set g = `echo $genome | sed 's/[0-9]\+$//'`
set gsize = `grep "^$g	" params/gsize.info.txt | cut -f2`
if ($gsize == "") then
  scripts-send2err "Error: genome effective size not found in params/gsize.info.txt!"
  exit 1
endif
scripts-send2err "- genome = $genome"
scripts-send2err "Running the macs2 callpeak command:"
scripts-send2err "macs2 callpeak -t $treatment_aln_ps $control_aln_ps --outdir=$outdir --name=macs $caller_params -g $gsize"
macs2 callpeak -t $treatment_aln_ps $control_aln_ps --outdir=$outdir --name=macs $caller_params -g $gsize

# create standardized output files
scripts-send2err "Standardizing output"
set ext = narrowPeak
if (`echo $caller_params | grep -cw '\--broad'` == 1) set ext = broadPeak
cat $outdir/macs_peaks.$ext | cut -f-3,7 >! $outdir/peak-scores.bed
(echo "PEAK-ID\tCHROMOSOME\tSTART\tEND\tFOLD-CHANGE\tPVALUE(-log10)\tQVALUE(-log10)"; cat $outdir/macs_peaks.$ext | sort -k7,7rg | cut -f-4,7- | tools-cols -t 3 0 1 2 4 5 6) >! $outdir/peaks.table.tsv


else if ("$caller" == 'homer') then
# run homer
scripts-send2err "Running the HOMER findPeaks command:"
scripts-send2err "findPeaks $treatment_aln_ps $control_aln_ps -style $caller_params -o ${outdir}/peaks.txt $extra_params"
findPeaks $treatment_aln_ps $control_aln_ps -style $caller_params -o ${outdir}/peaks.txt $extra_params
scripts-send2err "Standardizing output"
(echo "PEAK-ID\tCHROMOSOME\tSTART\tEND\tFOLD-CHANGE\tPVALUE(-log10)\tQVALUE(-log10)"; cat $outdir/peaks.txt | grep -v "^#" | cut -f-4,11,12 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"(0-log($6)/log(10))}' | tools-cols -t 0 1 2 3 4 5 5) >! $outdir/peaks.table.tsv



else if ("$caller" =~ SICER*) then
# run sicer
scripts-send2err "Running the SICER command:"
set CWD = `pwd`
scripts-send2err "cd `dirname $treatment_aln_ps`"
scripts-send2err "$CWD/code/SICER/$caller ./ $treatment_bed $control_bed ./ $genome $caller_params"
cd `dirname $treatment_aln_ps`
$CWD/code/SICER/$caller ./ $treatment_bed $control_bed ./ $genome $caller_params
cd $CWD
scripts-send2err "cd $CWD"
scripts-send2err "Standardizing output"
if ("$control_aln_ps" != "") then
  cat `dirname $treatment_aln_ps`/*-islands-summary-* >! $outdir/peaks.txt
else 
  cat `dirname $treatment_aln_ps`/*.scoreisland | awk '{print $1"\t"$2"\t"$3"\tNA\tNA\tNA\t"$4"\tNA"}' >! $outdir/peaks.txt
endif
set t = `mktemp`
cut -f1,3 $genome_dir/genome.bed >! $t
code/wigToBigWig `dirname $treatment_aln_ps | perl -lane 'system "ls $_/*.wig"' | grep -v islandfiltered` $t $outdir/track.bw
rm -rf $t
mv `dirname $treatment_aln_ps`/*.scoreisland $outdir/gap_size.info
(echo "PEAK-ID\tCHROMOSOME\tSTART\tEND\tFOLD-CHANGE\tPVALUE(-log10)\tQVALUE(-log10)"; cat $outdir/peaks.txt | grep -v "^#" | cut -f-3,6,7,8 | tools-cols -t 0 1 2 4 3 5 | tools-rows -number -pref sicer_peak_ | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"(0-log($6)/log(10))"\t"(0-log($7)/log(10))}') >! $outdir/peaks.table.tsv

# run genome_scan
else
scripts-send2err "Scanning Genome"
gtools-regions win $caller_params $ref | cut -f1-3 | awk '{print $0"\t"$1":"$2"-"$3}' >! $outdir/ref.bed
set n_inputs = `echo $treatment_aln_ps $control_aln_ps | wc -w`
scripts-send2err "gtools-threaded matrix -p $n_inputs -i -nbins 1 --overlap-op hits -rpkm -format '%.0f' $treatment_aln_ps,$control_aln_ps $outdir/ref.bed \>! $outdir/matrix.tsv"
gtools-threaded matrix -p $n_inputs -i -nbins 1 --overlap-op hits -rpkm $treatment_aln_ps,$control_aln_ps $outdir/ref.bed >! $outdir/matrix.tsv
( echo "PEAK-ID\tCHROMOSOME\tSTART\tEND\tFOLD-CHANGE\tPVALUE(-log10)\tQVALUE(-log10)" ; cat $outdir/matrix.tsv | scripts-skipn 1 | sed -r 's/\t/ /2' | tools-vectors div | tr -d ' ' | tr ':-' '\t' | awk '$4>1&&$4!~/nan|inf/{print}' | gtools-regions link --label-func max | cut -f1-4 | awk -v cutoff="$min_fc" '$4>cutoff{print $1":"$2"-"$3"\t"$1"\t"$2"\t"$3"\t"$4"\tNA\tNA\tNA"}' ) >! $outdir/peaks.table.tsv
endif

# annotate peaks
scripts-send2err "Annotate peaks"
cat $outdir/peaks.table.tsv | scripts-skipn 1 | tools-cols -t 1 2 3 0 4 | sort -u | sort -k1,1 -k2,2g | awk '{print $0"\t+"}' >! $outdir/peaks.bed
(echo "MERGED-PEAK-ID LOCUS PEAK-SIZE GENE-SYMBOL REGION DISTANCE-START DISTANCE-END" | tr ' ' '\t'; cat $outdir/peaks.bed | gtools-overlaps $annot_params | cut -f-7 | sort -u) >! $outdir/peaks.annotated.tsv
scripts-join-tables $outdir/peaks.table.tsv $outdir/peaks.annotated.tsv >! $outdir/peaks.tsv

# cleanup
rm -f $outdir/peaks.table.tsv $outdir/peaks.annotated.tsv
#rm -rf $outdir/{control,treatment}_preprocess


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



