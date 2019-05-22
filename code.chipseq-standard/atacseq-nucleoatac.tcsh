#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # shell settings (must be included in all scripts)

##
## USAGE: create-matrix.tcsh OUTPUT-DIR PARAMETER-SCRIPT PEAKS-BRANCH OBJECT(S)
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
set peak_branch = $branch
source ./code/code.main/scripts-read-job-vars $peak_branch "$objects" "genome genome_dir branch"
set aln_branch = $branch
set branch = $peak_branch

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- shift_dist = $shift_dist"
scripts-send2err "- nucleoatac = $nucleoatac"
scripts-send2err "- nucleoatac_path = $nucleoatac_path"

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
set alignments = `echo $samples | tr ' ' '\n' | awk -v d=$aln_branch '{print d"/"$0"/alignments.bam"}'`
set ref_regions = $branch/*/peaks.bed 

scripts-create-path $outdir/
# create padded broad peaks bed file
set padded_regions = $outdir/padded_peaks.bed
cat $branch/*/peaks.bed | sort -k1,1 -k2,2n | gtools-regions shiftp -5p -$shift_dist -3p $shift_dist | gtools-regions link > $outdir/padded_peaks.bed


# create matrix
set threads = $#alignments
set alignments = `echo $alignments | tr ' ' ','`
$nucleoatac_path run --cores $threads --fasta $genome_fa --bam $alignments --bed $padded_regions --out $outdir/nucleoatac

# cleanup
#rm -f $ref

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


