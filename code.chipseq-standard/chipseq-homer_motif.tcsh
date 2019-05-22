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

# http://homer.salk.edu/homer/ngs/peakMotifs.html
# examples:
# makeTagDirectory $outdir/ "${branch}/${objects}/alignments.bam"
# findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]

findMotifsGenome.pl ${branch}/${objects}/peaks.bed $genome $outdir/ -p $threads -preparsedDir $outdir/preparsed -dumpFasta $motif_options
# delete the preparsed
rm -rf $outdir/preparsed

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



