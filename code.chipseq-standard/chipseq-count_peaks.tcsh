#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: chipseq-count_peaks.tcsh OUTPUT-DIR PARAM-SCRIPT PEAKS-BRANCH [OBJECTS]
##

# process command-line inputs
if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

# inputs
set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if samples is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

echo "outdir is $outdir"
echo "params is $params"
echo "branch is $branch" # input directory
echo "objects is $objects" # samples to be processed


./code/chipseq-count_peaks.sh "$outdir" "$branch" 

# Rscript --vanilla ./code/chipseq-diffbind.R $outdir $diffbind_sample_sheet $genome $diffbind_blocking_factor


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



