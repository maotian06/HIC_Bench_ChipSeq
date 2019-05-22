#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: chipseq-summ-peakanno.tcsh OUTPUT-DIR PARAMETER-SCRIPT PEAKS-BRANCH [OBJECTS]
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
scripts-send2err "-- Parameters: "
scripts-send2err "- Extending promoter upstream and downstream by nt = $promoter_proximal"
scripts-send2err "- Genome = $genome"


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# filter out inputs
if ($include_input == 'false') set objects = `echo $objects | tr ' ' '\n' | grep -vi input`

# create merged peaks reference
set peaks = `echo $objects | tr ' ' '\n' | awk -v b=$branch '{print b"/"$1"/peaks.bed"}'`

# plot annotations

scripts-send2err "Rscript --vanilla code/chipseq-summ-peakanno.r -g $genome -d $promoter_proximal -o $outdir $peaks"
Rscript --vanilla code/chipseq-summ-peakanno.r -g $genome -d $promoter_proximal -o $outdir $peaks

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


