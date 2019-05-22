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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
scripts-send2err "Setting parameters..."
source $params
scripts-send2err "-- Parameters: "
scripts-send2err "- Extending promoter upstream and downstream by nt = $promoter_proximal"
scripts-send2err "- Genome = $genome"

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# filter out inputs
 if ($include_input == 'false') set objects = `echo $objects | tr ' ' '\n' | grep -vi input`

# determine input files
if ( $objects != "" ) then
  set peaks = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/peaks.bed"}'`

  scripts-send2err "Rscript --vanilla code/chipseq-peakanno.r -g $genome -d $promoter_proximal -o $outdir $peaks"
  Rscript --vanilla code/chipseq-peakanno.r -g $genome -d $promoter_proximal -o $outdir $peaks
endif

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



