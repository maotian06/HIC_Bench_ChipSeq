#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: chipseq-venn-diagram.tcsh OUTPUT-DIR PARAMETER-SCRIPT PEAKS-BRANCH [OBJECTS]
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


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------



# create merged peaks reference
set peaks = `echo $objects | tr ' ' '\n' | awk -v b=$branch '{print b"/"$1"/peaks.bed"}'`

# extract group level peak calling result
set group_branch = `echo $branch | sed -E 's/peaks.by_sample/peaks.by_group/'`
if ( $group_branch != "") then
  set group_name = `basename $outdir`
  set peaks = "$peaks $group_branch/$group_name/peaks.bed"
endif

scripts-send2err "-- CMD:"
scripts-send2err "Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir $peaks"
Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir $peaks

if ( $bed_files != "") then
  foreach file ($bed_files)
    set filename = `basename $file:r`.pdf
    scripts-send2err "Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir -f $filename $peaks $file"
    Rscript --vanilla code/chipseq-venn-diagram.r -o $outdir -f $filename $peaks $file
  end
endif
 
# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


