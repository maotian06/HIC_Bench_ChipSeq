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


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------
echo "outdir is $outdir"
echo "params is $params"
echo "branch is $branch"
echo "objects is $objects"
set var = `pwd`   
echo "pwd is $var"

# some notes:
# path to the input file
# echo "${branch}/${objects}/alignments.bam"
# makeTagDirectory $outdir/ "${branch}/${objects}/alignments.bam"
# findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
# i.e. findPeaks ERalpha-ChIP-Seq/ -style factor -o auto -i Control-ChIP-Seq/


# determine input control files
# # control tagdir entry empty by default
set control_tagdir = 
# check if params had use input true
if ($use_input != 'false') then
  # get the control sample from the sample sheet
  set control_samples = `./code/read-sample-sheet.tcsh $sheet "$objects" control | sort -u | grep -v '^NA$'`
  # check if the control sample entry is not empty
  if ("$control_samples" != "") then
    # get the branch for the control tagdir to use
    set control_tagdir = `echo $control_samples | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0}'`
    set control_tagdir = "-i $control_tagdir"
  endif
endif

echo "control_tagdir is $control_tagdir"
echo "outdir/HOMER_superenhancer.txt is ${outdir}/HOMER_superenhancer.txt"
echo "outdir/HOMER_typical_enhancer.txt is ${outdir}/HOMER_typical_enhancer.txt"

# run findPeaks
set my_command = "findPeaks ${branch}/${objects} $control_tagdir -style $peak_style  -o ${outdir}/peaks.txt $extra_params"
echo "Running the HOMER findPeaks command:"
$my_command

# convert to BED format
pos2bed.pl "${outdir}/peaks.txt" > "${outdir}/peaks.bed"

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."



