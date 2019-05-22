#!/bin/bash

##
## USAGE: setup-sample-files.sh <fastq source dir>
##
## FUNCTION: will scan given directory and generate proper sample structure in "./inputs/fastq/"
##


#########################


# check for correct number of arguments
if [ ! $# == 1 ] ; then
	echo -e "\n ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	exit 1
fi


#########################


# source directory
SOURCE_FASTQ_DIR=$1
SOURCE_FASTQ_DIR=$(readlink -f "$SOURCE_FASTQ_DIR")

# check that it exists
if [ ! -d $SOURCE_FASTQ_DIR ] || [ ! $SOURCE_FASTQ_DIR ]
then
	echo -e "\n ERROR! ERROR! FASTQ DIR $SOURCE_FASTQ_DIR DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


PROJECT_DIR=$(readlink -f .)
INPUTS_DIR="${PROJECT_DIR}/inputs"

# check that inputs directory exists
if [ ! -d $INPUTS_DIR ] || [ ! $INPUTS_DIR ]
then
	echo -e "\n ERROR! INPUTS DIR $INPUTS_DIR DOES NOT EXIST (THIS SHOULD BE RUN FROM ANALYSIS DIR) \n" >&2
	exit 1
fi

FASTQ_DEST_DIR="${INPUTS_DIR}/fastq"


#########################


# find FASTQs in source directory and create symlinks in ./inputs/fastqs
for SOURCE_FASTQ in $(find -L "$SOURCE_FASTQ_DIR" -type f -name "*.fastq.gz" | sort) ; do

	# for new fastq, remove flow cell index from the original filename (end it at R1/R2)
	DEST_FASTQ_BASE=$(basename "$SOURCE_FASTQ")
	DEST_FASTQ_BASE=${DEST_FASTQ_BASE/_R1_001.fastq.gz/_R1.fastq.gz}
	DEST_FASTQ_BASE=${DEST_FASTQ_BASE/_R2_001.fastq.gz/_R2.fastq.gz}

	# determine sample name (multiple steps for easy commenting and reading)
	SAMPLE=$DEST_FASTQ_BASE
	# 2 barcodes
	SAMPLE=$(echo "$SAMPLE" | perl -pe 's/_[ACTG]{6,}-[ACTG]{6,}_L00[0-9]_R[12].*//')
	# 1 barcode
	SAMPLE=$(echo "$SAMPLE" | perl -pe 's/_[ACTG]{4,}_L00[0-9]_R[12].*//')
	# no bar codes
	SAMPLE=$(echo "$SAMPLE" | perl -pe 's/_L00[0-9]_R[12].*//')
	# no bar codes or lane
	SAMPLE=$(echo "$SAMPLE" | perl -pe 's/_R[12].fastq.gz//')
	# improper name
	SAMPLE=$(echo "$SAMPLE" | perl -pe 's/.fastq.gz//')

	# new fastq full path
	DEST_FASTQ="${FASTQ_DEST_DIR}/${SAMPLE}/${DEST_FASTQ_BASE}"
	mkdir -p "${FASTQ_DEST_DIR}/${SAMPLE}"

	# print names
	echo " SAMPLE:    $SAMPLE"
	echo " ORIGINAL:  $SOURCE_FASTQ"
	echo " LINK:      $DEST_FASTQ"

	# create the symlink
	ln -s "$SOURCE_FASTQ" "$DEST_FASTQ"
	echo ""

done


#########################



# end
