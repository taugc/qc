#!/bin/bash

#data atual no formato YYYY-MM-DD
filedate=`date "+%F"`


## FIX PAR ##############################
SUBSETSCRIPT='/home/user/tools'

## INPUT PAR ############################
#nome do PROJECT
PROJECT=$1

samples=$PROJECT'_sample.txt'

FRACTION=$2
#########################################
echo 'Check par...'
# chk output
dir=${PROJECT}/'fastq/trim/subset'
if ! [ -d ${dir} ]
then
     mkdir ${dir}
fi

echo 'subsetting reads...'
# spades

while IFS= read -r line
do
	'python' ${SUBSETSCRIPT}'/subsetfastq.py' ${FRACTION} ${PROJECT}'/fastq/trim/'${line}'_L001_R1_001_paired.fq.gz' ${PROJECT}'/fastq/trim/'${line}'_L001_R2_001_paired.fq.gz' ${PROJECT}/'fastq/trim/subset/'${line}'_L001_R1_001_paired.fq.gz' ${PROJECT}'/fastq/trim/subset/'${line}'_L001_R2_001_paired.fq.gz'
done <"$samples"

echo 'DONE '${PROJECT}
