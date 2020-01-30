#!/bin/bash

#data atual no formato YYYY-MM-DD
filedate=`date "+%F"`


## FIX PAR ##############################
TRIMDIR='/home/user/tools/Trimmomatic-0.38/'

## INPUT PAR ############################
#nome do PROJECT
PROJECT=$1

samples=$PROJECT'_sample.txt' 

DIR_SEQ=$2

#########################################

## Pre Trim FastQC ##
# chk output
dir=${PROJECT}/'fastq/fastqc'
if ! [ -d ${dir} ]
then
     mkdir ${dir}
fi
# run fastqc
for file in $(ls ${PROJECT}/'fastq/'*.fastq.gz); do /home/jonas/tools/FastQC/fastqc $file -o ${PROJECT}/'fastq/fastqc' ; done;

## Trimmomatic ##
# chk output 
dir=${PROJECT}/'fastq/trim'
if ! [ =d ${dir} ]
then
	mkdir ${dir}
fi
# run trimmomatic
while IFS= read -r line
do
	java -jar ${TRIMDIR}trimmomatic-0.38.jar PE -phred33 ${PROJECT}/'fastq/'${line}'_L001_R1_001.fastq.gz' ${PROJECT}/'fastq/'${line}'_L001_R2_001.fastq.gz' ${PROJECT}/'fastq/trim/'${line}'_L001_R1_001_paired.fq.gz' ${PROJECT}/'fastq/trim/'${line}'_L001_R1_001_unpaired.fq.gz' ${PROJECT}/'fastq/trim/'${line}'_L001_R2_001_paired.fq.gz' ${PROJECT}/'fastq/trim/'${line}'_L001_R2_001_unpaired.fq.gz' ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35
done <"$samples"

## Pos Trim FastQC
 chk output
dir=${PROJECT}/'fastq/trim/fastqc'
if ! [ =d ${dir} ]
then
	mkdir ${dir}
fi
#run fastqc
for file in $(ls ${PROJECT}/'fastq/trim/'*_paired*); do /home/jonas/tools/FastQC/fastqc $file -o ${PROJECT}/'fastq/trim/fastqc'; done;

echo 'DONE '${PROJECT}

