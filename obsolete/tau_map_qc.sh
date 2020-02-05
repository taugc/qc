#!/bin/bash

#data atual no formato YYYY-MM-DD
filedate=`date "+%F"`

## FIX PAR ##############################
TRIMDIR='/home/user/tools/Trimmomatic-0.38/'

## INPUT PAR ############################
#nome do PROJECT
PROJECT=$1

samples=$PROJECT'_sample.txt'

THREADS=$2

REFERENCE='/home/user/references/Sequence/BWAIndex/genome.fa'

#########################################
echo 'BWA mapping...'
## BWA mapping ##
# chk output
dir=${PROJECT}'/fastq/mapped'
if ! [ -d ${dir} ]
then
     mkdir ${dir}
fi
## run bwa
while IFS= read -r line
do
	bwa mem -t ${THREADS} -M ${REFERENCE} ${PROJECT}'/fastq/'${line}'_L001_R1_001.fastq.gz' ${PROJECT}'/fastq/'${line}'_L001_R2_001.fastq.gz' >$PROJECT'/fastq/mapped/'${line}'.aligned.sam'

done <"$samples"

echo 'samtools processing...'
# samtools sam to bam to fq

while IFS= read -r line
do
samtools view -Sb ${PROJECT}/'fastq/mapped/'${line}'.aligned.sam' >${PROJECT}/'fastq/mapped/'${line}'.aligned.bam'
samtools flagstat ${PROJECT}/'fastq/mapped/'${line}'.aligned.bam' >${PROJECT}/'fastq/mapped/'${line}'.aligned.bam.stats'
samtools view -b -F 2 ${PROJECT}/'fastq/mapped/'${line}'.aligned.bam' >${PROJECT}/'fastq/mapped/'${line}'.aligned.unmapped.bam'
samtools bam2fq -1 ${PROJECT}/'fastq/mapped/'${line}'.aligned.unmapped_R1.fq' -2 ${PROJECT}/'fastq/mapped/'${line}'.aligned.unmapped_R2.fq' ${PROJECT}/'fastq/mapped/'${line}'.aligned.unmapped.bam'

done <"$samples"

echo 'DONE '${PROJECT}
