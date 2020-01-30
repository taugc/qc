#!/bin/bash

#data atual no formato YYYY-MM-DD
filedate=`date "+%F"`


## FIX PAR ##############################
TOOLDIR='/home/user/tools/sortmerna-3.0.3-Linux'

## INPUT PAR ############################

#numero de threads

THREADS=$1


#########################################

## Merge-paired reads, do sortmerna and unmerge ##

for file in $(ls *_R1_001_paired.trim.fastq*);do gunzip ${file%_R1_paired.trim.fastq*}'_R1_paired.trim.fastq.gz';gunzip ${file%_R1_paired.trim.fastq.gz}'_R2_paired.trim.fastq.gz';${TOOLDIR}/scripts/merge-paired-reads.sh ${file%_R1_paired.trim.fastq*}'_R1_paired.trim.fastq' ${file%_R1_paired.trim.fastq*}'_R2_paired.trim.fastq' ${file%_R1_paired.trim.fastq*}'_paired.trim.fastq';time ${TOOLDIR}/bin/sortmerna --ref ${TOOLDIR}'/rRNA_databases/silva-euk-18s-id95.fasta,'${TOOLDIR}'/index/silva-euk-18s-db:'${TOOLDIR}'/rRNA_databases/silva-euk-28s-id98.fasta,'${TOOLDIR}'/index/silva-euk-28s-db --reads ${file%_R1_paired.trim.fastq*}'_paired.trim.fastq' --sam --num_alignments 1 --fastx --paired_in --log -a ${THREADS} -v; rm ${file%_R1_paired.trim.fastq*}'_paired.trim_rRNA.fastq'; rm *.sam;gzip ${file%_R1_paired.trim.fastq*}'_R2_paired.trim.fastq';gzip ${file%_R1_paired.trim.fastq*}'_R1_paired.trim.fastq'; done;

##

echo 'DONE'
