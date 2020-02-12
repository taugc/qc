# -*- coding: utf-8 -*-

import argparse
import os
import sys
from src import qc_helper

parser = argparse.ArgumentParser(description='QC TauGC: used to make the quality control from NGS fastq files and generate an automatic report.')

parser.add_argument("id",
                    nargs=1,
                    action="store",
                    type=str,
                    help="Sample name")
parser.add_argument("t",
                    nargs=1,
                    action="store",
                    type=str,
                    choices=['SE', 'PE'],
                    help="Choose the type of sequencing: (SE): Single-end (PE): Pair-end.")
parser.add_argument("e",
                    nargs=1,
                    type=int,
                    action="store",
                    choices=[0, 1],
                    help="Choose experiment: (0): RNAseq (1): Others")
parser.add_argument("f",
                    nargs=1,
                    type=str,
                    help="Forward fastq.gz file.")
parser.add_argument("-r",
                    nargs=1,
                    type=str,
                    help="Reverse fastq.gz file.")

args = parser.parse_args()

qc_rnaseq = qc_helper.QualityControl()
experiment = args.e[0]
type_seq = args.t[0]
sample = args.id[0]
out_dir = '{}_out_dir'.format(sample)
adapters = os.path.abspath('src/adapters.fa')

print(adapters)

if not os.path.exists(out_dir):
    out_dir = '{}_out_dir'.format(sample)
    os.makedirs(out_dir)

print(experiment,type_seq,sample,out_dir)

############################ RNAseq ######################################
if experiment == 0:
  if type_seq == 'PE' or type_seq == 'pe':
    r1=args.f[0]
    r2=args.r[0]
    os.system('trimmomatic PE -phred33 {} {} {}/{}_L001_R1_001_paired.fq.gz {}/{}_L001_R1_001_unpaired.fq.gz {}/{}_L001_R2_001_paired.fq.gz {}/{}_L001_R2_001_unpaired.fq.gz \
               ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35'
               .format(r1,r2,out_dir,sample,out_dir,sample,out_dir,sample,out_dir,sample,adapters))
    
    os.system('fastqc {}/*_paired*' .format(out_dir))
   
  elif type_seq == 'SE' or type_seq == 'se':
    r1=args.f[0]  
    os.system('trimmomatic SE -phred33 {} {}/{}_L001_R1_001_paired.fq.gz \
               ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35'
               .format(r1,out_dir,sample,adapters))

    os.system('fastqc {}/*_paired*' .format(out_dir))
    
  else:
    print('Error: Invalid parameter!')

else:
  print('Sorry, there is just RNAseq module for while!')











