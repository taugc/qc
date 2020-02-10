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
                    type=argparse.FileType('r'),
                    help="Forward fastq.gz file.")
parser.add_argument("-r",
                    nargs=1,
                    type=argparse.FileType('r'),
                    help="Reverse fastq.gz file.")

args = parser.parse_args()

qc_rnaseq = qc_helper.QualityControlRNAseq()
out_dir = '{}_out_dir'.format(args.id)

if not os.path.exists(out_dir):
    out_dir = '{}_out_dir'.format(args.id)
    os.makedirs(out_dir)

# RNAseq
if args.e == 0:

  if args.t == 'PE' or args.t == 'pe':
    qc_rnaseq.trimmomaticPE(args.f,args.r,args.id,out_dir)
    #run fastqc

  elif args.t == 'SE' or args.t == 'se':
    qc_rnaseq.trimmomaticSE(args.f,args.id,out_dir)
    #run fastqc

  #Create report

else:
  print('Sorry, there is just RNAseq module for while!')











