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
                    choices=['S', 'P'],
                    help="Choose the type of sequencing: (S): Single-end (P): Pair-end.")
parser.add_argument("s",
                    nargs=1,
                    type=int,
                    action="store",
                    choices=[0, 1],
                    help="Choose one software option to trimming: (0): Trimmomatic (1): Trim_galore.")
parser.add_argument("f",
                    nargs=1,
                    type=argparse.FileType('r'),
                    help="Forward fastq.gz file.")
parser.add_argument("-r",
                    nargs=1,
                    type=argparse.FileType('r'),
                    help="Reverse fastq.gz file.")


args = parser.parse_args()

qc = qc_helper.QualityControl()
out_dir = '{}_out_dir'.format(args.id)

if not os.path.exists(out_dir):
    out_dir = '{}_out_dir'.format(args.id)
    os.makedirs(out_dir)








