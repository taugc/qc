# -*- coding: utf-8 -*-
import os

class ProcessResult:

  def run_trimPE(self,r1,r2,out_dir,sample,adapters):
    os.system('trimmomatic PE -phred33 {} {} {}/{}_L001_R1_001_paired.fq.gz {}/{}_L001_R1_001_unpaired.fq.gz {}/{}_L001_R2_001_paired.fq.gz {}/{}_L001_R2_001_unpaired.fq.gz \
              ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35'
              .format(r1,r2,out_dir,sample,out_dir,sample,out_dir,sample,out_dir,sample,adapters))

  def run_trimSE(self, r1,out_dir,sample,adapters):
    os.system('trimmomatic SE -phred33 {} {}/{}_L001_R1_001_paired.fq.gz \
              ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35'
              .format(r1,out_dir,sample,adapters))

  def run_fastQC(self, out_dir):
    os.system('fastqc {}/*_paired*' .format(out_dir))

  def get_typeseq(self, typeseq):
    if typeseq == 'PE' or typeseq == 'pe':
      return 'Pair-End'
    elif typeseq == 'SE' or typeseq == 'se':
      return 'Single-End'
    else:
      return 'Unknow'

  def get_n_reads(self, path):
    lines = os.system('zcat {} | wc -l'.format(path))
    total_reads = int(lines) / 4
    return total_reads

  









