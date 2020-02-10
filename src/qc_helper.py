import os

class QualityControlRNAseq:

  def trimmomaticPE(fastqr1,fastqr2,sample_id,out_dir):
    os.system('trimmomatic PE -phred33 {} {} {}/{}_L001_R1_001_paired.fq.gz {}/{}_L001_R1_001_unpaired.fq.gz {}/{}_L001_R2_001_paired.fq.gz {}/{}_L001_R2_001_unpaired.fq.gz \
               ILLUMINACLIP:/src/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35'
               .format(fastqr1,fastqr2,out_dir,sample_id,out_dir,sample_id,out_dir,sample_id,out_dir,sample_id))

  def trimmomaticSE(fastq,sample_id,out_dir):
    os.system('trimmomatic SE -phred33 {} {}/{}_L001_R1_001_paired.fq.gz {}/{}_L001_R1_001_unpaired.fq.gz \
               ILLUMINACLIP:/src/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35'
               .format(fastq,out_dir,sample_id,out_dir,sample_id))

  def run_pair():
    pass

  def run_fastqc():
    pass







