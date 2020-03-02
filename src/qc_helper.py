# -*- coding: utf-8 -*-
import os
import subprocess

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
    total_reads = subprocess.check_output(' expr $(zcat {} | wc -l) / 4'.format(path), shell=True)
    return int(total_reads.strip())

  def get_text_quality_base(self):
    text = "No eixo X está a posição de cada base na 'read'. O eixo Y mostra o 'score' de qualidade. Quanto maior o 'score', melhor. O plano de fundo do gráfico divide o eixo Y em 'scores' bons (verde), razoáveis (laranja) e de baixa qualidade (vermelho)."
    return text

  def get_text_quality_seq(self):
    text = "No eixo X temos a média de qualidade ao logo das sequências de todas as leituras (reads) e no eixo Y o número total de 'reads'. É importante notar que a qualidade média das reads deve apresentar um pico na faixa superior do gráfico."
    return text

  def get_text_seq_content(self):
    text = "Porcentagem de cada um dos quatro nucleotídeos (T, C, A, G) em cada posição em todas as reads. O esperado é que haja pouca ou nenhuma diferença entre as quatro bases. A proporção de cada uma deve ser relativamente constante ao longo do comprimento da read com %A = %T e %G = %C, sendo que as linhas neste gráfico devem ser paralelas entre si."
    return text

  def get_text_duplication(self):
    text = "Espera-se que em uma biblioteca diversificada, a maioria das sequências ocorra apenas uma vez no conjunto final. Um baixo nível de duplicação pode indicar um nível alto de cobertura da sequência alvo, mas um alto nível de duplicação provavelmente indica algum viés de enriquecimento. A linha azul mostra a distribuição dos níveis de duplicação para o conjunto completo de sequências, já a linha vermelha apresenta a distribuição para as sequências desduplicadas. A maioria das sequências deve cair na extremidade esquerda do gráfico nas linhas vermelha e azul. Isso indica uma biblioteca altamente diversificada que não foi sequenciada em excesso."
    return text

  def get_text_adapters(self):
    text = "No eixo X está a posição nas sequências em pares de base e no eixo Y a porcentagem de adaptador."
    return text









