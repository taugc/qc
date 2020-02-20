# -*- coding: utf-8 -*-

import argparse
import os
import sys
from datetime import date
from src import qc_helper
from reportlab.pdfgen.canvas import Canvas
from reportlab.graphics.shapes import *
from reportlab.graphics import renderPDF
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.platypus import Image, Paragraph, Table, BaseDocTemplate, TableStyle, PageTemplate, SimpleDocTemplate, Spacer, PageBreak, Flowable, NextPageTemplate
from reportlab.lib.pagesizes import A4,landscape, letter
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.units import inch, cm, mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT, TA_CENTER, TA_RIGHT

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

pr = qc_helper.ProcessResult()
experiment = args.e[0]
type_seq = args.t[0]
sample = args.id[0]
r1 = args.f[0]
out_dir = '{}_out_dir'.format(sample)
adapters = os.path.abspath('src/adapters.fa')

if not os.path.exists(out_dir):
    out_dir = '{}_out_dir'.format(sample)
    os.makedirs(out_dir)

n_reads = pr.get_n_reads(r1)
sequencing = pr.get_typeseq(type_seq)

##########################################################################
############################ RNAseq ######################################
##########################################################################
if experiment == 0:
  if type_seq == 'PE' or type_seq == 'pe':
    r1=args.f[0]
    r2=args.r[0]
    pr.run_trimPE(r1,r2,out_dir,sample,adapters)
    pr.run_fastQC(out_dir)
    os.system('unzip {}/{}_L001_R1_001_paired_fastqc.zip'.format(out_dir, sample))
    os.system('unzip {}/{}_L001_R2_001_paired_fastqc.zip'.format(out_dir, sample))

  elif type_seq == 'SE' or type_seq == 'se':
    r1=args.f[0]
    pr.run_trimSE(r1,out_dir,sample,adapters)
    pr.run_fastQC(out_dir)
    os.system('unzip {}/{}_L001_R1_001_paired_fastqc.zip'.format(out_dir, sample))

  else:
    print('Error: Invalid parameter!')

else:
  print('Sorry, there is just RNAseq module for while!')

fastq_trimmed = os.path.join('{}/{}_L001_R1_001_paired.fq.gz' .format(out_dir, sample))
n_reads_after_trim = pr.get_n_reads(fastq_trimmed)
percent = (n_reads_after_trim / n_reads) * 100
n_reads_trimmed_perc = '{} ({}%)'.format(n_reads_after_trim, round(percent, 2))


##########################################################################
############################ REPORT ######################################
##########################################################################
doc = SimpleDocTemplate('QCReport_{}.pdf'.format(sample))
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))
styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
styles.add(ParagraphStyle(name='Left', alignment=TA_LEFT))
styles.add(ParagraphStyle(name='LeftBold', alignment=TA_LEFT, fontName='Helvetica-Bold'))
styles.add(ParagraphStyle(name='Right', alignment=TA_RIGHT))
styles.add(ParagraphStyle(name='Infos', alignment=TA_RIGHT, fontName='Helvetica-Bold'))
styles.add(ParagraphStyle(name='Disclaimer', alignment=TA_CENTER, fontName='Helvetica-Bold'))
styles.add(ParagraphStyle(name='DisclaimerText', alignment=TA_JUSTIFY, fontName='Helvetica', firstLineIndent=40))

spacer = Spacer(0, 0.7*inch)
Story = []
######################################
#################### First page
######################################
# logo
logo = os.path.abspath('images/taugc.png')
im = Image(logo, 3*inch,3*inch)
Story.append(im)
Story.append(spacer)

# rect blue
d = Drawing(500, 200)
d.add(Rect(-75, 120, 280, 40, fillColor=colors.darkblue))
d.add(String(-68,132, 'Quality Control Report', fontSize=25, fontName='Helvetica-Bold', fillColor=colors.white))
Story.append(d)

# rect blue
today = date.today()
x = today.strftime('%B, %d %Y')
e = Drawing(500, 200)
e.add(Rect(233, 320, 280, 40, fillColor=colors.darkgreen))
e.add(String(260,330, x, fontSize=25, fontName='Helvetica-Bold', fillColor=colors.white))
Story.append(e)

# Contact
Story.append(Spacer(0, -3*inch))
email = '<font size=12>%s</font>' % "Contact us: info@taugc.com"
Story.append(Paragraph(email, styles["Infos"]))
Story.append(Spacer(0, 0.2*inch))
site = '<font size=12>%s</font>' % "About us: www.taugc.com"
Story.append(Paragraph(site, styles["Infos"]))
Story.append(Spacer(0, 1*inch))

#Disclaimer first page
disclaimer_title = '<font size=13 textColor="black">%s</font>' % "DISCLAIMER"
Story.append(Paragraph(disclaimer_title, styles["Disclaimer"]))
Story.append(Spacer(0, 0.2*inch))
Story.append(Spacer(5, 0))
disclaimer_text = '<font size=10 textColor="grey">%s</font>' % "  This report is intended for general guidance and information purposes only. This report is under no circumstances intended to be used or considered directly for any kind of scientific or academic publication and clinical diagnosis.  Please note that this is reference material for the original authors of the study. The report is to be considered only as a tool for hypotheses and data interpretation by the authors along with all data provided."
Story.append(Paragraph(disclaimer_text, styles["DisclaimerText"]))
Story.append(PageBreak())
######################################
#################### Second page
######################################
logo = os.path.abspath('images/taugc.png')
im = Image(logo, 0.8*inch,0.8*inch)
im.hAlign = 'LEFT'
im.vAlign = 'TOP'
Story.append(im)

disclaimer_title = '<font size=13 textColor="black">%s</font>' % "DISCLAIMER"
Story.append(Paragraph(disclaimer_title, styles["Disclaimer"]))
Story.append(Spacer(0, 0.2*inch))
disclaimer_2 = '<font size=10 textColor="grey">%s</font>' % "The data used for constructing this material in the report is obtained directly from the customer (customer). We have taken reasonable care to ensure that, and to the best of our knowledge, material information contained herein is in accordance with information and guidelines provided by the customer. All data presented are prepared by TauGC Bioinformatics in accordance with best practices from third parties, as part of the service included in the report. The estimates are subject to uncertainties, misinterpretations and other factors that may cause actual events to differ materially from any anticipated hypotheses. Please note that we make no assurance that the underlying statements are free from errors. Neither the authors nor TauGC Bioinformatics is making any representation or warranty, express or implied, as to the accuracy or completeness of data presented in this report and none of the authors or TauGC Bioinformatics will have any liability towards any other person resulting in incorrect use this report. TauGC do not intend, and do not assume any obligation to update or correct the information included in this report."
Story.append(Paragraph(disclaimer_2, styles["DisclaimerText"]))
Story.append(Spacer(0, 0.1*inch))
disclaimer_3 = '<font size=10 textColor="grey">%s</font>' % "The information contained herein may be subject to interpretation by authors. TauGC does not accept any form of liability, neither legally nor financially, for loss (direct or indirect) caused by the understanding and/or use of this report or its content."
Story.append(Paragraph(disclaimer_3, styles["DisclaimerText"]))
Story.append(Spacer(0, 0.1*inch))
disclaimer_4 = '<font size=10 textColor="grey">%s</font>' % "This report is only intended for the recipients, and should not be copied or otherwise distributed, in whole or in part, to any other person. This report is subject to Brazilian law, and any dispute arising in respect of this report is subject to the exclusive jurisdiction of Brazilian courts."
Story.append(Paragraph(disclaimer_4, styles["DisclaimerText"]))

Story.append(Spacer(0, 0.3*inch))
disclaimer_5 = '<font size=12 textColor="grey">%s</font>' % "CONSIDERATIONS ABOUT AUTHORSHIP & LIABILITY"
Story.append(Paragraph(disclaimer_5, styles["Disclaimer"]))
Story.append(Spacer(0, 0.2*inch))
disclaimer_6 = '<font size=10 textColor="grey">%s</font>' % "TauGC Bioinformatics is not to be considered as author or co-author of any kind publication solely based on direct use of this report. Any consideration for authorship by TauGC must be discussed further by authors and TauGC Bioinformatics."
Story.append(Paragraph(disclaimer_6, styles["DisclaimerText"]))
Story.append(Spacer(0, 0.1*inch))
disclaimer_7 = '<font size=10 textColor="grey">%s</font>' % "TauGC and its members follow the Nature Journals Policy for Authorship (https://www.nature.com/authors/policies/authorship.html) adapted from McNutt et al. 2018 (DOI: 10.1073/pnas.1715374115)."
Story.append(Paragraph(disclaimer_7, styles["DisclaimerText"]))

# footer
Story.append(Spacer(0, 2.7*inch))
footer_1 = '<font size=8 textColor="grey">%s</font>' % "TauGC Bioinformatics - Rua Apiacas 886, São Paulo - SP"
Story.append(Paragraph(footer_1, styles["Center"]))
Story.append(Spacer(0, 0.1*inch))
footer_2 = '<font size=8 textColor="grey">%s</font>' % "This report is only intended for the recipients, and should not be copied or otherwise distributed, in whole or in part, to any other person."
Story.append(Paragraph(footer_2, styles["Center"]))
Story.append(PageBreak())

######################################
#################### Third page
######################################
logo = os.path.abspath('images/taugc.png')
im = Image(logo, 0.8*inch,0.8*inch)
im.hAlign = 'LEFT'
im.vAlign = 'TOP'
Story.append(im)

title = '<font size=14 >%s</font>' % "Controle de Qualidade"
Story.append(Paragraph(title, styles["Disclaimer"]))
Story.append(Spacer(0, 0.7*inch))

text = '<font size=12 >%s</font>' % "Tab. 1- Informações gerais sobre os dados."
Story.append(Paragraph(text, styles["Left"]))
Story.append(Spacer(0, 0.1*inch))

header = ['Amostra', 'Sequenciamento', 'No. reads', 'No. reads após trimming']
data = [header, [sample, sequencing, n_reads, n_reads_trimmed_perc]]

t = Table(data)
t.setStyle(TableStyle([
      ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
      ('BOX', (0,0), (-1,-1), 0.25, colors.black),
      ('ALIGN',(0,-1),(-1,-1),'CENTER'),
      ('VALIGN',(0,-1),(-1,-1),'MIDDLE'),
      ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
      ('TEXTCOLOR',(0,0),(-1,0),colors.white),
      ('FONTSIZE',(0,0),(-1,-1), 12),
      ('TEXTFONT',(0,0),(-1,0), 'Helvetica-Bold')
  ]))
Story.append(t)

chart_style = TableStyle([('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                          ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                          ('BOX', (0,0), (-1,-1), 0.25, colors.black)])
if type_seq == 'PE' or type_seq == 'pe':
  ###### Sequence lenght distribution
  Story.append(Spacer(0, 1.5*inch))
  text = '<font size=12 >%s</font>' % "Fig. 1- Distribuição do tamanho das reads (R1 e R2)."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/sequence_length_distribution.png' .format(sample))
  qper_base2 = os.path.join('{}_L001_R2_001_paired_fastqc/Images/sequence_length_distribution.png' .format(sample))
  I1 = Image(qper_base1)
  I2 = Image(qper_base2)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  I2.drawHeight = 3.0*inch*I2.drawHeight / I2.drawWidth
  I2.drawWidth = 3.0*inch
  data = [[I1, I2]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))



  ###### quality per base
  Story.append(Spacer(0, 2.0*inch))
  logo = os.path.abspath('images/taugc.png')
  im = Image(logo, 0.8*inch,0.8*inch)
  im.hAlign = 'LEFT'
  im.vAlign = 'TOP'
  Story.append(im)

  Story.append(Spacer(0, 0.3*inch))
  text = '<font size=12 >%s</font>' % "Fig. 2- Qualidade por base das sequências (R1 e R2)."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/per_base_quality.png' .format(sample))
  qper_base2 = os.path.join('{}_L001_R2_001_paired_fastqc/Images/per_base_quality.png' .format(sample))
  I1 = Image(qper_base1)
  I2 = Image(qper_base2)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  I2.drawHeight = 3.0*inch*I2.drawHeight / I2.drawWidth
  I2.drawWidth = 3.0*inch
  data = [[I1, I2]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_quality_base()
  Story.append(Paragraph(text, styles["Justify"]))

  ###### quality per seq
  Story.append(Spacer(0, 0.5*inch))
  text = '<font size=12 >%s</font>' % "Fig. 3- Média de qualidade por sequência (R1 e R2)."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/per_sequence_quality.png' .format(sample))
  qper_base2 = os.path.join('{}_L001_R2_001_paired_fastqc/Images/per_sequence_quality.png' .format(sample))
  I1 = Image(qper_base1)
  I2 = Image(qper_base2)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  I2.drawHeight = 3.0*inch*I2.drawHeight / I2.drawWidth
  I2.drawWidth = 3.0*inch
  data = [[I1, I2]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_quality_seq()
  Story.append(Paragraph(text, styles["Justify"]))

###### per seq content
  Story.append(Spacer(0, 0.5*inch))
  logo = os.path.abspath('images/taugc.png')
  im = Image(logo, 0.8*inch,0.8*inch)
  im.hAlign = 'LEFT'
  im.vAlign = 'TOP'
  Story.append(im)
  Story.append(Spacer(0, 0.3*inch))

  text = '<font size=12 >%s</font>' % "Fig. 3- Conteúdo das reads por base (R1 e R2)."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/per_base_sequence_content.png' .format(sample))
  qper_base2 = os.path.join('{}_L001_R2_001_paired_fastqc/Images/per_base_sequence_content.png' .format(sample))
  I1 = Image(qper_base1)
  I2 = Image(qper_base2)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  I2.drawHeight = 3.0*inch*I2.drawHeight / I2.drawWidth
  I2.drawWidth = 3.0*inch
  data = [[I1, I2]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_seq_content()
  Story.append(Paragraph(text, styles["Justify"]))

  ############ Level duplication
  Story.append(Spacer(0, 0.4*inch))
  text = '<font size=12 >%s</font>' % "Fig. 4- Nível de duplicação das sequências (R1 e R2)."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/duplication_levels.png' .format(sample))
  qper_base2 = os.path.join('{}_L001_R2_001_paired_fastqc/Images/duplication_levels.png' .format(sample))
  I1 = Image(qper_base1)
  I2 = Image(qper_base2)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  I2.drawHeight = 3.0*inch*I2.drawHeight / I2.drawWidth
  I2.drawWidth = 3.0*inch
  data = [[I1, I2]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_duplication()
  Story.append(Paragraph(text, styles["Justify"]))

  ###### Adapters
  Story.append(Spacer(0, 0.5*inch))
  logo = os.path.abspath('images/taugc.png')
  im = Image(logo, 0.8*inch,0.8*inch)
  im.hAlign = 'LEFT'
  im.vAlign = 'TOP'
  Story.append(im)
  Story.append(Spacer(0, 0.3*inch))

  text = '<font size=12 >%s</font>' % "Fig. 5- Adaptadores nas sequências analisadas (R1 e R2)."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/adapter_content.png' .format(sample))
  qper_base2 = os.path.join('{}_L001_R2_001_paired_fastqc/Images/adapter_content.png' .format(sample))
  I1 = Image(qper_base1)
  I2 = Image(qper_base2)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  I2.drawHeight = 3.0*inch*I2.drawHeight / I2.drawWidth
  I2.drawWidth = 3.0*inch
  data = [[I1, I2]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_adapters()
  Story.append(Paragraph(text, styles["Justify"]))

elif type_seq == 'SE' or type_seq == 'se':
  ###### Sequence lenght distribution
  Story.append(Spacer(0, 1.5*inch))
  text = '<font size=10 >%s</font>' % "Fig. 1- Distribuição do tamanho das reads."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/sequence_length_distribution.png' .format(sample))
  I1 = Image(qper_base1)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  data = [[I1]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))


  ###### quality per base
  Story.append(Spacer(0, 2.0*inch))
  logo = os.path.abspath('images/taugc.png')
  im = Image(logo, 0.8*inch,0.8*inch)
  im.hAlign = 'LEFT'
  im.vAlign = 'TOP'
  Story.append(im)

  Story.append(Spacer(0, 0.3*inch))
  text = '<font size=12 >%s</font>' % "Fig. 2- Qualidade por base das sequências."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/per_base_quality.png' .format(sample))
  I1 = Image(qper_base1)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  data = [[I1]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_quality_base()
  Story.append(Paragraph(text, styles["Justify"]))

  ###### quality per seq
  Story.append(Spacer(0, 0.5*inch))
  text = '<font size=10 >%s</font>' % "Fig. 3- Média de qualidade por sequência."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/per_sequence_quality.png' .format(sample))
  I1 = Image(qper_base1)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  data = [[I1]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_quality_seq()
  Story.append(Paragraph(text, styles["Justify"]))

###### per seq content
  Story.append(Spacer(0, 0.5*inch))
  logo = os.path.abspath('images/taugc.png')
  im = Image(logo, 0.8*inch,0.8*inch)
  im.hAlign = 'LEFT'
  im.vAlign = 'TOP'
  Story.append(im)
  Story.append(Spacer(0, 0.3*inch))

  text = '<font size=12 >%s</font>' % "Fig. 3- Conteúdo das reads por base."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/per_base_sequence_content.png' .format(sample))
  I1 = Image(qper_base1)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  data = [[I1]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_seq_content()
  Story.append(Paragraph(text, styles["Justify"]))

  ############ Level duplication
  Story.append(Spacer(0, 0.4*inch))
  text = '<font size=12 >%s</font>' % "Fig. 4- Nível de duplicação das sequências."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/duplication_levels.png' .format(sample))
  I1 = Image(qper_base1)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  data = [[I1]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_duplication()
  Story.append(Paragraph(text, styles["Justify"]))

  ###### Adapters
  Story.append(Spacer(0, 0.5*inch))
  logo = os.path.abspath('images/taugc.png')
  im = Image(logo, 0.8*inch,0.8*inch)
  im.hAlign = 'LEFT'
  im.vAlign = 'TOP'
  Story.append(im)
  Story.append(Spacer(0, 0.3*inch))

  text = '<font size=12 >%s</font>' % "Fig. 5- Adaptadores nas sequências analisadas."
  Story.append(Paragraph(text, styles["Left"]))
  Story.append(Spacer(0, 0.1*inch))
  qper_base1 = os.path.join('{}_L001_R1_001_paired_fastqc/Images/adapter_content.png' .format(sample))
  I1 = Image(qper_base1)
  I1.drawHeight = 3.0*inch*I1.drawHeight / I1.drawWidth
  I1.drawWidth = 3.0*inch
  data = [[I1]]
  t= Table(data)
  t.setStyle(chart_style)
  Story.append(t)
  Story.append(Spacer(0, 0.2*inch))
  text = '<font size=12 >%s</font>' % pr.get_text_adapters()
  Story.append(Paragraph(text, styles["Justify"]))

# References
Story.append(Spacer(0, 0.9*inch))
title = '<font size=12 >%s</font>' % "References"
Story.append(Paragraph(title, styles["LeftBold"]))
Story.append(Spacer(0, 0.2*inch))
ref1 = '<font size=10 >%s</font>' % "Batut et al., 2018. Community-Driven Data Analysis Training for Biology. Cell Systems, 6(6):752-758."
Story.append(Paragraph(ref1, styles["Left"]))
Story.append(Spacer(0, 0.2*inch))
ref2 = '<font size=10 >%s</font>' % "Bérénice Batut, 2019. Quality Control. /training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html. Online; accessed 2020-02-20."
Story.append(Paragraph(ref2, styles["Left"]))

doc.build(Story)

os.system(" rm -fr {}_*_fastqc" .format(sample))

print("################################################\n")
print("#### The Quality Control Report was created! ###\n")
print("################################################\n")
