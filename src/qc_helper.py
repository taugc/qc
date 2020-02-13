# -*- coding: utf-8 -*-
from datetime import date
import os
from reportlab.pdfgen import canvas
from textwrap import wrap


class QCReport:

  def firstPage(canvas,doc):
    logo = os.path.abspath('images/taugc.png')
    doc.drawImage(logo, 200, 500)

    doc.setStrokeColorRGB(0,0,0.5)
    doc.setFillColorRGB(0,0,0.5)
    doc.rect(0,400, 290,40, fill=1)

    doc.setStrokeColorRGB(0,0.5,0)
    doc.setFillColorRGB(0,0.5,0)
    doc.rect(310,400, 285 ,40, fill=1)

    doc.setFillColorRGB(255,255,255)
    doc.setFont('Helvetica-Bold', 25)
    doc.drawString(12,410, "Quality Control Report")

    today = date.today()
    x = today.strftime('%B, %d %Y')
    doc.setFillColorRGB(255,255,255)
    doc.setFont('Helvetica-Bold', 25)
    doc.drawString(350,410, x)

    doc.setFillColorRGB(0,0,0)
    doc.setFont('Helvetica-Bold', 20)
    doc.drawString(265,270, "Contact:")

    doc.setFillColorRGB(0,0,0)
    doc.setFont('Helvetica', 20)
    doc.drawString(350,270, "info@taugc.com")

    doc.setFillColorRGB(0.5,0.5,0.5)
    doc.setFont('Helvetica-Bold', 14)
    doc.drawString(255,180, "Disclaimer")

    t = doc.beginText()
    t.setFont('Helvetica', 10)
    t.setFillColorRGB(0.5,0.5,0.5)
    t.setCharSpace(2)
    t.setTextOrigin(30, 150)
    text = "This report is intended for general guidance and information purposes only. This report is under no circumstances intended to be used or considered directly for any kind of scientific or academic publication and clinical diagnosis. Please note that this is reference material for the original authors of the study. The report is to be considered only as a tool for hypotheses and data interpretation by the authors along with all data provided."
    wraped_text = "\n".join(wrap(text, 85)) # 85 is line width
    t.textLines(wraped_text)
    doc.drawText(t)












