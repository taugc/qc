# -*- coding: utf-8 -*-
from datetime import date
import os
from reportlab.pdfgen import canvas

def firstPage(doc):
  logo = os.path.abspath('images/taugc.png')
  doc.drawImage(logo, 200, 450)

  doc.setStrokeColorRGB(0,0,0.5)
  doc.setFillColorRGB(0,0,0.5)
  doc.rect(0,300, 290,40, fill=1)

  doc.setStrokeColorRGB(0,0.5,0)
  doc.setFillColorRGB(0,0.5,0)
  doc.rect(310,300, 285 ,40, fill=1)

  doc.setFillColorRGB(255,255,255)
  doc.setFont('Helvetica-Bold', 25)
  doc.drawString(12,310, "Quality Control Report")
  
  today = date.today()
  x = today.strftime('%B, %d %Y')
  doc.setFillColorRGB(255,255,255)
  doc.setFont('Helvetica-Bold', 25)
  doc.drawString(350,310, x)

  doc.setFillColorRGB(0,0,0)
  doc.setFont('Helvetica-Bold', 20)
  doc.drawString(265,200, "Contact:")

  doc.setFillColorRGB(0,0,0)
  doc.setFont('Helvetica', 20)
  doc.drawString(350,200, "info@taugc.com")
 
  



sample = 'TESTE'
doc = canvas.Canvas('QCReport_{}.pdf'.format(sample))   

firstPage(doc)
doc.showPage()
doc.save()



