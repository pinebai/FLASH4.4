TARGET=flash_chombo_applications

latex_source=$(TARGET).tex

default: pdf

pdf: $(TARGET).pdf

$(TARGET).pdf: $(latex_source)
	pdflatex $(TARGET).tex
	pdflatex $(TARGET).tex
	pdflatex $(TARGET).tex

clean:
	rm -f ${TARGET}.{log,aux,bbl,blg,log,nav,toc,out,snm,vrb}

realclean: clean
	rm -f $(TARGET).pdf
