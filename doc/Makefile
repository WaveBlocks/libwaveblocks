%.svg.pdf~: %.svg
	touch $*.svg.pdf~
	inkscape -z -T -A $*.pdf $*.svg

svgfiles: $(patsubst %.svg,%.svg.pdf~,$(wildcard *.svg))

manual.pdf: manual.tex svgfiles
	pdflatex manual.tex
#	pdflatex --shell-escape manual.tex	

pdf: manual.pdf

.PHONY clean:
	rm -f *.pdf
	rm -f *.log
	rm -f *.pdf_tex
	rm -f *.aux
	rm -f *.pdf~
