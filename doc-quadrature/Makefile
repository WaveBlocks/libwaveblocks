PDFLATEX = pdflatex
BIBTEX = bibtex
SECTIONS = $(wildcard sections/*.tex)
NAME = manual-quadrature

all: $(NAME).pdf

$(NAME).pdf: $(NAME).tex $(SECTIONS) figs
	$(PDFLATEX) $(NAME).tex
	$(BIBTEX) $(NAME).aux
	$(PDFLATEX) $(NAME).tex
	$(PDFLATEX) $(NAME).tex

figs:
	$(MAKE) -C figures

clean:
	$(MAKE) clean -C figures
	rm -f $(NAME).aux
	rm -f $(NAME).bbl
	rm -f $(NAME).blg
	rm -f $(NAME).log
	rm -f $(NAME).out
	rm -f $(NAME).toc
