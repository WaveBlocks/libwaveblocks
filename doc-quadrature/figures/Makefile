SRCS = $(wildcard *.gnuplot)
EPSS = $(SRCS:%.gnuplot=%.eps)
TEXS = $(SRCS:%.gnuplot=%.tex)
PDFS = $(SRCS:%.gnuplot=%-eps-converted-to.pdf)

all: $(EPSS)

clean:
	rm -f $(EPSS)
	rm -f $(TEXS)
	rm -f $(PDFS)

%.eps: %.gnuplot
	gnuplot $<
