SRCS = $(wildcard *.gnuplot)
EPSS = $(SRCS:%.gnuplot=%.eps)

all: $(EPSS)

clean:
	rm -f $(EPSS)

%.eps: %.gnuplot
	gnuplot $<
