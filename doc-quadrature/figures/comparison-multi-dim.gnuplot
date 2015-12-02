set term epslatex color
set key top left
set xlabel "number of dimensions"
set ylabel "runtime [ms]"
set log y
set grid
set xtics 1
set output "comparison-multi-dim.tex"
plot "-" with lp pt 7 title "Tensor-Product order 3", \
     "-" with lp pt 2 lc 3 title "Genz-Keister level 3", \
     "-" with lp pt 3 lc 4 title "Genz-Keister level 2"

  1 0.00831543
  2 0.0187738
  3 0.058378
  4 0.396888
  5 5.12689
  6 96.8282
  7 2711.88
  8 74816.3
e
  1 0.00803576
  2 0.018205
  3 0.04863
  4 0.206399
  5 1.11746
  6 11.8937
  7 145.012
  8 1889.25
  9 19316.7
e
  1 0.00853511
  2 0.0178002
  3 0.0346079
  4 0.109373
  5 0.480391
  6 8.54215
  7 71.4104
  8 605.528
  9 6871.61
e

set output
