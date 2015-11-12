set term epslatex color
set key top left
set xlabel "number of components"
set ylabel "runtime [ms]" rotate by 0
set log x 2
set log y
set grid
set output "speedup-multi-comp.tex"
plot "-" with lp pt 7 title "Python", "-" with lp pt 2 lc 3 title "C++"
  1 42.4715995789
  2 170.664310455
  3 382.325196266
  4 681.163907051
  5 1066.02609158
  6 1533.49120617
  7 2087.76860237
  8 2726.53651237
e
  1 0.572397
  2 2.21927
  3 5.3071
  4 8.90839
  5 14.4957
  6 20.5578
  7 28.391
  8 36.5239
e

set output
