set term epslatex color
set key top left
set xlabel "number of coefficients"
set ylabel "runtime [ms]" rotate by 0
set log x 2
set log y
set grid
set output "speedup-one-dim.tex"
plot "-" with lp pt 7 title "Python", "-" with lp pt 2 lc 3 title "C++"
  64   10.2753400803
  128  20.5673599243
  256  42.4106502533
  512  90.8460903168
  1024 207.515687943
  2048 530.075509548
e
  64   0.205799
  128  0.298305
  256  0.819514
  512  4.44001
  1024 16.5571
  2048 57.2452
e

set output
