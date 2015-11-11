set term epslatex color
set key top left
set xlabel "number of dimensions"
set ylabel "runtime [ms]" rotate by 0
set log y
set grid
set xtics 1
set output "speedup-multi-dim.tex"
plot "-" with lp pt 7 title "Python", "-" with lp pt 2 lc 3 title "C++"
  1 1.75928559303
  2 22.5061798096
  3 2763.9549017
  4 2051378.05891
e
  1 0.0206635
  2 0.400831
  3 120.666
  4 101581
e

set output
