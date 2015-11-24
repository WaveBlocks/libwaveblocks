set term epslatex color
set key top left
set xlabel "number of dimensions"
set ylabel "runtime [ms]"
set log y
set grid
set xtics 1
set output "speedup-genz-keister.tex"
plot "-" with lp pt 7 title "Python", "-" with lp pt 2 lc 3 title "C++"
  1 1.81172499657
  2 21.8601160049
  3 1012.47501373
  4 221381.198883
e
  1 0.0218638
  2 0.323561
  3 46.9158
  4 10843.9
e

set output
