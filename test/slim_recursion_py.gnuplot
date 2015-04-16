# transparent
set terminal pngcairo size 1024,768

set title "diagram"
#set samples 1000
set hidden3d
set dgrid3d 100,100 splines

set output "slim_recursion_py_real.png"
splot "slim_recursion_py.csv" using 1:2:3 with lines

set output "slim_recursion_py_imag.png"
splot "slim_recursion_py.csv" using 1:2:4 with lines