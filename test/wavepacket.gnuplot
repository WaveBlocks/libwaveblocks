# transparent
set terminal pngcairo size 1024,768

set title "diagram"
#set samples 1000
set hidden3d
set dgrid3d 100,100 splines

set output "wavepacket_real.png"
splot "wavepacket.csv" using 1:2:3 with lines

set output "wavepacket_imag.png"
splot "wavepacket.csv" using 1:2:4 with lines