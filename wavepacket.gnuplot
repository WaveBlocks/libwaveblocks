# transparent
set terminal pngcairo size 1024,768
set output "output.png"

set title "diagram"
#set samples 1000
set hidden3d
set dgrid3d 100,100 splines
splot "wavepacket.csv" with lines
print "all done"