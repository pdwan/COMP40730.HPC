#!/bin/bash 

gnuPlotTitleMatrix="Matrix size"
gnuPlotTitleBlock="Block size"
localSrcFileName="'demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat"
localGraphFileName="graphDir/pdwan-20140708.163544-graph-A1-Sijk-1D-0.png"

/usr/bin/gnuplot  << __EOF
unset log 
unset label 
set xtic auto 
set ytic auto 
set size 1,1 
set grid 
set key outside 
set output 'graphDir/pdwan-20140708.163544-graph-A1-Sijk-1D-0.png'

# 
set title '\${gnuPlotTitleMatrix} - comparison' 
set ylabel 'Time taken / ms' 
set xlabel 'Matrix size' 
set size 1,1 
set origin 0,0 
plot 'demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat' u 1:3 t 'simple' w l lw 0.5 lc rgb 'blue',  'demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat' u 1:4 t 'complex' w l lw 0.5 lc rgb 'red', 'demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat' u 1:5 t 'dgemm' w l lw 0.5 lc rgb 'black' 
# 
#set output '\${2.localGraphFileName}'
#set title '\${gnuPlotTitleblock} - comparison' 
#set title 'BLOCK - comparison' 
#set ylabel 'Time taken / ms' 
#set xlabel 'Block Size' 
#set size 0.5,0.5 
#set origin 0.5,0 
#set key outside 
#plot 'demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat' u 2:3 t 'simple' w l lw 0.5 lc rgb 'blue', 'demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat' u 2:4 t 'complex' w l lw 0.5 lc rgb 'red', 'demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat' u 2:5 t 'dgemm' w l lw 0.5 lc rgb 'black' 
# 
unset multiplot 
reset
pause -1
__EOF

