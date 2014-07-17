# Paula Dwan : Assignment 1
reset 
set xtic auto 
set ytic auto
set size 1,1 
set grid 
set key outside 
# 
set title 'Straight-forward IJK Comparison : Matrix size v Time taken' 
set ylabel 'Time taken / seconds' 
set xlabel 'Matrix size : NxN' 
#set xrange [40:110]
#set yrange [0:0.005]
#set ytics (0.001,0.002,0.003,0.004,0.005)
# set xtics (10,20,30,40,50,60,70,80,90,100,110)
set xrange [0:110]
#set yrange [0:0.05]

set origin 0,0 
set key outside 
plot 'logDir/Sijk.dat' u 1:2 t 'manual simple' w l lw 0.8 lc rgb 'blue', 'logDir/Sijk.dat' u 1:3 t 'manual complex' w l lw 0.8 lc rgb 'black', 'logDir/Sijk.dat' u 1:4 t 'dgemm' w l lw 0.8 lc rgb 'red' 
# 
pause -1
