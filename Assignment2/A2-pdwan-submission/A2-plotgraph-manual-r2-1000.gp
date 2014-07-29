# Paula Dwan : Assignment 2
reset 
set xrange [0:1000]
set xtic (0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
# set yrange [0:0.01]
#set ytic (0, 0.02, 0.04, 0.06, 0.08, 0.1)
set ytic auto
set size 1,1 
set grid 
set key outside 
# 
set title 'pThread Comparison : Matrix Size -v- Time Taken' 
set ylabel 'Time taken / seconds' 
set xlabel 'Matrix size / NxN' 
set origin 0,0 
plot 'logDir/manual-1000-r.dat' u 1:2 t 'manual' w l lw 0.5 lc rgb 'blue', 'logDir/manual-1000-r.dat' u 1:4 t 'dgemm' w l lw 0.5 lc rgb 'red'
# 
pause -1
