# Paula Dwan : Assignment 4
reset 
set xrange [0:100]
set xtic (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
set yrange [0:0.01]
#set ytic (0, 0.02, 0.04, 0.06, 0.08, 0.1)
# set xtic auto
#set ytic auto
set size 1,1 
set grid 
set key outside 
# 
set title 'Comparison : Matrix Size v time taken' 
set ylabel 'Time taken / s' 
set xlabel 'Matrix size' 

set origin 0,0 
plot 'logDir/manual-i-100.dat' u 1:2 t 'manual' w l lw 0.5 lc rgb 'blue', 'logDir/manual-i-100.dat' u 1:4 t 'dgemm' w l lw 0.5 lc rgb 'red'

# 
pause -1
