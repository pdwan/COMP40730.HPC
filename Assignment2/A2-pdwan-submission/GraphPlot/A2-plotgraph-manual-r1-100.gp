# Paula Dwan : Assignment 2
reset 
set xrange [0:100]
set xtic (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
 set yrange [0:0.012]
set ytic (0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012)
set size 1,1 
set grid 
set key outside 
# 
set title 'pThread Comparison : Matrix Size -v- Time Taken' 
set ylabel 'Time taken / seconds' 
set xlabel 'Matrix size / NxN' 
set origin 0,0 
plot '../logDir/manual-100-i.dat' u 1:2 t 'manual' w l lw 0.5 lc rgb 'blue', '../logDir/manual-100-i.dat' u 1:4 t 'dgemm' w l lw 0.5 lc rgb 'red'
# 
pause -1
