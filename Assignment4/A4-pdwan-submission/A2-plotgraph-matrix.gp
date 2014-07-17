# Paula Dwan : Assignment 4
reset 
set xtic auto 
set ytic auto
set size 1,1 
set grid 
set key outside 
# 
set title 'Comparison : Matrix Size v time taken' 
set ylabel 'Time taken / s' 
set xlabel 'Matrix size' 
set origin 0,0 
plot 'logDir/<filename>.dat' u 1:2 t 'mpi' w l lw 0.5 lc rgb 'blue', 'logDir/<filename>.dat' u 1:4 t 'dgemm' w l lw 0.5 lc rgb 'red', 'logDir/<filename>.dat' u 1:6 t 'manual' w l lw 0.5 lc rgb 'green'  
# 
pause -1
