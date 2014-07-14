# Paula Dwan : Assignment 3
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
plot 'logDir/pdwan-20140714.052102-data-A3-omp-1D.dat' u 1:3 t 'simple' w l lw 0.5 lc rgb 'blue', 'logDir/pdwan-20140714.052102-data-A3-omp-1D.dat' u 1:5 t 'dgemm' w l lw 0.5 lc rgb 'red' 
# 
pause -1
