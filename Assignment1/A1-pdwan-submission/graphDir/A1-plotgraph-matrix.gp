# Paula Dwan : Assignment 1
reset 
set xtic auto 
set ytic auto
set size 1,1 
set grid 
set key outside 
# 
set title 'Comparison : Matrix size v Time taken' 
set ylabel 'Time taken / s' 
set xlabel 'Matrix size : NxN' 
set origin 0,0 
set key outside 
plot 'logDir/pdwan-<filename>.dat' u 2:3 t 'manual simple' w l lw 0.8 lc rgb 'blue', 'logDir/pdwan-<filename>.dat' u 2:4 t 'manual complex' w l lw 0.8 lc rgb 'black', 'logDir/pdwan-<filename>.dat' u 2:5 t 'dgemm' w l lw 0.8 lc rgb 'red' 
# 
pause -1
