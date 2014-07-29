# Paula Dwan : Assignment 2
reset 
set xrange [0:1000]
set xtic (0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
set yrange [0:25]
set ytic (0, 5, 10,15, 20, 25)
set size 1,1 
set grid 
set key outside 
# 
set title 'pThread Comparison : Matrix Size -v- Time Taken' 
set ylabel 'Time taken / seconds' 
set xlabel 'Matrix size / NxN' 
set origin 0,0 
plot '../logDir/manual-1000-i.dat' u 1:2 t 'manual' w l lw 0.5 lc rgb 'blue', '../logDir/manual-1000-i.dat' u 1:4 t 'dgemm' w l lw 0.5 lc rgb 'red'
# 
pause -1
