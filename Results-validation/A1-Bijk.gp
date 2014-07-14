# To execute, launch GNUpopt and run :
# gnuplot> load <filename.gp>
# making sure that <data file>.dat used is updated as required.

unset log
unset label
set xtic auto
set ytic auto
set size 1,1
set grid
set multiplot
set key outside
#
set title 'MATRIX - comparison'
set ylabel 'Time taken / ms'
set xlabel 'Matrix size'
set size 0.5,0.5
set origin 0,0
plot 'sample-A1-Bijk.dat' u 1:3 t 'simple' w l lw 0.5 lc rgb 'blue', 'sample-A1-Bijk.dat' u 1:4 t 'complex' w l lw 0.5 lc rgb 'red',  'sample-A1-Bijk.dat' u 1:5 t 'dgemm' w l lw 0.5 lc rgb 'black'
#
set title 'BLOCK - comparison'
set ylabel 'Time taken / ms'
set xlabel 'Block Size'
set size 0.5,0.5
set origin 0.5,0
set key outside
plot 'sample-A1-Bijk.dat' u 2:3 t 'simple' w l lw 0.5 lc rgb 'blue', 'sample-A1-Bijk.dat' u 2:4 t 'complex' w l lw 0.5 lc rgb 'red', 'sample-A1-Bijk.dat' u 2:5 t 'dgemm' w l lw 0.5 lc rgb 'black'
#
unset multiplot
pause -1
