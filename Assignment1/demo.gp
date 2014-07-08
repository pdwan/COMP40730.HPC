set title "Demo data" 
set xlabel "Message size"
set ylabel "Time (sec)"
#set legend position
set key left box

plot "datafile" using 1:2 title 'Title1' with lines, "datafile2.dat" using 1:2 title 'Title2' with lines
#it will close immediately without pause
pause -1

===========
gnuPlotTitleMatrix="Matrix size"
gnuPlotTitleBlock="Block size"
localSrcFileName="demoDir/pdwan-20140708.163544-demo-A1-Sijk-1D-0.dat"
localGraphFileName="graphDir/pdwan-20140708.163544-graph-A1-Sijk-1D-0.png"

unset log 
unset label 
set xtic auto 
set ytic auto 
set size 1,1 
set grid 
set multiplot 
set key outside 
# 
set title '${gnuPlotTitleMatrix} - comparison' 
set ylabel 'Time taken / ms' 
set xlabel 'Matrix size' 
set size 0.5,0.5 
set origin 0,0 
plot '${localSrcFileName}' u 1:3 t 'simple' w l lw 0.5 lc rgb 'blue', \ 
  '${localSrcFileName}' u 1:4 t 'complex' w l lw 0.5 lc rgb 'red', \
  '${localSrcFileName}' u 1:5 t 'dgemm' w l lw 0.5 lc rgb 'black' 
# 
set title '${gnuPlotTitleblock} - comparison' 
set title 'BLOCK - comparison' 
set ylabel 'Time taken / ms' 
set xlabel 'Block Size' 
set size 0.5,0.5 
set origin 0.5,0 
set key outside 
plot '${localSrcFileName}' u 2:3 t 'simple' w l lw 0.5 lc rgb 'blue', \
  '${localSrcFileName}' u 2:4 t 'complex' w l lw 0.5 lc rgb 'red', \
  '${localSrcFileName}' u 2:5 t 'dgemm' w l lw 0.5 lc rgb 'black' 
# 
unset multiplot 
pause -1
