# Paula Dwan : Assignment 1
reset 
set xtic auto 
set ytic auto
set size 1,1 
set grid 
set key outside 
# 
set title 'Blocked IJK Comparison : Block size -v- Time taken' 
set ylabel 'Time taken / seconds' 
set xlabel 'Block size : BxB' 
# set xrange [0:25]
#set yrange [0:0.01]
#set ytics (0.001,0.002,0.003,0.004,0.005)
# set xtics (10,20,30,40,50,60,70,80,90,100,110)
#set xrange [0:1100]
#set xtics (100,200,300,400,500,600,700,800,900,1000,1100)
# set yrange [0:25]
#set ytics (5,10,15,20,25,30,35,40)
set xrange [0:30]
#set yrange [0:0.05]
set origin 0,0 
set key outside 
plot 'logDir/Bijk.dat' u 2:3 t 'manual simple' w l lw 0.8 lc rgb 'blue', 'logDir/Bijk.dat' u 2:4 t 'manual complex' w l lw 0.8 lc rgb 'black', 'logDir/Bijk.dat' u 2:5 t 'dgemm' w l lw 0.8 lc rgb 'red' 
# 
pause -1
