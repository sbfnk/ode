set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "time" ; #0.000000,0.000000
set ylabel "" ; #0.000000,0.000000

set style line 1 lt 1 lc rgb "blue"  lw 2
set style line 2 lt 1 lc rgb "red"   lw 2
set style line 3 lt 1 lc rgb "green" lw 2
set style line 4 lt 2 lc rgb "blue"  lw 3
set style line 5 lt 2 lc rgb "red"   lw 3
set style line 6 lt 2 lc rgb "green" lw 3

set format y "%3.1e"

set key right top samplen 2 spacing 1.5
set xrange [0:20]

#####################################

set output 'ode.ps'
set title "plot"
plot 'FILE_ID.dat' u 1:2 title 'S' w l ls 1 \
   , 'FILE_ID.dat' u 1:3 title 'I' w l ls 2 \
   , 'FILE_ID.dat' u 1:4 title 'J' w l ls 3 
#   , 'FILE_ID.dat' u 1:5 title 'Screened' w l ls 4 \

#####################################

quit



