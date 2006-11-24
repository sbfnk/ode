set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "time" 0.000000,0.000000
set ylabel "" 0.000000,0.000000

#####################################

set output 'mf.ps'
set title "Mean Filed"
plot 'try1.mf.dat' u 1:2 title 'S' w lp \
   , 'try1.mf.dat' u 1:3 title 'I' w lp \
   , 'try1.mf.dat' u 1:4 title 'R' w lp \
   , 'try1.mf.dat' u 1:5 title 's' w lp \
   , 'try1.mf.dat' u 1:6 title 'i' w lp \
   , 'try1.mf.dat' u 1:7 title 'r' w lp

#####################################

set output 'pa.ps'
set title "Pair Approx"
plot 'try1.pa.dat' u 1:2 title 'S' w lp \
   , 'try1.pa.dat' u 1:3 title 'I' w lp \
   , 'try1.pa.dat' u 1:4 title 'R' w lp \
   , 'try1.pa.dat' u 1:5 title 's' w lp \
   , 'try1.pa.dat' u 1:6 title 'i' w lp \
   , 'try1.pa.dat' u 1:7 title 'r' w lp

quit


#plot 'panelA/ee1.dat' u 1:4 w l ls 1 \
#   , 'panelA/wdots_lin.dat' u 1:2 w p 7 7 \
#   , x/4 w l ls 1
#
#####################################

set noytics
set noxtics
set noy2tics
set noxlabel
set noylabel
set key 6.3,.75 samplen 2 spacing 1.3


### panel a ##################################
set size .45,0.6
set origin 0.15,0.6
set ytics ("0" 0)

set xlabel "xx" 0.000000,0.000000
set xrange [0:7.5]
set yrange [-.05:0.8]

set nolabel
set label 'a)' at 0.6,0.75 center

plot 'p.0.25.m2.6.take1.LG.cross.dat' u 1:2 title 'b1' w l ls 1 \
   , 'p.0.25.m2.6.take1.LG.cross.dat' u 1:3 title 'b3' w l ls 2 

### panel b ##################################
set size .45,0.6
set origin 0.6,0.6
set noytics

set nolabel
set label 'b)' at 0.6,0.75 center

plot 'p.0.6.m2.6.take1.LG.cross.dat' u 1:2 title 'b1' w l ls 1 \
   , 'p.0.6.m2.6.take1.LG.cross.dat' u 1:3 title 'b3' w l ls 2 

### panel c ##################################
set size .45,0.6
set origin 0.15,0.0
set ytics ("0" 0)

set nolabel
set label 'c)' at 0.6,0.75 center

plot 'p.0.25.m2.4.1.take1.LS.cross.dat' u 1:2 title 'b1' w l ls 1 \
   , 'p.0.25.m2.4.1.take1.LS.cross.dat' u 1:3 title 'b4' w l ls 2

### panel d ##################################
set size .45,0.6
set origin 0.6,0.0
set noytics

set nolabel
set label 'd)' at 0.6,0.75 center

plot 'p.0.6.m2.4.1.take1.LS.cross.dat' u 1:2 title 'b1' w l ls 1 \
   , 'p.0.6.m2.4.1.take1.LS.cross.dat' u 1:3 title 'b4' w l ls 2

quit


