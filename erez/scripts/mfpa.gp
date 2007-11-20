set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "time" 0.000000,0.000000
set ylabel "" 0.000000,0.000000

set style line 1 lt 1 lc rgb "blue"  lw 2
set style line 2 lt 1 lc rgb "red"   lw 2
set style line 3 lt 1 lc rgb "green" lw 2
set style line 4 lt 2 lc rgb "blue"  lw 3
set style line 5 lt 2 lc rgb "red"   lw 3
set style line 6 lt 2 lc rgb "green" lw 3

set format y "% g"

set key right top samplen 3 spacing 1.5
set xrange [0:Tmax]

#####################################

set output 'sim.ps'
set title "Simulation"
plot 'FILE_ID.sim.dat' u 1:2 title 'S-' w l ls 1 \
   , 'FILE_ID.sim.dat' u 1:3 title 'I-' w l ls 2 \
   , 'FILE_ID.sim.dat' u 1:4 title 'R-' w l ls 3 \
   , 'FILE_ID.sim.dat' u 1:5 title 'S+' w l ls 4 \
   , 'FILE_ID.sim.dat' u 1:6 title 'I+' w l ls 5 \
   , 'FILE_ID.sim.dat' u 1:7 title 'R+' w l ls 6

#####################################

set output 'mf.ps'
set title "Mean Field"
plot 'FILE_ID.mf.dat' u 1:2 title 'S-' w l ls 1 \
   , 'FILE_ID.mf.dat' u 1:3 title 'I-' w l ls 2 \
   , 'FILE_ID.mf.dat' u 1:4 title 'R-' w l ls 3 \
   , 'FILE_ID.mf.dat' u 1:5 title 'S+' w l ls 4 \
   , 'FILE_ID.mf.dat' u 1:6 title 'I+' w l ls 5 \
   , 'FILE_ID.mf.dat' u 1:7 title 'R+' w l ls 6

#####################################

set output 'pa-di.ps'
set title "Pair Approx di"
plot 'FILE_ID.di.dat' u 1:2 title 'S-' w l ls 1 \
   , 'FILE_ID.di.dat' u 1:3 title 'I-' w l ls 2 \
   , 'FILE_ID.di.dat' u 1:4 title 'R-' w l ls 3 \
   , 'FILE_ID.di.dat' u 1:5 title 'S+' w l ls 4 \
   , 'FILE_ID.di.dat' u 1:6 title 'I+' w l ls 5 \
   , 'FILE_ID.di.dat' u 1:7 title 'R+' w l ls 6

#####################################

set output 'pa-dib.ps'
set title "Pair Approx dib"
plot 'FILE_ID.dib.dat' u 1:2 title 'S-' w l ls 1 \
   , 'FILE_ID.dib.dat' u 1:3 title 'I-' w l ls 2 \
   , 'FILE_ID.dib.dat' u 1:4 title 'R-' w l ls 3 \
   , 'FILE_ID.dib.dat' u 1:5 title 'S+' w l ls 4 \
   , 'FILE_ID.dib.dat' u 1:6 title 'I+' w l ls 5 \
   , 'FILE_ID.dib.dat' u 1:7 title 'R+' w l ls 6

#####################################

# correlations di model

#####################################

set output 'Cxx1.ps'	
set title "Correlation di model"
set xrange [0:Tmax]
set yrange [0:2.5]
set key right top samplen 3 spacing 2.5
	
plot 'FILE_ID.di.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w l \
   , 'FILE_ID.di.dat' u 1:(N*$12/($2*$6)/Qd) title 'CSid' w l \
   , 'FILE_ID.di.dat' u 1:(N*$16/($5*$3)/Qd) title 'CsId' w l \
   , 'FILE_ID.di.dat' u 1:(N*$24/($5*$6)/Qd) title 'Csid' w l \

#####################################

set output 'Cxx2.ps'
set title "Correlation di model"
set xrange [0:Tmax]
#set yrange [0:2.5]
#set key right top samplen 2 spacing 2.5

plot 'FILE_ID.di.dat' u 1:(N*$30/($2*$3)/Qi) title 'CSIi' w l \
   , 'FILE_ID.di.dat' u 1:(N*$33/($2*$6)/Qi) title 'CSii' w l \
   , 'FILE_ID.di.dat' u 1:(N*$37/($5*$3)/Qi) title 'CsIi' w l \
   , 'FILE_ID.di.dat' u 1:(N*$45/($5*$6)/Qi) title 'Csii' w l \

#####################################

set output 'Cxx3.ps'
set title "Correlation di model"
set xrange [0:Tmax]
#set yrange [0:2.5]
#set key right top samplen 2 spacing 2.5

Qd2=Qd/2
Qi2=Qi/2

plot 'FILE_ID.di.dat' u 1:(N*$8/($2*$2)/Qd2) title 'CSSd' w l \
   , 'FILE_ID.di.dat' u 1:(N*$11/($2*$5)/Qd) title 'CSsd' w l \
   , 'FILE_ID.di.dat' u 1:(N*$23/($5*$5)/Qd2) title 'Cssd' w l \
   , 'FILE_ID.di.dat' u 1:(N*$29/($2*$2)/Qi2) title 'CSSi' w l \
   , 'FILE_ID.di.dat' u 1:(N*$32/($2*$5)/Qi) title 'CSsi' w l \
   , 'FILE_ID.di.dat' u 1:(N*$44/($5*$5)/Qi2) title 'Cssi' w l \

#####################################

set output 'Cxx4.ps'
set title "Correlation di model"
set xrange [0:Tmax]
#set yrange [0:2.5]
#set key right top samplen 2 spacing 2.5

plot 'FILE_ID.di.dat' u 1:(N*$14/($3*$3)/Qd2) title 'CIId' w l \
   , 'FILE_ID.di.dat' u 1:(N*$17/($3*$6)/Qd) title 'CIid' w l \
   , 'FILE_ID.di.dat' u 1:(N*$26/($6*$6)/Qd2) title 'Ciid' w l \
   , 'FILE_ID.di.dat' u 1:(N*$35/($3*$3)/Qi2) title 'CIIi' w l \
   , 'FILE_ID.di.dat' u 1:(N*$38/($3*$6)/Qi) title 'CIii' w l \
   , 'FILE_ID.di.dat' u 1:(N*$47/($6*$6)/Qi2) title 'Ciii' w l

#####################################

# correlations dib model

#####################################

set output 'Cxx5.ps'	
set title "Correlation dib model"
set xrange [0:Tmax]
set yrange [0:2.5]
#set key right top samplen 3 spacing 2.5

Qb2=Qb/2
	
plot 'FILE_ID.dib.dat' u 1:(N*$51/($2*$3)/Qb) title 'CSIb' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$54/($2*$6)/Qb) title 'CSib' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$58/($5*$3)/Qb) title 'CsIb' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$66/($5*$6)/Qb) title 'Csib' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSIdsim' w l \

#####################################

set output 'Cxx6.ps'
set title "Correlation dib model"
set xrange [0:Tmax]
#set yrange [0:2.5]
#set key right top samplen 2 spacing 2.5

plot 'FILE_ID.dib.dat' u 1:(N*$50/($2*$2)/Qb2) title 'CSSb' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$53/($2*$5)/Qb) title 'CSsb' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$65/($5*$5)/Qb2) title 'Cssb' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$56/($3*$3)/Qb2) title 'CIIb' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$59/($3*$6)/Qb) title 'CIib' w l \
   , 'FILE_ID.dib.dat' u 1:(N*$68/($6*$6)/Qb2) title 'Ciib' w l \

#####################################

# correlations simulation

#####################################

set output 'Cxx7.ps'	
set title "Correlation sim"
set xrange [0:Tmax]
set yrange [0:2.5]
set key right top samplen 3 spacing 2.5
	
plot 'FILE_ID.sim.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$12/($2*$6)/Qd) title 'CSid' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$16/($5*$3)/Qd) title 'CsId' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$24/($5*$6)/Qd) title 'Csid' w l \

#####################################

set output 'Cxx8.ps'
set title "Correlation sim"

plot 'FILE_ID.sim.dat' u 1:(N*$30/($2*$3)/Qi) title 'CSIi' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$33/($2*$6)/Qi) title 'CSii' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$37/($5*$3)/Qi) title 'CsIi' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$45/($5*$6)/Qi) title 'Csii' w l \

#####################################

set output 'Cxx9.ps'
set title "Correlation sim"

plot 'FILE_ID.sim.dat' u 1:(N*$8/($2*$2)/Qd2) title 'CSSd' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$11/($2*$5)/Qd) title 'CSsd' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$23/($5*$5)/Qd2) title 'Cssd' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$29/($2*$2)/Qi2) title 'CSSi' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$32/($2*$5)/Qi) title 'CSsi' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$44/($5*$5)/Qi2) title 'Cssi' w l \

#####################################

set output 'Cxx10.ps'
set title "Correlation sim"

plot 'FILE_ID.sim.dat' u 1:(N*$14/($3*$3)/Qd2) title 'CIId' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$17/($3*$6)/Qd) title 'CIid' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$26/($6*$6)/Qd2) title 'Ciid' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$35/($3*$3)/Qi2) title 'CIIi' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$38/($3*$6)/Qi) title 'CIii' w l \
   , 'FILE_ID.sim.dat' u 1:(N*$47/($6*$6)/Qi2) title 'Ciii' w l

#####################################



quit




