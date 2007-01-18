set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "time" 0.000000,0.000000
set ylabel "" 0.000000,0.000000

set format y "%3.1e"

#####################################

set output 'mf.ps'
set title "Mean Field"
plot 'FILE_ID.mf.dat' u 1:2 title 'S' w lp \
   , 'FILE_ID.mf.dat' u 1:3 title 'I' w lp \
   , 'FILE_ID.mf.dat' u 1:4 title 'R' w lp \
   , 'FILE_ID.mf.dat' u 1:5 title 's' w lp \
   , 'FILE_ID.mf.dat' u 1:6 title 'i' w lp \
   , 'FILE_ID.mf.dat' u 1:7 title 'r' w lp

#####################################

set output 'pa.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:2 title 'S' w lp \
   , 'FILE_ID.pa.dat' u 1:3 title 'I' w lp \
   , 'FILE_ID.pa.dat' u 1:4 title 'R' w lp \
   , 'FILE_ID.pa.dat' u 1:5 title 's' w lp \
   , 'FILE_ID.pa.dat' u 1:6 title 'i' w lp \
   , 'FILE_ID.pa.dat' u 1:7 title 'r' w lp

#####################################

set output 'pa1.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:8 title 'SSd' w lp \
   , 'FILE_ID.pa.dat' u 1:9 title 'SId' w lp \
   , 'FILE_ID.pa.dat' u 1:10 title 'SRd' w lp \
   , 'FILE_ID.pa.dat' u 1:14 title 'IId' w lp \
   , 'FILE_ID.pa.dat' u 1:15 title 'IRd' w lp \
   , 'FILE_ID.pa.dat' u 1:19 title 'RRd' w lp

#####################################

set output 'pa2.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:29 title 'SSi' w lp \
   , 'FILE_ID.pa.dat' u 1:30 title 'SIi' w lp \
   , 'FILE_ID.pa.dat' u 1:31 title 'SRi' w lp \
   , 'FILE_ID.pa.dat' u 1:35 title 'IIi' w lp \
   , 'FILE_ID.pa.dat' u 1:36 title 'IRi' w lp \
   , 'FILE_ID.pa.dat' u 1:40 title 'RRi' w lp

#####################################

set output 'pa3.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:23 title 'ssd' w lp \
   , 'FILE_ID.pa.dat' u 1:24 title 'sid' w lp \
   , 'FILE_ID.pa.dat' u 1:25 title 'srd' w lp \
   , 'FILE_ID.pa.dat' u 1:26 title 'iid' w lp \
   , 'FILE_ID.pa.dat' u 1:27 title 'ird' w lp \
   , 'FILE_ID.pa.dat' u 1:28 title 'rrd' w lp

#####################################

set output 'pa4.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:44 title 'ssi' w lp \
   , 'FILE_ID.pa.dat' u 1:45 title 'sii' w lp \
   , 'FILE_ID.pa.dat' u 1:46 title 'sri' w lp \
   , 'FILE_ID.pa.dat' u 1:47 title 'iii' w lp \
   , 'FILE_ID.pa.dat' u 1:48 title 'iri' w lp \
   , 'FILE_ID.pa.dat' u 1:49 title 'rri' w lp

#####################################

set output 'pa5.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:11 title 'Ssd' w lp \
   , 'FILE_ID.pa.dat' u 1:12 title 'Sid' w lp \
   , 'FILE_ID.pa.dat' u 1:13 title 'Srd' w lp \
   , 'FILE_ID.pa.dat' u 1:17 title 'Iid' w lp \
   , 'FILE_ID.pa.dat' u 1:18 title 'Ird' w lp \
   , 'FILE_ID.pa.dat' u 1:22 title 'Rrd' w lp

#####################################

set output 'pa6.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:32 title 'Ssi' w lp \
   , 'FILE_ID.pa.dat' u 1:33 title 'Sii' w lp \
   , 'FILE_ID.pa.dat' u 1:34 title 'Sri' w lp \
   , 'FILE_ID.pa.dat' u 1:38 title 'Iii' w lp \
   , 'FILE_ID.pa.dat' u 1:39 title 'Iri' w lp \
   , 'FILE_ID.pa.dat' u 1:43 title 'Rri' w lp

#####################################

set output 'pa7.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:16 title 'sId' w lp \
   , 'FILE_ID.pa.dat' u 1:20 title 'sRd' w lp \
   , 'FILE_ID.pa.dat' u 1:21 title 'iRd' w lp \

#####################################

set output 'pa8.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:37 title 'sIi' w lp \
   , 'FILE_ID.pa.dat' u 1:41 title 'sRi' w lp \
   , 'FILE_ID.pa.dat' u 1:42 title 'iRi' w lp \

#####################################

set output 'Cxx1.ps'
set title "Correlation"
set xrange [0:20]
#set yrange [0:1]
set key spacing 2

# try1:
#plot 'FILE_ID.pa.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$30/($2*$3)/Qi) title 'CSIi' w lp \
#   , 'FILE_ID.pa.dat' u 1:($3/N) title 'I' w lp \

# try2:
#plot 'FILE_ID.pa.dat' u 1:(N*$11/($2*$5)/Qd) title 'CSsd' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$32/($2*$5)/Qi) title 'CSsi' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w lp \

# try3:
#plot 'FILE_ID.pa.dat' u 1:(N*$11/($2*$5)/Qd) title 'CSsd' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$32/($2*$5)/Qi) title 'CSsi' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w lp \
#   , 'FILE_ID.pa.dat' u 1:($5/N) title 's' w lp \

# try4:
#plot 'FILE_ID.pa.dat' u 1:(N*$16/($5*$3)/Qd) title 'CsId' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$32/($2*$5)/Qi) title 'CSsi' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w lp \
#   , 'FILE_ID.pa.dat' u 1:($3/N) title 'I' w lp \

# try5
#plot 'FILE_ID.pa.dat' u 1:(N*$24/($5*$6)/Qd) title 'Csid' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$45/($5*$6)/Qi) title 'Csii' w lp \
#   , 'FILE_ID.pa.dat' u 1:($5/N) title 's' w lp \
#   , 'FILE_ID.pa.dat' u 1:($6/N) title 'i' w lp \

# try6
#plot 'FILE_ID.pa.dat' u 1:(N*$24/($5*$6)/Qd) title 'Csid' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$45/($5*$6)/Qi) title 'Csii' w lp \
#   , 'FILE_ID.pa.dat' u 1:($3/N) title 'I' w lp \
#   , 'FILE_ID.pa.dat' u 1:($5/N) title 's' w lp \
#   , 'FILE_ID.pa.dat' u 1:($6/N) title 'i' w lp \

# try7
#plot 'FILE_ID.pa.dat' u 1:(N*$24/($5*$6)/Qd) title 'Csid' w lp \
#   , 'FILE_ID.pa.dat' u 1:(N*$45/($5*$3)/Qd) title 'CsId' w lp \
#   , 'FILE_ID.pa.dat' u 1:(($3+$6)/N) title 'I+i' w lp \

# try10
plot 'FILE_ID.pa.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$12/($2*$6)/Qd) title 'CSid' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$16/($5*$3)/Qd) title 'CsId' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$24/($5*$6)/Qd) title 'Csid' w lp \
   , 'FILE_ID.pa.dat' u 1:(($3+$6)/N) title 'I+i' w lp \

set output 'Cxx2.ps'
set title "Correlation"
set xrange [0:20]
set yrange [0:2.5]
set key spacing 2

plot 'FILE_ID.pa.dat' u 1:(N*$30/($2*$3)/Qi) title 'CSIi' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$33/($2*$6)/Qi) title 'CSii' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$37/($5*$3)/Qi) title 'CsIi' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$45/($5*$6)/Qi) title 'Csii' w lp \
   , 'FILE_ID.pa.dat' u 1:(($3+$6)/N) title 'I+i' w lp \

set output 'Cxx3.ps'
set title "Correlation"
set xrange [0:20]
#set yrange [0:2.5]
set key spacing 2

plot 'FILE_ID.pa.dat' u 1:(N*$8/($2*$2)/Qd) title 'CSSd' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$11/($2*$5)/Qd) title 'CSsd' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$23/($5*$5)/Qd) title 'Cssd' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$32/($2*$5)/Qi) title 'CSsi' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$44/($5*$5)/Qi) title 'Cssi' w lp \

set output 'Cxx4.ps'
set title "Correlation"
set xrange [0:20]
set yrange [0:2.5]
set key spacing 2

plot 'FILE_ID.pa.dat' u 1:(N*$14/($3*$3)/Qd) title 'CIId' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$17/($3*$6)/Qd) title 'CIid' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$26/($6*$6)/Qd) title 'Ciid' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$35/($3*$3)/Qi) title 'CIIi' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$38/($3*$6)/Qi) title 'CIii' w lp \
   , 'FILE_ID.pa.dat' u 1:(N*$47/($6*$6)/Qi) title 'Ciii' w lp \

quit



