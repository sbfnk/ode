#!/bin/bash

# make figs
gnuplot mfpa.gp

# ps -> eps
ps2epsi mf.ps mf.eps
ps2epsi pa.ps pa.eps

# fix fonts
latex fix1.tex; dvips -o mfpa.ps fix1.dvi

# ps -> eps
ps2epsi mfpa.ps mfpa.eps

# group figures
#latex group_figs.tex
#dvips -o $fname.ps group_figs.dvi

# clean
rm -f mfpa.ps mf.ps pa.ps mf.eps pa.eps fix1.log fix1.aux fix1.dvi

# ghostview
gv mfpa.eps 
