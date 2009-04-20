#!/bin/bash

# check command-line args
if [[ $# != 1 ]]; then
    echo "Usage: ./makeplot.bash file_id"
    echo "       file_id is missing."
    exit 1
fi

gp_script=plot.gp

# make figs
echo '... running gnuplot'
sed -e s:FILE_ID:$1:g $gp_script > tmp.gp
gnuplot tmp.gp > /dev/null
rm tmp.gp

# ps -> eps
echo '... converting figures to eps'
ps2epsi ode.ps ode.eps

# fix fonts
echo '... creating ps file'
latex fix1.tex
dvips -q -o mfpa.ps fix1.dvi

# fix lines
echo '... fix lines'
sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" mfpa.ps > tmp.ps
mv -f tmp.ps ode.ps

# clean
echo '... cleaning'
rm -f fix1.log fix1.aux fix1.dvi 
rm -f ode.eps mfpa.ps
