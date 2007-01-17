#!/bin/bash

# check command-line args
if [[ $# != 1 ]]; then
    echo "Usage: ./makeplot.bash file_id"
    echo "       file_id is missing."
    exit 1
fi

gp_script=mfpa.gp
fixf=fix21.tex

# creating gnuplot script
echo '... creating gnuplot script'
sed -e s:FILE_ID:$1:g $gp_script > tmp.gp

# make figs
echo '... running gnuplot'
gnuplot tmp.gp > /dev/null
rm -f tmp.gp

# ps -> eps
echo '... converting figures to eps'
ps2epsi mf.ps mf.eps
ps2epsi pa.ps pa.eps
ps2epsi pa1.ps pa1.eps
ps2epsi pa2.ps pa2.eps
ps2epsi pa3.ps pa3.eps
ps2epsi pa4.ps pa4.eps
ps2epsi pa5.ps pa5.eps
ps2epsi pa6.ps pa6.eps
ps2epsi pa7.ps pa7.eps
ps2epsi pa8.ps pa8.eps
ps2epsi Cxx1.ps Cxx1.eps
ps2epsi Cxx2.ps Cxx2.eps
ps2epsi Cxx3.ps Cxx3.eps
ps2epsi Cxx4.ps Cxx4.eps

# fix fonts
echo '... creating ps file'
latex fix1.tex > /dev/null 
dvips -q -o mfpa.ps fix1.dvi > /dev/null

# clean
echo '... cleaning'
rm -f fix1.log fix1.aux fix1.dvi 
rm -f mf.ps  mf.eps 
rm -f pa.ps  pa.eps
rm -f pa1.ps pa1.eps 
rm -f pa2.ps pa2.eps 
rm -f pa3.ps pa3.eps 
rm -f pa4.ps pa4.eps
rm -f pa5.ps pa5.eps 
rm -f pa6.ps pa6.eps 
rm -f pa7.ps pa7.eps 
rm -f pa8.ps pa8.eps
rm -f Cxx1.ps Cxx1.eps
rm -f Cxx2.ps Cxx2.eps
rm -f Cxx3.ps Cxx3.eps
rm -f Cxx4.ps Cxx4.eps

# extracting single figures
#for fig in mf.eps pa.eps pa5.eps pa6.eps Cxx.eps
#do
#    echo "... $fig"
#    sed -e s/FIGURE/$fig/g $fixf > tmp.tex
#    latex tmp.tex > /dev/null  
#    dvips -q -o tmp.ps tmp.dvi > /dev/null
#    ps2epsi tmp.ps $1.$fig
#    rm -f tmp.dvi tmp.log tmp.aux tmp.ps tmp.tex $fig
#done	

# renaming
mv mfpa.ps $1.mfpa.ps
echo '... done'

# ghostview
gv $1.mfpa.ps &
