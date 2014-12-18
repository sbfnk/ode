#!/bin/sh

sed 's/TEST/'"$1"'.eps/' fix1.tex | latex
dvips -o $1.ps texput.dvi
ps2epsi $1.ps
mv $1.epsi $1.eps
rm $1.ps texput.dvi texput.log texput.aux
#rm fix1.log fix1.dvi fix1.aux
