#!/bin/bash

# check command-line args
if [[ $# != 1 ]] && [[ $# != 2 ]]; then
    echo "Usage: ./makeplot.bash file_id [nosim]"
    exit 1
fi

# set path
file_id=$1
output_dir=$DATADIR/$file_id
gp_script=mfpa
tex_script=fix1
nosim=n

# check if nosim
if [[ $# == 2 ]]; then
    if [ $2 == "nosim" ]; then
	nosim="y"
    fi	
fi

if [ $nosim == "y" ]; then
    gp_script="$gp_script-no-sim"
    tex_script="$tex_script-no-sim"
fi

# creating gnuplot script
echo '... creating gnuplot script'
cat $output_dir/$file_id.gp > tmp.gp
sed -e s:FILE_ID:$output_dir/$file_id:g $gp_script.gp >> tmp.gp

# make figs
echo '... running gnuplot'
gnuplot tmp.gp > /dev/null
mv -f tmp.gp $output_dir/$file_id.gp

# ps -> eps
echo '... converting figures to eps'
if [ $nosim != "y" ]; then
    (ps2epsi sim.ps sim.eps)
fi
(ps2epsi mf.ps mf.eps) &
(ps2epsi pa-di.ps pa-di.eps) &
(ps2epsi pa-dib.ps pa-dib.eps) &

# number of figs
if [ $nosim == "y" ]; then
    nfigs=6
else
    nfigs=10
fi

for i in $(seq 1 1 $nfigs)
  do  
  (ps2epsi Cxx$i.ps Cxx$i.eps) &
done

wait

# fix fonts
echo '... creating ps file'
latex $tex_script.tex > /dev/null 
dvips -q -o mfpa.ps $tex_script.dvi > /dev/null

# fix lines
echo '... fix lines'
sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" mfpa.ps > tmp.ps
mv -f tmp.ps mfpa.ps

# copy files
mv -f mfpa.ps $output_dir/$file_id.ps

# clean
echo '... cleaning'
rm -f $tex_script.log $tex_script.aux $tex_script.dvi 
#if [ $nosim != "y" ]; then
#    rm -f sim.ps  sim.eps 
#fi
#rm -f mf.ps  mf.eps 
#rm -f pa-di.ps  pa-di.eps
#%rm -f pa-dib.ps  pa-dib.eps
#for i in $(seq 1 1 $nfigs)
#  do  
#  rm -f Cxx$i.ps Cxx$i.eps
#done

echo 'done'
echo

#gv $output_dir/$file_id.ps &

exit 0
