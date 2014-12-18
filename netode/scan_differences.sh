#!/bin/zsh

degree_dist=$1

for ((i=1;$i<100;i=$i+1)) 
do
  T=$(echo "scale=2;$i/100" | bc)
  beta=$(echo "scale=2;1/(1-$T)" | bc)
  compfull=$(./compare_nb.sh -o params/ode.prm  -T $T -d $degree_dist --ic-file init/net_ode.ic --check-convergence 2> /dev/null)
  comp=$(echo $compfull | awk '{print $1}')
  compscaled=$(echo "$comp*$beta" | bc 2> /dev/null)
  echo $(../percolation/bin/percolation -t single -d $degree_dist -T $T | grep "R_0" | awk '{print $3}') $compscaled $compfull
done
