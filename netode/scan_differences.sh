#!/bin/zsh

degree_dist=$1

for ((i=1;$i<100;i=$i+1)) 
do
  T=$(echo "scale=2;$i/100" | bc)
  echo $(../percolation/bin/percolation -t single -d $degree_dist -T $T | grep "R_0" | awk '{print $3}') $(./compare_nb.sh -o params/ode.prm  -T $T -d $degree_dist --ic-file init/net_ode.ic --check-convergence 2> /dev/null | awk '{print $1}')
done
