#!/bin/sh

./bin/solve_ode.x --susfunc 0 $* > /dev/null
first=$(awk -f maxtime.awk try1.dat)
./bin/solve_ode.x --susfunc 1 $* > /dev/null
second=$(awk -f maxtime.awk try1.dat)
echo $(echo "$first - $second" | bc) $first $second

