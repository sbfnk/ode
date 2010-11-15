#!/bin/sh

./bin/net_ode --susfunc 0 $* > /dev/null
first=$(awk -f maxtime.awk try1.dat)
./bin/net_ode --susfunc 1 $* > /dev/null
second=$(awk -f maxtime.awk try1.dat)
echo $(echo "$first - $second" | bc) $first $second

