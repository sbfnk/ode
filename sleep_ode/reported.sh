#!/bin/sh

cat $1 | awk 'BEGIN{prev=0;diff=0;curyear=0}{if (curyear == 0) {curyear = int($1);} else while ($1 >= (curyear + 1)) {if (prev>0) {diff=$4-prev} else {prev=$4}; print curyear" "$2" "$3" "$4" "diff;curyear = curyear + 1}}'
