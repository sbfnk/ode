#!/bin/sh

cat $1 | awk 'BEGIN{prev=0;prev2=0;diff=0;diff2=0;curyear=0}{if (curyear == 0) {curyear = int($1);} else while ($1 >= (curyear + 1)) {diff=$4-prev; prev=$4; diff2=$5-prev2;prev2=$5; print curyear" "$2" "$3" "$4" "diff" "$5" "diff2; curyear=curyear+1}}'
