#!/bin/sh

firstfile=$1;
secondfile=$2;

lqstring="0"

if [ ! -e "$firstfile" -o ! -e "$secondfile" ];
then
  echo Must specify two files
  exit 0
fi

while read line
do
  index=$(echo $line | awk '{print $1}')
  value=$(echo $line | awk '{print $5}')
  value2=$(grep $index $secondfile | awk '{print $2}')
  lqstring="$lqstring + ($value - $value2)^2"
done < $firstfile

echo $lqstring | bc
