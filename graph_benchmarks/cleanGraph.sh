#!/bin/bash
cat $1 | while read line
do
	what=`echo $line | awk '{split($0,a," "); print a[1]}'`
	what2=`echo $line | awk '{split($0,a," "); print a[2]}'`
	echo $what $what2 >> $2
done
