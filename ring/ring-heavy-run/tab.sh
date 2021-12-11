#!/bin/bash

for (( i = 1; i <= $1; i++ ));
do
	for (( j = 1 ; j <= $2 ; j++ ));
	do
		num=$(($i*$j))
		echo "numproc,exec-time,init-time" > ring-node-$i-core-$num/summary.csv
		for (( k = 0 ; k < 10 ; k++));
		do
			cat ring-node-$i-core-$num/output-${k}.out | tail -1 | tr -s '\t ' ',' >> ring-node-$i-core-$num/summary.csv
		done
	done
done
