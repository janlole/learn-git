#!/bin/bash
folder=./matrix-check
for(( j = 0; j < 10; j++));
do
	subfolder=${folder}/matrix-$j
	echo "pre.send,main.matrix,total" > ${subfolder}/summary.csv
	for (( k = 0 ; k < 10 ; k++));
	do
		cat ${subfolder}/output-${k}.out | tail -1 | tr -s '\t ' ',' >>  ${subfolder}/summary.csv
	done
done