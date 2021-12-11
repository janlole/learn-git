#!/bin/bash

topol=('core' 'socket')
decom=( 6 10 18 )

for ((i = 1; i <= 3; i++));
do
	num_proc=$((${i}*4))
	deco_for=${decom[$i]}

	for ((k = 0; k < ${deco_for}; k++));
	do
		for j in "${topol[@]}"
		do
			qsub -l nodes=1:ppn=24 -l walltime=0:10:00 -v CORE=${num_proc},TOPO=${j},DEC=${k} -q dssc ./jacobi-one.sh
		done
	done
done