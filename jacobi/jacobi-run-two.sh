#!/bin/bash
# This script submit a job calling the jacobi-one.sh script 

topol=node
decom=( 18 30 45 )
number_processes=( 12 24 48 )

for ((i = 1; i <= 3; i++));
do
	num_proc=${number_processes[$i]}
	deco_for=${decom[$i]}

	for ((k = 0; k < ${deco_for}; k++));
	do
		qsub -l nodes=2:ppn=24 -l walltime=0:10:00 -v CORE=${num_proc},TOPO=${topol},DEC=${k} -q dssc ./jacobi-one.sh

	done
done