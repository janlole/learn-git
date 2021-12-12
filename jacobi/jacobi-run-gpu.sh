#!/bin/bash
# This script submit a job calling the jacobi-gpu.sh script 

decom=( 18 30 45 )
number_processes=( 12 24 48 )

for ((i = 1; i <= 3; i++));
do
	num_proc=${number_processes[$i]}
	deco_for=${decom[$i]}

	for ((k = 0; k < ${deco_for}; k++));
	do
		qsub -l nodes=1:ppn=48 -l walltime=0:10:00 -v CORE=${num_proc},DEC=${k} -q dssc_gpu ./jacobi-gpu.sh

	done
done