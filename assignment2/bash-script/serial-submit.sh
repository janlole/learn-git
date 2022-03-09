#!/bin/bash
# This script submit a job calling the serial-test.sh script
# for each number of dimension it will submit a job that
# 1. compile the code 14_select_test.cc with the right number of dimnesions ndim
# 2. run that code collecting the data in "time-serial_ndim.csv" 

dimension=( 2 3 4 5 )

for ((i = 0; i < 4; i++));
do
	ndim=${dimension[$i]}

	for (( prec = 0; prec < 2; prec++ ))
	do
		qsub -l nodes=1:ppn=2 -l walltime=4:00:00 -v PREC=${prec},DIMENSION=${ndim} -q dssc_gpu ./serial-test.sh
	done
done
