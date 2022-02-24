#!/bin/bash
# This script submit a job calling the mpi-test.sh script 
# the if-else block is needed to NOT submit jobs that requires more than 30 minutes to complete

dimension=( 2 3 4 5 )

numberprocesses=( 2 4 8 16 32)

numpoints=( 134217727 268435455 536870911 1073741823 )

for ((k = 0; k < 5; k++));
do
	num_proc=${numberprocesses[$k]}

	for ((i = 0; i < 4; i++));
	do
		n_dim=${dimension[$i]}

		for (( j = 0; j < 4; j++ ))
		do
			num_points=${numpoints[$j]}

			if [[ ${num_proc} -ge 8 ]] || [[ ${num_proc} -eq 2 && ${num_points} -lt 2147483647 ]] || [[ ${num_proc} -eq 4 && ${num_points} -le 2147483647 ]] 
			then
				for (( prec = 0; prec < 2; prec++ ))
				do
					qsub -l nodes=1:ppn=${num_proc} -l walltime=1:00:00 -v PREC=${prec},DIMENSION=${n_dim},NUMP=${num_points},NUMPROC=${num_proc} -q dssc_gpu ./omp-test.sh
				done
			fi
		done
	done
done
