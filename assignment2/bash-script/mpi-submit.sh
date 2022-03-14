#!/bin/bash
# This script submit a job calling the mpi-test.sh script 

# WARNINGS: present configuration to collect information for 2^31  kpoints but just for numberprocesses = 32

dimension=( 2 3 4 5 )

numberprocesses=( 2 4 8 16 32 )

numpoints=( 134217727 268435455 536870911 1073741823 2147483647 4294967295 )

for ((k = 4; k < 5; k++)); # numberprocesses
do
	num_proc=${numberprocesses[$k]}

	for ((i = 0; i < 4; i++)); # dimension
	do
		n_dim=${dimension[$i]}

		for (( j = 4; j < 5; j++ )) # numpoints
		do
			num_points=${numpoints[$j]}
			for (( prec = 0; prec < 2; prec++ )) # prec
			do
				qsub -l nodes=1:ppn=${num_proc} -l walltime=4:00:00 -v PREC=${prec},DIMENSION=${n_dim},NUMP=${num_points},NUMPROC=${num_proc} -q dssc_gpu ./mpi-test.sh
			done
		done
	done
done
