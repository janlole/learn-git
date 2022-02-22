#!/bin/bash
# This script submit a job calling the mpi-test.sh script 

dimension=( 2 3 4 5 )

numberprocesses=( 2 4 8 16 32 )

numpoints=( 134217727 268435455 536870911 1073741823 2147483647 4294967295 )

for ((k = 0; k < 5; k++));
do
	num_proc=${numberprocesses[$k]}

	for ((i = 0; i < 1; i++));
	do
		n_dim=${dimension[$i]}

		for (( j = 0; j < 6; j++ ))
		do
			num_points=${numpoints[$j]}

			if [[ ${num_proc} -ge 8 ]] || [[ ${num_proc} -eq 2 && ${num_points} -lt 2147483647 ]] || [[ ${num_proc} -eq 4 && ${num_points} -le 2147483647 ]] 
			then
				echo ${num_proc} - ${num_points}
				# for (( prec = 0; prec < 2; prec++ ))
				# do
				# 	# qsub -l nodes=1:ppn=${num_proc} -l walltime=0:30:00 -v PREC=${prec},DIMENSION=${n_dim},NUMP=${num_points},NUMPROC=${num_proc} -q dssc_gpu ./mpi-test.sh
				# done
			# elif [[ ${num_proc} -eq 2 && ${num_points} -lt 2147483647 ]];
			# then
			# 	echo ${num_proc} - ${num_points}
			# 	# for (( prec = 0; prec < 2; prec++ ))
			# 	# do
			# 	# 	# qsub -l nodes=1:ppn=${num_proc} -l walltime=0:30:00 -v PREC=${prec},DIMENSION=${n_dim},NUMP=${num_points},NUMPROC=${num_proc} -q dssc_gpu ./mpi-test.sh
			# 	# done

			# elif [[ ${num_proc} -eq 4 && ${num_points} -le 2147483647 ]];
			# then
			# 	echo ${num_proc} - ${num_points}
				# for (( prec = 0; prec < 2; prec++ ))
				# do
				# 	# qsub -l nodes=1:ppn=${num_proc} -l walltime=0:30:00 -v PREC=${prec},DIMENSION=${n_dim},NUMP=${num_points},NUMPROC=${num_proc} -q dssc_gpu ./mpi-test.sh
				# done

			fi

		done
	done
done