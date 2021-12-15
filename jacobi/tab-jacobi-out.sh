#!/bin/bash
######################################
# Extract from the "output-*.out" files the last ten rows with all the times 
# given by Jacoby-3D.x and make theme readable as .csv file
#
# Variable to set
#	- working_folder 
# 	- howmany runs were completed
######################################
scan_sumarize ()
{
for top in ${topol[@]};
do
	for (( i = 0; i < 3; i++ ));
	do
		proc=${proc_num[$i]}
		deco=${deco_num[$i]}
		for (( j = 0 ; j < ${deco} ; j++ ));
		do
			fold=${working_folder}/jacobi-${proc}-${top}/deco-${j}
			for (( k =  0; k < ${howmany} ; k++));
			do
				echo '1,2,3,4,5,6,7,8,9,10,11,12,13,14' > ${fold}/output-time-${k}.csv
				cat ${fold}/output-${k}.out | tail -10 | tr -s ' ' ',' | sed 's/.\{1\}$//' |  sed 's/,//'>> ${fold}/output-time-${k}.csv
			done
		done
	done
done
}

working_folder=test-data
howmany=10

topol=( 'core' 'socket' )
proc_num=( 4 8 12 )
deco_num=( 6 10 18 )

scan_sumarize

topol=( 'node' )
proc_num=( 12 24 48 )
deco_num=( 18 30 45 )

scan_sumarize

topol=( 'gpu' )
proc_num=( 12 24 )
deco_num=( 18 30 )

scan_sumarize

topol=( 'gpu' )
proc_num=( 48 )
deco_num=( 45 )
# Due to errors in the computation of the expected execution time
# the case "gpu with 48 cores" was NOT able to complete all the 10 runs.
# Despite the incompleted submit, I will not collect new data because 
# I consider 7 runs enough
howmany=7

scan_sumarize


# NOTE: ALL the runs are executed in a safe "communication-space"
# 		in the sense that for every test the nodes requested were 
# 		busy with just Jacoby program to execute