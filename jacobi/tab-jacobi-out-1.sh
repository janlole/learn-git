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

	fold=${working_folder}/jacobi-1-proc
	for (( k =  0; k < ${howmany} ; k++));
	do
		echo '1,2,3,4,5,6,7,8,9,10,11,12,13,14' > ${fold}/output-time-${k}.csv
		cat ${fold}/output-${k}.out | tail -10 | tr -s ' ' ',' | sed 's/.\{1\}$//' |  sed 's/,//'>> ${fold}/output-time-${k}.csv
	done

}

working_folder=.
howmany=10

scan_sumarize


# NOTE: ALL the runs are executed in a safe "communication-space"
# 		in the sense that for every test the requested nodes were 
# 		busy with just Jacoby program to execute