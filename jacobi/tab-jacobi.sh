#!/bin/bash
#######################################
# make an .csv file in each deco-* folder
# this .csv file contains the time of all the runs for that decomposition
#######################################

working_folder=test-data

topol=( 'core' 'socket' )
num_proc=( 4 8 12 )
num_deco=( 6 10 18 )
howmany=10

cd ${working_folder}

for (( i = 0; i < 3; i++ ));
do
	proc=${num_proc[$i]}
	deco=${num_deco[$i]}
	for top in "${topol[@]}"
	do
		for (( j = 0; j < ${deco} ; j++ ));
		do
			subfolder="jacobi-${proc}-${top}/deco-${j}"
			echo "real,user,sys" > ${subfolder}/summary.csv
			
			for (( k = 0; k < ${howmany} ; k++ ));
			do
				cat ${subfolder}/timing-${k}.txt | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' > tmp 
				echo '' > tmp1
				while read line1 ;
				do
					read line2
					echo ${line1}*60 +  ${line2} | bc >> tmp1
				done < tmp
				cat tmp1 | tr -s '\n' ',' | sed 's/,//' | sed 's/.\{1\}$//' >> ${subfolder}/summary.csv
				echo '' >> ${subfolder}/summary.csv
			done
		done
	done
done
