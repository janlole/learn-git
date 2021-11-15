#!/bin/bash

cd /u/dssc/lorenzo/try/learn-shell

rm hello.sh.*
rm -r data
mkdir data

# type of nodes involved
topol=('core' 'socket' 'node')
# how many iteration for measure
howmany=10
# saving folder
sav_fol=./data

running()
{
	# save different information in specific files
	local topo=$1	
	local folder=$2
	local num=$3
	mpirun --report-bindings --map-by $topo -np 2 ./IMB-MPI1 PingPong -msglog 28 2> error-tmp > results-tmp

	echo \#header_line 1: mpirun --report-bindings --map-by $topo -np 2 ./IMB-MPI1 PingPong > $folder/involving-$topo-$num.csv
	echo \#header_line 2: >> $folder/involving-$topo-$num.csv

	cat error-tmp | grep rank >> $folder/involving-$topo-$num.csv
	cat results-tmp | grep -v '^#' | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' > $folder/results-$topo-$num.csv
}



module load openmpi-4.1.1+gnu-9.3.0

for i in "${topol[@]}"
do
	foo=${sav_fol}/topo-$i
	mkdir $foo
	for (( j = 0 ; j < $howmany ; j++ ));
	do
		running $i $foo $j
	done
done

rm error-tmp results-tmp 
