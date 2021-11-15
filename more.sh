#!/bin/bash

cd /u/dssc/lorenzo/try/learn-shell

working_folder=data_1

rm hello.sh.*
rm -r ${working_folder}
mkdir ${working_folder}

# type of nodes involved
topol=('core' 'socket' 'node')
# array of pml protocol
list_pml=('ob1' 'cm' 'ucx')
# array of btl protocol
list_btl=('vader' 'tcp' 'uct' 'openib')
# array of network
list_net=('' '--mca btl_tcp_if_include br0')

# how many iteration for measure
howmany=1
# saving folder
sav_fol=./${working_folder}

running()
{
	# save different information in specific files
	local topo=$1	
	local folder=$2
	local num=$3
	local pml_l=$4
	local btl_l=$5
	local net_l=$6

	mpirun --report-bindings --mca pml ${pml_l} --mca btl ${btl_l},self ${net_l} --map-by $topo -np 2 ./IMB-MPI1 PingPong -msglog 28 2> error-tmp > results-tmp

	echo \#header_line 1: mpirun --report-bindings --mca pml ${pml_l} --mca btl ${btl_l},self ${net_l} --map-by $topo -np 2 ./IMB-MPI1 PingPong > $folder/involving-$topo-$num.csv
	echo \#header_line 2: >> $folder/involving-$topo-$num.csv

	# cat error-tmp | grep rank >> $folder/involving-$topo-$num.csv
	# following lines to understand which combinations works
	echo "---------------------" >> summary.csv
	echo "---------------------" >> summary.csv
	echo "combination: " >> summary.csv
	echo "top:		${topo}" >> summary.csv
	echo "pml:		${pml_l}" >> summary.csv
	echo "btl:		${btl_l}" >> summary.csv
	echo "net:		${net_l}" >> summary.csv
	cat error-tmp >> summary.csv
	echo "---------------------" >> summary.csv
	echo "---------------------" >> summary.csv

	cat results-tmp | grep -v '^#' | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' > $folder/results-$topo-$num.csv
}



module load openmpi-4.1.1+gnu-9.3.0

for i in "${topol[@]}"
do
	foo=${sav_fol}/topo-$i
	mkdir $foo
	for net_i in "${list_net[@]}"
	do
		foo=${sav_fol}/net-${net_i}
		mkdir $foo
		for pml_i in "${list_pml[@]}"
		do
			foo=${foo}/pml-${pml_i}
			mkdir $foo
			for btl_i in "${list_btl[@]}"
			do
				foo=${foo}/btl-${btl_i}
				mkdir foo
				for (( j = 0 ; j < $howmany ; j++ ));
				do
					running $i $foo $j ${pml_i} ${btl_i} ${net_i}
				done
			done
		done
	done
done

rm error-tmp results-tmp 
