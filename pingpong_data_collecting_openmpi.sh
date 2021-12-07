#!/bin/bash

cd $PBS_O_WORKDIR

working_folder=data_col_mpi


rm -r ${working_folder}
mkdir ${working_folder}

# type of nodes involved
topol=('core' 'socket' 'node')
# array of pml protocol
list_pml=('ob1' 'ucx')
# array of btl protocol
list_btl=('vader' 'tcp' 'openib')
# array of network
list_net=('gig' 'infin')
# --mca btl_tcp_if_include br0
infin_check=infin


# how many iteration for measurement
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


	if [ ${net_l} == ${infin_check} ]
	then
		mpirun --report-bindings --mca pml ${pml_l} --mca btl ${btl_l},self --mca btl_tcp_if_include br0 --map-by $topo -np 2 ./IMB-MPI1 PingPong -msglog 28 2> error-tmp > results-tmp
		echo \#header_line 1: mpirun --report-bindings --mca pml ${pml_l} --mca btl ${btl_l},self --mca btl_tcp_if_include br0 --map-by $topo -np 2 ./IMB-MPI1 PingPong > $folder/involving-$num.csv
	else
		mpirun --report-bindings --mca pml ${pml_l} --mca btl ${btl_l},self --map-by $topo -np 2 ./IMB-MPI1 PingPong -msglog 28 2> error-tmp > results-tmp
		echo \#header_line 1: mpirun --report-bindings --mca pml ${pml_l} --mca btl ${btl_l},self --map-by $topo -np 2 ./IMB-MPI1 PingPong > $folder/involving-$num.csv
	fi
	
	echo \#header_line 2: >> $folder/involving-$num.csv
	cat error-tmp | grep rank >> $folder/involving-$num.csv

	cat results-tmp | grep -v '^#' | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' > $folder/results-$num.csv
}



module load openmpi-4.1.1+gnu-9.3.0

for i in "${topol[@]}"
do
	foo_top=${sav_fol}/topo-$i
	mkdir ${foo_top}
	for net_i in "${list_net[@]}"
	do
		foo_net=${foo_top}/net-${net_i}
		mkdir ${foo_net}
		for pml_i in "${list_pml[@]}"
		do
			foo_pml=${foo_net}/pml-${pml_i}
			mkdir ${foo_pml}
			for btl_i in "${list_btl[@]}"
			do
				if [[ ((${pml_i} == 'ob1') && (${btl_i} == 'openib')) || ((${i} == 'node') && (${pml_i} == 'ob1') && (${btl_i} == 'vader')) ]];
				then
					continue
				else
					foo_btl=${foo_pml}/btl-${btl_i}
					mkdir ${foo_btl}
					for (( j = 0 ; j < $howmany ; j++ ));
					do
						foo=${foo_btl}
						running $i $foo $j ${pml_i} ${btl_i} ${net_i}
					done
				fi
			done
		done
	done
done

rm error-tmp results-tmp 
