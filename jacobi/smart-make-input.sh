#!/bin/bash


declare -A P4

P4[0,0]=4 P4[0,1]=1 P4[0,2]=1
P4[1,0]=1 P4[1,1]=4 P4[1,2]=1
P4[2,0]=1 P4[2,1]=1 P4[2,2]=4
P4[3,0]=2 P4[3,1]=2 P4[3,2]=1
P4[4,0]=1 P4[4,1]=2 P4[4,2]=2
P4[5,0]=2 P4[5,1]=1 P4[5,2]=2

check_ex()
{
	declare -n array=$1
	num_triplet=$2
	elem_0=$3
	elem_1=$4
	elem_2=$5
	check=1
	if [[ ${num_triplet} -eq 0 ]];
	then
		return
	fi
	for (( j = 0 ; j < ${num_triplet}; j++));
	do
		if  [[ ${elem_0} -eq ${array[${j},0]} &&  ${elem_1} -eq ${array[${j},1]} && ${elem_2} -eq ${array[${j},2]} ]]
		then
			check=0
			return
		fi
	done
}

new_list()
{
	# WARINING: use only with a prime moltiplicator
	declare -n old_list="$4"
	multiplicator=$2
	num_triplet_old=$3
	declare -n new_array="$1"

	counter=0
	for (( i = 0 ; i < ${num_triplet_old} ; i++ ));
	do
		place=$(( $i * 3 ))

		new_elem=$(( ${old_list[$i,0]} * $multiplicator))
		check_ex new_array ${counter} ${new_elem} ${old_list[$i,1]} ${old_list[$i,2]}

		if [[ $check -eq 1 ]];
		then
			new_array[$counter,0]=${new_elem}
			new_array[$counter,1]=${old_list[$i,1]}
			new_array[$counter,2]=${old_list[$i,2]}
			counter=$(( ${counter} + 1 ))
		fi

		new_elem=$(( ${old_list[$i,1]} * $multiplicator))
		check_ex new_array ${counter} ${old_list[$i,0]} ${new_elem} ${old_list[$i,2]}
		if [[ $check -eq 1 ]];
		then
			new_array[$counter,0]=${old_list[$i,0]}
			new_array[$counter,1]=${new_elem}
			new_array[$counter,2]=${old_list[$i,2]}
			counter=$(( ${counter} + 1 ))
		fi

		new_elem=$(( ${old_list[$i,2]} * $multiplicator))
		check_ex new_array ${counter} ${old_list[$i,0]} ${old_list[$i,1]} ${new_elem}
		if [[ $check -eq 1 ]];
		then
			new_array[$counter,0]=${old_list[$i,0]}
			new_array[$counter,1]=${old_list[$i,1]}
			new_array[$counter,2]=${new_elem}
			counter=$(( ${counter} + 1 ))
		fi
		
	done

}

input_maker()
{
	local -n array=$1
	local saving_folder=$2
	local num_proce=$3
	local num_side=$4
	local num_decompo=$5

	for ((i = 0; i < ${num_decompo}; i++));
	do
		echo $((${array[$i,0]}*${num_side}))","${array[$i,0]}",T" >  ${saving_folder}/input.${num_proce}.${i}
		echo $((${array[$i,1]}*${num_side}))","${array[$i,1]}",T" >> ${saving_folder}/input.${num_proce}.${i}
		echo $((${array[$i,2]}*${num_side}))","${array[$i,2]}",T" >> ${saving_folder}/input.${num_proce}.${i}
		echo >> ${saving_folder}/input.${num_proce}.${i}
	done
}

# 8 processes, 10 decompositions
declare -A PN8 
# 12 processes, 18 decompositions
declare -A PN12
# 24 processes, 30 decompositions
declare -A PN24
# 48 processes, 45 decompositions
declare -A PN48

new_list PN8 2 6 P4
new_list PN12 3 6 P4
new_list PN24 2 18 PN12
new_list PN48 2 30 PN24

folder=input

mkdir ${folder}

side=600

input_maker P4   ${folder} 4  ${side} 6
input_maker PN8  ${folder} 8  ${side} 10
input_maker PN12 ${folder} 12 ${side} 18
input_maker PN24 ${folder} 24 ${side} 30
input_maker PN48 ${folder} 48 ${side} 45
