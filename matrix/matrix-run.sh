#!/bin/bash
declare -A matrix

matrix[0,0]=2400
matrix[0,1]=100
matrix[0,2]=100
matrix[1,0]=1200
matrix[1,1]=200
matrix[1,2]=100
matrix[2,0]=800
matrix[2,1]=300
matrix[2,2]=100

cd $PBS_O_WORKDIR

# loading openmpi module
module load openmpi-4.1.1+gnu-9.3.0
# compilation of the matrix.cc file
mpic++ matrix.cc -o matrix.o


working_folder=matrix-$DUNO-$DDUE-$DTRE

rm -r ${working_folder}
mkdir ${working_folder}


# how many iteration for measurement
howmany=10

for ((i = 0; i < 3; i++));
do
	# saving folder
	sav_fol=./${working_folder}/matrix-${i}
	mkdir ${sav_fol}
	for (( j = 0 ; j < $howmany ; j++ ));
	do
		mpirun -np 24 ./matrix.o $DUNO $DDUE $DTRE ${matrix[$i,0]} ${matrix[$i,1]} ${matrix[$i,2]} 2>>/dev/null >${sav_fol}/output-${j}.out
	done
done
