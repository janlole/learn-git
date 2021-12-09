#!/bin/bash

cd $PBS_O_WORKDIR

# loading openmpi module
module load openmpi-4.1.1+gnu-9.3.0
# compilation of the ring.cc file
mpic++ ring.cc -o ring.o

num=$(($NODE*$CORE))

working_folder=ring-node-$NODE-core-${num}

rm -r ${working_folder}
mkdir ${working_folder}


# how many iteration for measurement
howmany=10
# saving folder
sav_fol=./${working_folder}



for (( j = 0 ; j < $howmany ; j++ ));
do
	mpirun -np ${num} --map-by node ./ring.o > ${sav_fol}/output-${j}.out
done

