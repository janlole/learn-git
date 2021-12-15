#!/bin/bash

cd $PBS_O_WORKDIR

# loading openmpi module
module load openmpi-4.1.1+gnu-9.3.0

sav_fol=jacobi-1-proc
mkdir ${sav_fol}

# how many iteration for measurement
howmany=10

for (( j = 0 ; j < $howmany ; j++ ));
do
	out=${sav_fol}/output-${j}.out
	timing=${sav_fol}/timing-${j}.txt

	(time mpirun -np 1 ./jacoby3D.x <input/input.1.0 2>/dev/null > ${out} ) 2> ${timing}
done