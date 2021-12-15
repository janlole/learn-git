#!/bin/bash

cd $PBS_O_WORKDIR

# loading openmpi module
module load openmpi-4.1.1+gnu-9.3.0

working_folder=jacobi-${CORE}-${TOPO}
mkdir ${working_folder}

# how many iteration for measurement
howmany=10

# saving folder
sav_fol=${working_folder}/deco-${DECO}
mkdir ${sav_fol}

for (( j = 0 ; j < $howmany ; j++ ));
do
	out=${sav_fol}/output-${j}.out
	timing=${sav_fol}/timing-${j}.txt

	(time mpirun -np ${CORE} --map-by ${TOPO} ./jacoby3D.x <input/input.${CORE}.${DECO} 2>/dev/null > ${out} ) 2> ${timing}
done