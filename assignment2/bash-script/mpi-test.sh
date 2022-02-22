#!/bin/bash
# This script compile and run 11_parallel_mpi-01.cc according to these variables
# NUMPROC = number of processor executing the code
# DIMENSION = number of dimensions
# NUMP = number of kpoints
# PREC = to make all in double precision

cd $PBS_O_WORKDIR
cd ..

# loading openmpi module
module load openmpi-4.1.1+gnu-9.3.0

saving_folder=data/parallel-mpi

mkdir data
mkdir data/parallel-mpi

if [[ $PREC == 1 ]];
then
	mpic++ 11_parallel_mpi-01.cc -o 11_parallel_mpi-01.o -O3 -DNUMPOINTS=${NUMP} -DDOUBLE_PRECISION_KPOINT -DNDIM=${DIMENSION}
else
	mpic++ 11_parallel_mpi-01.cc -o 11_parallel_mpi-01.o -O3 -DNUMPOINTS=${NUMP}  -DNDIM=${DIMENSION}
fi

mpirun -np ${NUMPROC} ./11_parallel_mpi-01.o > ${saving_folder}/time_${NUMPROC}-${NUMP}-${DIMENSION}-${PREC}.csv
