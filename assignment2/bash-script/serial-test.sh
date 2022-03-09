#!/bin/bash
# This script compile and run the serial test for select and the kdtree building function
# varying the number of dimensions received by the calling script

# note: the bash scripts are in a subfolder called "bash-script" so it is neccessary to
# 		move in the parent directory
cd $PBS_O_WORKDIR
cd ..

saving_folder=data/serial/lastValues

mkdir data
mkdir data/serial
mkdir data/serial/lastValues

module load gnu/9.3.0

if [[ $PREC == 1 ]];
then
	g++ 14_select_test.cc -o 14_select_test.o -std=c++11 -O3 -DNDIM=${DIMENSION} -DDOUBLE_PRECISION_KPOINT
else
	g++ 14_select_test.cc -o 14_select_test.o -std=c++11 -O3 -DNDIM=${DIMENSION}
fi


./14_select_test.o > ${saving_folder}/time-serial-double_${DIMENSION}-${PREC}.csv


