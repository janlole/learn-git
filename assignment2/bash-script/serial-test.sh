#!/bin/bash
# This script compile and run the serial test for select and the kdtree building function
# varying the number of dimensions received by the calling script


cd $PBS_O_WORKDIR
cd ..

saving_folder=data/serial

mkdir data
mkdir data/serial

module load gnu/9.3.0

g++ 14_select_test.cc -o 14_select_test.o -std=c++11 -O3 -DNDIM=${DIMENSION}

./14_select_test.o > ${saving_folder}/time-serial_${DIMENSION}.csv


