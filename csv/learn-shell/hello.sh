#!/bin/bash
# questo è un commento

cd /u/dssc/lorenzo/try/learn-shell
rm resu*
rm hello.sh.*

module load openmpi-4.1.1+gnu-9.3.0

mpirun --report-bindings --map-by core -np 2 -npernode 2  ./IMB-MPI1 PingPong 2>error > results-tmp
cat error | grep rank > results-0.csv
cat results-tmp | grep -v  ^# | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' >> results-0.csv
#cat results-tmp | grep -v ^# | grep -v '^$' >> results-0.csv
cat results-tmp | grep -v  ^# | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' >> results.csv

mpirun --report-bindings --map-by socket -np 2 -npernode 2  ./IMB-MPI1 PingPong 2>error > results-tmp
cat error | grep rank > results-1.csv
cat results-tmp | grep -v  ^# | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' >> results-1.csv
#cat results-tmp | grep -v ^# | grep -v '^$' >> results-1.csv

mpirun --report-bindings --map-by node -np 2 -npernode 1  ./IMB-MPI1 PingPong 2>error > results-tmp
cat error | grep rank > results-2.csv
cat results-tmp | grep -v  ^# | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' >> results-2.csv
#cat results-tmp | grep -v ^# | grep -v '^$' >> results-2.csv



#mpirun --report-bindings --map-by core -np 2 -npernode 2  ./IMB-MPI1 PingPong 2>/dev/null | grep -v  ^# | grep -v '^$' | tr -s '\t ' ',' | sed 's/,//' > results.txt

