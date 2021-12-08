#!/bin/bash

D1_vec=24

qsub -Roe -l nodes=1:ppn=24 -l walltime=0:10:00 -q dssc -v DUNO=${D1_vec},DDUE=1,DTRE=1 ./matrix-run.sh

declare -A D2_vec

D2_vec[0,0]=12
D2_vec[0,1]=2
D2_vec[1,0]=2
D2_vec[1,1]=12
D2_vec[2,0]=6
D2_vec[2,1]=4
D2_vec[3,0]=4
D2_vec[3,1]=6
D2_vec[4,0]=3
D2_vec[4,1]=8
D2_vec[5,0]=8
D2_vec[5,1]=3

for ((i = 0; i < 6; i++));
do
	qsub -Roe -l nodes=1:ppn=24 -l walltime=0:10:00 -q dssc -v DUNO=${D2_vec[$i,0]},DDUE=${D2_vec[$i,1]},DTRE=1 ./matrix-run.sh
done

declare -A D3_vec

D3_vec[0,0]=2
D3_vec[0,1]=2
D3_vec[0,2]=6
D3_vec[1,0]=2
D3_vec[1,1]=6
D3_vec[1,2]=2
D3_vec[2,0]=6
D3_vec[2,1]=2
D3_vec[2,2]=2
D3_vec[3,0]=4
D3_vec[3,1]=2
D3_vec[3,2]=3
D3_vec[4,0]=4
D3_vec[4,1]=3
D3_vec[4,2]=2
D3_vec[5,0]=2
D3_vec[5,1]=4
D3_vec[5,2]=3
D3_vec[6,0]=2
D3_vec[6,1]=3
D3_vec[6,2]=4
D3_vec[7,0]=3
D3_vec[7,1]=4
D3_vec[7,2]=2
D3_vec[8,0]=3
D3_vec[8,1]=2
D3_vec[8,2]=4


for ((i = 0; i < 9; i++));
do
	qsub -Roe -l nodes=1:ppn=24 -l walltime=0:10:00 -q dssc -v DUNO=${D3_vec[$i,0]},DDUE=${D3_vec[$i,1]},DTRE=${D3_vec[$i,2]} ./matrix-run.sh
done
