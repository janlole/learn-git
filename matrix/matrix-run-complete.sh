#!/bin/bash 
D1_vec=24 

qsub -Roe -l nodes=1:ppn=24 -l walltime=0:10:00 -q dssc -v DUNO=1,DDUE=${D1_vec},DTRE=1 ./matrix-run.sh 
qsub -Roe -l nodes=1:ppn=24 -l walltime=0:10:00 -q dssc -v DUNO=1,DDUE=1,DTRE=${D1_vec} ./matrix-run.sh 
 
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
	qsub -Roe -l nodes=1:ppn=24 -l walltime=0:10:00 -q dssc -v DUNO=${D2_vec[$i,0]},DDUE=1,DTRE=${D2_vec[$i,1]} ./matrix-run.sh 
	qsub -Roe -l nodes=1:ppn=24 -l walltime=0:10:00 -q dssc -v DUNO=1,DDUE=${D2_vec[$i,1]},DTRE=${D2_vec[$i,0]} ./matrix-run.sh 
done
         
