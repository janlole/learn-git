#!/bin/bash

folder=input

mkdir ${folder}

side=600

declare -A P4

P4[0,1]=4 P4[0,0]=1 P4[0,2]=1
P4[1,0]=1 P4[1,1]=4 P4[1,2]=1
P4[2,0]=1 P4[2,1]=1 P4[2,2]=4
P4[3,0]=2 P4[3,1]=2 P4[3,2]=1
P4[4,0]=1 P4[4,1]=2 P4[4,2]=2
P4[5,0]=2 P4[5,1]=1 P4[5,2]=2

declare -A P8

P8[0,0]=2 P8[0,1]=2 P8[0,2]=2
P8[1,0]=4 P8[1,1]=2 P8[1,2]=1
P8[2,0]=4 P8[2,1]=1 P8[2,2]=2
P8[3,0]=2 P8[3,1]=4 P8[3,2]=1
P8[4,0]=2 P8[4,1]=1 P8[4,2]=4
P8[5,0]=1 P8[5,1]=4 P8[5,2]=2
P8[6,0]=1 P8[6,1]=2 P8[6,2]=4
P8[7,0]=1 P8[7,1]=1 P8[7,2]=8
P8[8,0]=1 P8[8,1]=8 P8[8,2]=1
P8[9,0]=8 P8[9,1]=1 P8[9,2]=1

declare -A P12

P12[0,0]=3 P12[0,1]=4 P12[0,2]=1		P12[6,0]=2 P12[6,1]=6 P12[6,2]=1
P12[1,0]=3 P12[1,1]=1 P12[1,2]=4		P12[7,0]=2 P12[7,1]=1 P12[7,2]=6
P12[2,0]=4 P12[2,1]=1 P12[2,2]=3		P12[8,0]=6 P12[8,1]=1 P12[8,2]=2
P12[3,0]=4 P12[3,1]=3 P12[3,2]=1		P12[9,0]=6 P12[9,1]=2 P12[9,2]=1
P12[4,0]=1 P12[4,1]=3 P12[4,2]=4		P12[10,1]=1 P12[10,0]=2 P12[10,2]=6
P12[5,0]=1 P12[5,1]=4 P12[5,2]=3		P12[11,0]=1 P12[11,1]=6 P12[11,2]=2

P12[12,0]=3 P12[12,1]=2 P12[12,2]=2		P12[15,0]=1 P12[15,1]=1 P12[15,2]=12
P12[13,0]=2 P12[13,1]=3 P12[13,2]=2		P12[16,0]=1 P12[16,1]=12 P12[16,2]=1
P12[14,0]=2 P12[14,1]=2 P12[14,2]=3		P12[17,0]=12 P12[17,1]=1 P12[17,2]=1


for ((i = 0; i < 6; i++));
do
	echo $((${P4[$i,0]}*${side}))","${P4[$i,0]}",T" >  ${folder}/input.4.${i}
	echo $((${P4[$i,1]}*${side}))","${P4[$i,1]}",T" >> ${folder}/input.4.${i}
	echo $((${P4[$i,2]}*${side}))","${P4[$i,2]}",T" >> ${folder}/input.4.${i}
	echo >> ${folder}/input.4.${i}
done

for ((i = 0; i < 10; i++));
do
	echo $((${P8[$i,0]}*${side}))","${P8[$i,0]}",T" >  ${folder}/input.8.${i}
	echo $((${P8[$i,1]}*${side}))","${P8[$i,1]}",T" >> ${folder}/input.8.${i}
	echo $((${P8[$i,2]}*${side}))","${P8[$i,2]}",T" >> ${folder}/input.8.${i}
	echo >> ${folder}/input.8.${i}
done

for ((i = 0; i < 18; i++));
do
	echo $((${P12[$i,0]}*${side}))","${P12[$i,0]}",T" >  ${folder}/input.12.${i}
	echo $((${P12[$i,1]}*${side}))","${P12[$i,1]}",T" >> ${folder}/input.12.${i}
	echo $((${P12[$i,2]}*${side}))","${P12[$i,2]}",T" >> ${folder}/input.12.${i}
	echo >> ${folder}/input.12.${i}
done


	