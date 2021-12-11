#!/bin/bash

max_time=$(($1 * 10))

qsub -Roe -l nodes=$1:ppn=24 -l walltime=0:${max_time}:00 -q dssc -v NODE=$1,CORE=$2 ./ring.sh
