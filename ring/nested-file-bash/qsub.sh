#!/bin/bash

qsub -Roe -l nodes=$1:ppn=$2 -l walltime=0:10:00 -q dssc -v NODE=$1,CORE=$2 ./ring.sh
