#!/bin/bash

qsub -Roe -l nodes=$1:ppn=24 -l walltime=0:10:00 -q dssc -v NODE=$1,CORE=$2 ./ring.sh
