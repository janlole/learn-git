#!/bin/bash


qsub -l nodes=2:ppn=3 -l walltime=0:00:05 -q dssc hello.sh
