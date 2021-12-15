#!/bin/bash


qsub -Roe -l nodes=2:ppn=24 -l walltime=1:00:00 -q dssc pingpong_noflag.sh




