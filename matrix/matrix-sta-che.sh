#!/bin/bash

qsub -Roe -l nodes=1:ppn=24 -l walltime=0:20:00 -q dssc ./matrix-che.sh

