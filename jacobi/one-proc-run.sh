#!/bin/bash

qsub -Roe -l nodes=1:ppn=24 -l walltime=1:00:00 -q dssc ./one-proc-jacobi.sh