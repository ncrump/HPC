#!/bin/bash

# Issue directives to scheduler:

#$ -N ompRun
#$ -cwd
#$ -o output.txt
#$ -M ncrump@gmu.edu
#$ -j y
#$ -pe shared 16
#$ -l mf=1G
#$ -q all.q,all-HiPri.q,all-LoPri.q

# Run commands needed here:
. /etc/profile
OMP_NUM_THREADS=16
time ./2dvectorfield
