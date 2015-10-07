#!/bin/bash

# Issue directives to scheduler:

#$ -N mpiRun
#$ -cwd
#$ -o output.txt
#$ -M ncrump@gmu.edu
#$ -j y
#$ -pe openmpi 4
#$ -v MPI_HOME
#$ -l mf=1G
#$ -q all.q,all-HiPri.q,all-LoPri.q

# Run commands needed here:

. /etc/profile
module load openmpi/gcc/64/1.8.1
CMD="$MPI_HOME/bin/mpirun -np 4 ./2dParticles_mpi"
time $CMD
