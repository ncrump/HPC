#!/bin/bash

# Issue directives to scheduler:

#$ -N pyRun
#$ -cwd
#$ -o outputPy.txt
#$ -M ncrump@gmu.edu
#$ -j y
#$ -pe shared 16
#$ -l mf=1G
#$ -l m_arch=INTEL
#$ -q all.q

# Run commands needed here:

. /etc/profile
export PYTHONPATH=$PYTHONPATH:$PWD
module load python/python2.7
module load numpy/1.8.0

python calcCheck.py
