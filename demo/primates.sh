#!/bin/bash
## Use the current working directory
#$ -cwd
## use bash commands
#$ -S /bin/bash
## Launch parallel mpi threads 
#$ -pe mpi 64
#$ -o primates2.out
# combine error and output files
#$ -j y
#$ -N primates2

module load gcc openmpi R Rmpi
R -f primates2.R
~                                
