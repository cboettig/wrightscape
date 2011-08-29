#!/bin/bash
## Use the current working directory
#$ -cwd
## use bash commands
#$ -S /bin/bash
## Launch parallel mpi threads 
#$ -pe mpi 16
#$ -o primates_faster.out
# combine error and output files
#$ -j y


module load gcc openmpi R Rmpi
R -f primates.R
~                                
