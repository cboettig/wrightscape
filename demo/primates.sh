#!/bin/bash
## Use the current working directory
#$ -cwd
## use bash commands
#$ -S /bin/bash
## Launch parallel mpi threads 
#$ -pe mpi 161
#$ -o primates.sungrid.out
# combine error and output files
#$ -j y


module load gcc openmpi R Rmpi
setenv R_HOME ~/R/x86_64-redhat-linux-gnu-library/2.13
R -f primates.R
~                                
