#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pi mpi 11
module load gcc openmpi R Rmpi
R -f primates.R
~                                
