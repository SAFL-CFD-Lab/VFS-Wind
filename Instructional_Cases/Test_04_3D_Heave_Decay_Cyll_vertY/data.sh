#!/bin/bash
### Job name
#PS -N Charles_Test_Prime
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:10:00
#PBS -j oe

cd $PBS_O_WORKDIR

mpiexec ./data -tis 0 -tie 9700 -ts 100  > dataerr




