#!/bin/bash
### Job name
#PS -N Charles_Test_Prime
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=1:ppn=16
#PBS -l walltime=5:00:00
#PBS -j oe

cd /$PBS_O_WORKDIR

/safl/software/aegean/openmpi/1.5.5/gcc/4.7.0-fix/bin/mpirun --bind-to-core ./vwis > err1
