#!/bin/bash

#SBATCH --nodes=1               # Number of nodes
#SBATCH --time=0:59:00            # Wall clock time
#SBATCH --account=FY140378         # WC ID
#SBATCH --job-name=MonochromaticWave     # Name of job
#SBATCH --partition=ec            # partition name

nodes=$SLURM_JOB_NUM_NODES          # Number of nodes
cores=8                            # Number MPI processes to run on each node (a.k.a. PPN)

cd /gscratch1/ancalde/code_validation_test/Test_01_2D_Monochromatic_wave

module purge
module load gnu/4.4.6
module load openmpi-gnu/1.6
mpiexec --bind-to-core --npernode $cores --n $(($cores*$nodes)) ./test_v0.00>&err
