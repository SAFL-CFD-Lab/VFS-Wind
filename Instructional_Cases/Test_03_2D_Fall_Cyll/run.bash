#!/bin/bash

#SBATCH --nodes=2               # Number of nodes
#SBATCH --time=5:59:00            # Wall clock time
#SBATCH --account=FY140378         # WC ID
#SBATCH --job-name=MonochromaticWave     # Name of job
#SBATCH --partition=ec            # partition name

nodes=$SLURM_JOB_NUM_NODES          # Number of nodes
cores=8                            # Number MPI processes to run on each node (a.k.a. PPN)

cd /gscratch1/ancalde/Validation_Cases/03_it100tau001

module purge
module load gnu/4.4.6
module load openmpi-gnu/1.6
mpiexec --bind-to-core --npernode $cores --n $(($cores*$nodes)) ./test_v1.00>&err
