#!/bin/bash

#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=testmpi%j.out

module load nrn/nrn-mpi-7.4
module load openmpi/openmpi-2.0.0

mpirun nrniv -mpi testmpi.hoc
