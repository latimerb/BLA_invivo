#!/bin/bash

#SBATCH -p Lewis,RC
#SBATCH -N 4
#SBATCH -A rc
#SBATCH -n 80
#SBATCH --qos=normal
#SBATCH --job-name=ca1
#SBATCH --output=ca1%j.out
#SBATCH --time 0-00:30

module load intel/intel-2016-update2
module load nrn/nrn-mpi-7.4
module load openmpi/openmpi-2.0.0

module list
echo "Starting model at $(date)"

mpirun nrniv -mpi try_all_loadfromfile.hoc

echo "Simulation over at $(date)"



