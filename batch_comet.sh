#!/bin/bash

#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A TG-DBS180005
#SBATCH --job-name=ca1
#SBATCH --output=ca1%j.out
#SBATCH --time 0-00:10
#SBATCH --qos=normal


echo "Starting model at $(date)"

mpirun nrniv -mpi main.hoc
