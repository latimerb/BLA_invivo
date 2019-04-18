#!/bin/bash

#SBATCH --partition compute
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=20
#SBATCH -A TG-DBS180005
#SBATCH --job-name=ca1
#SBATCH --output=ca1%j.out
#SBATCH --time 0-01:00
#SBATCH --qos=normal


echo "Starting model at $(date)"

mpirun nrniv -mpi main.hoc
