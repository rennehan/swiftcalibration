#!/bin/bash -l
#########################################################
#SBATCH -J JOB_NAME
#SBATCH -p cca
#SBATCH --mail-user=drennehan@flatironinstitute.org
#SBATCH --mail-type=FAIL 
#SBATCH -C rome
#SBATCH -o ./slurm-%j.out
#########################################################
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#########################################################

module purge
module load modules/2.3-20240529
module load openblas/single-0.3.26
module load gcc/11.4.0
module load openmpi/4.0.7
module load hdf5/1.12.3
module load gsl/2.7.1
module load fftw/mpi-3.3.10

./swift -r --pin --cosmology --simba --threads=${SLURM_CPUS_PER_TASK} YML_FILE

#########################################################
