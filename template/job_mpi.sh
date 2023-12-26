#!/bin/bash -l
#########################################################
#SBATCH -J JOB_NAME_s50n512_simba
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o ./slurm-%j.out
#########################################################
#SBATCH --time=24:0:0
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=20
#########################################################

## Each Niagara node has 2 sockets (NUMA regions/chips) each with 20 cores --> so 2 tasks per node, and 20 cpus per task, and a total of np=2*8=16 processes
## This gives a total of 16*20=320 core processors for the simulation

#module load intel/2019u4 intelmpi/2019u4 hdf5/1.8.21 fftw/3.3.8 gsl/2.5 metis/5.1.0-shared parmetis/4.0.3-shared autotools tbb/2019u9
#module load intel/2019u4 intelmpi/2019u4 hdf5/1.8.21 fftw/3.3.8 gsl/2.5 metis/5.1.0-shared parmetis/4.0.3-shared
module load NiaEnv/.2022a intel/2022u2 intelmpi/2022u2+ucx-1.11.2 hdf5/1.10.9 fftw/3.3.10 gsl/2.7 parmetis/4.0.3-shared autotools

mpirun -np 4 ./swift_mpi_ps2020_sphenix_v2 --pin --cosmology --simba --threads=20 YML_FILE

#########################################################
