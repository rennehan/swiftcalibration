#!/bin/bash -l
#########################################################
#SBATCH -J JOB_NAME_s25n256_simba
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o ./slurm-%j.out
#########################################################
#SBATCH --time=24:0:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#########################################################

#module load intel/2019u4 intelmpi/2019u4 hdf5/1.8.21 fftw/3.3.8 gsl/2.5 metis/5.1.0-shared parmetis/4.0.3-shared
module load NiaEnv/.2022a intel/2022u2 intelmpi/2022u2+ucx-1.11.2 hdf5/1.10.9 fftw/3.3.10 gsl/2.7 parmetis/4.0.3-shared autotools

./swift --pin --cosmology --simba --threads=40 YML_FILE
#./swift --pin --cosmology --simba --threads=40 swimba_s12.5n128.yml

#########################################################
