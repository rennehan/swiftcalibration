#!/bin/bash -l
#########################################################
#SBATCH -J JOB_NAME
#SBATCH --mail-user=douglas.rennehan@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/rennehan/swift_tests/cali_simba/JOB_NAME/slurm-%j.out
#########################################################
#SBATCH --time=24:0:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#########################################################

module load NiaEnv/.2022a intel/2022u2 intelmpi/2022u2+ucx-1.11.2 hdf5/1.10.9 fftw/3.3.10 gsl/2.7 parmetis/4.0.3-shared
./swift --pin --cosmology --simba --threads=40 YML_FILE

#########################################################
