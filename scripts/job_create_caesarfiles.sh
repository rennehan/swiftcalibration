#!/bin/bash -l
#########################################################
#SBATCH -J create_caesarfiles
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/s18n128_simba/slurm-%j.out
#########################################################
#SBATCH --time=4:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar

for FILE in /scratch/b/babul/aspadawe/swift_tests/cali_simba/cali_*; do
	python create_caesarfiles.py /scratch/b/babul/aspadawe/swift_tests/cali_simba/ $(basename $FILE)
done

#########################################################
