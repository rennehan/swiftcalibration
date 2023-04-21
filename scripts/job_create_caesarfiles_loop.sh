#!/bin/bash -l
#########################################################
#SBATCH -J create_caesarfiles
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/slurm-%j.out
#########################################################
#SBATCH --time=3:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar

for i in 0{0..9}
do
	python create_caesarfiles.py /scratch/b/babul/aspadawe/swift_tests/cali_simba/ cali_00${i}
done

for i in {10..15}
do
        python create_caesarfiles.py /scratch/b/babul/aspadawe/swift_tests/cali_simba/ cali_00${i}
done

#########################################################
