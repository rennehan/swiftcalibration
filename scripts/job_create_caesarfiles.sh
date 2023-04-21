#!/bin/bash -l
#########################################################
#SBATCH -J create_caesarfiles
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/slurm-%j.out
#########################################################
#SBATCH --time=2:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar
#python create_caesarfiles.py /scratch/b/babul/aspadawe/swift_tests/cali_simba/ cali_0000 cali_0001 cali_0002 cali_0003 cali_0004 cali_0005 cali_0006 cali_0007 cali_0008 cali_0009 cali_0010 cali_0011 cali_0012 cali_0013 cali_0014 cali_0015
python create_caesarfiles.py /scratch/b/babul/aspadawe/swift_tests/cali_simba/ cali_0013 cali_0014 cali_0015

#########################################################
