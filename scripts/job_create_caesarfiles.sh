#!/bin/bash -l
#########################################################
#SBATCH -J create_caesarfiles
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations/swimba_s25n256_ps2020_64x_v2/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/s50n512_simba_ps2020/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/s50n512_simba_ps2020_calibration6_v2_mpi/slurm_files/slurm-%j.out
#SBATCH -o /scratch/b/babul/aspadawe/swift_tests/s50n512_simba_ps2020_calibration6_v2_mpi/slurm_files/slurm-%j.out
#########################################################
#SBATCH --time=12:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar

#root_dir=/scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations
root_dir=/scratch/b/babul/aspadawe/swift_tests
#cali_dir=/swimba_s25n256_ps2020_64x_v2/calibrations
#cali_dir=/s50n512_simba_ps2020_sphenix_calibration7_mpi
cali_dir=/s50n512_simba_ps2020_calibration6_v2_mpi
#cali_dir=/s50n512_simba_ps2020
#cali_subdir=/cali_*
cali_subdir=/.
sim_name=simba_s50n512

full_dir=$root_dir$cali_dir


for FILE in $full_dir$cali_subdir; do
	echo $FILE
	python create_caesarfiles.py --modeldir=$full_dir --simdir=$(basename $FILE) --sim=$sim_name
done

#########################################################
