#!/bin/bash -l
#########################################################
#SBATCH -J gen_observables
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations/swimba_s50n512_ps2020_sphenix_32x/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/s18n128_simba/slurm-%j.out
#########################################################
#SBATCH --time=10:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar

root_dir=/scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations
cali_dir=/swimba_s50n512_ps2020_sphenix_32x/calibrations
cali_subdir=/cali_*
#cali_subdir=(/cali_0006 /cali_0012 /cali_0025 /cali_0030 /cali_0031)
# cali_subdir=(/cali_0007)
sim_name=simba_s50n512

full_dir=$root_dir$cali_dir


# for FILE in $full_dir${cali_subdir[@]}; do
for FILE in $full_dir$cali_subdir; do
	echo $FILE
	python create_caesarfiles.py --modeldir=$full_dir --simdir=$(basename $FILE) --sim=$sim_name
	python gen_gsmf.py --modeldir=$full_dir --simdir=$(basename $FILE) --sim=$sim_name
	python gen_ssfr.py --modeldir=$full_dir --simdir=$(basename $FILE) --sim=$sim_name
	python gen_bhmsm.py --modeldir=$full_dir --simdir=$(basename $FILE) --sim=$sim_name
	python gen_csfh.py --modeldir=$full_dir --simdir=$(basename $FILE) --sim=$sim_name
done

#########################################################
