#!/bin/bash -l
#########################################################
#SBATCH -J gen_observables
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o ../template/slurm_files/slurm-%j.out
#########################################################
#SBATCH --time=10:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar

root_dir=../template
cali_dir=/calibrations
cali_subdir=/cali_*
#cali_subdir=(/cali_0006 /cali_0012 /cali_0025 /cali_0030 /cali_0031)
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
