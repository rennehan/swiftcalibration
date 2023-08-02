#!/bin/bash -l
#########################################################
#SBATCH -J create_caesarfiles
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations/swimba_s25n256_ps2020_64x/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/s18n128_simba/slurm-%j.out
#########################################################
#SBATCH --time=6:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar

root_dir=/scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations
cali_dir=/swimba_s25n256_ps2020_64x
calis=/cali_*

full_dir=$root_dir$cali_dir

for FILE in $full_dir$calis; do
	echo $FILE
	python create_caesarfiles.py $full_dir $(basename $FILE)
done

for FILE in $full_dir$calis; do
	echo $FILE
        python gen_gsmf.py $full_dir $(basename $FILE)
done

for FILE in $full_dir$calis; do
	echo $FILE
        python gen_bhmsm.py $full_dir $(basename $FILE)
done

#########################################################
