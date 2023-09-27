#!/bin/bash -l
#########################################################
#SBATCH -J gen_ssfr
#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=FAIL
#SBATCH --account=rrg-babul-ad
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/slurm_files/slurm-%j.out
#SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations/swimba_s25n256_ps2020_64x/slurm_files/slurm-%j.out
##SBATCH -o /scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations/swimba_calibration/slurm_files/slurm-%j.out
#########################################################
#SBATCH --time=2:00:0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#########################################################

conda activate yt_swift_caesar

root_dir=/scratch/b/babul/aspadawe/swift_tests/cali_simba/past_calibrations
#cali_dir=/swimba_s25n256_ps2020_64x/calibrations
cali_dir=/swimba_s25n256_ps2020_64x
cali_subdir=/cali_*
sim_name=simba_s25n256

full_dir=$root_dir$cali_dir

for FILE in $full_dir$cali_subdir; do
        echo $FILE
        python gen_ssfr.py --modeldir=$full_dir --simdir=$(basename $FILE) --sim=$sim_name
done

#########################################################
