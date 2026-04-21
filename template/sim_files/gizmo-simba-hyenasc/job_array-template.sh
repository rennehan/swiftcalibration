#!/bin/bash -l
#########################################################
#SBATCH -J JOB_NAME
#SBATCH -o ./slurm_files/slurm-%j.out

#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=ALL 
#SBATCH --account=rrg-babul-ad

##SBATCH --nodes=1 # 128/12.5
##SBATCH --nodes=1 # 256/25 + Hyenas L0
##SBATCH --nodes=6 # 512/50 + Hyenas L1
#SBATCH --nodes=NUMBER_NODES

#SBATCH --ntasks-per-node=192

##SBATCH --mem=0

#SBATCH --array=1-NUM_JOBS%1 # Run a N-job array, 1 job at a time
#SBATCH --time=24:0:0
#########################################################


# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "Job task $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_COUNT"
echo ""
# ---------------------------------------------------------------------


CURRPATH=.
EXEFILE=GIZMO_EXE
PARAMFILE=PARAMETER_FILE


RESTART_FLAG=0

if (( SLURM_ARRAY_TASK_ID > 1)); then	## CHANGE NUMBER!!!! 0=run from restart files; 1=run from ICs
    echo ""
    echo "Starting from restart files"
    echo ""
    RESTART_FLAG=1
fi


#########################################################
#module purge
#module load NiaEnv/2019b intel/2019u3 hdf5/1.10.5 gsl/2.5 openmpi/4.0.1
#module load NiaEnv/2019b intel/2019u3 hdf5/1.8.21 gsl/2.5 openmpi/4.0.1
#module load StdEnv/2023 intel/2023.2.1 openmpi/4.1.5 hdf5-mpi/1.14.2 gsl/2.7 fftw-mpi/3.3.10
#module load StdEnv/2023 intel/2023.2.1 openmpi/4.1.5 hdf5/1.14.2 hdf5-mpi/1.14.2 gsl/2.7 fftw/3.3.10 fftw-mpi/3.3.10
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 parmetis/4.0.3 fftw-mpi/3.3.10 gsl/2.7 hdf5-mpi/1.14.4 python/3.12.4	# Weiguang
which mpirun
#########################################################


mpirun -np $((SLURM_NTASKS_PER_NODE*SLURM_JOB_NUM_NODES)) --map-by core --report-bindings --mca orte_base_help_aggregate 0 --mca btl_openib_max_eager_rdma 0 --mca mpi_leave_pinned 0 --mca coll_hcoll_enable 0 ${CURRPATH}/${EXEFILE} ${CURRPATH}/${PARAMFILE} ${RESTART_FLAG}
