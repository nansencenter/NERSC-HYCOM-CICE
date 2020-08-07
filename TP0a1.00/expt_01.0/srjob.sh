#!/bin/bash -l

#SBATCH --qos=devel
#SBATCH --account=nn9481k
#SBATCH --job-name=TP0a10
#SBATCH -t 00:20:00
#SBATCH -N 1   # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH  --mail-type=END
#SBATCH --mail-user=achref.othmani@nersc.no

#SBATCH -o log/HY_CICE.out
#SBATCH -e log/HY_CICE.err

#
#  Give the job a name
#
#         
module restore system
module load NCL/6.6.2-intel-2018b
module load FFTW/3.3.8-intel-2018b
module load Python/2.7.15-intel-2018b

export NMPI=4
export SLURM_SUBMIT_DIR=$(pwd)
# Enter directory from where the job was submitted
cd $SLURM_O_WORKDIR       ||  { echo "Could not go to dir $SLURM_O_WORKDIR  "; exit 1; }

# ------------------- Fetch Environment ------------------------------
# -------- these are needed in preprocess scripts---------------------
echo "SLURM_JOBID    = $SLURM_JOBID     "
echo "SLURM_JOBNAME  = $SLURM_JOBNAME   "
echo "SLURM_SUBMIT_DIR= $SLURM_SUBMIT_KDIR "
echo "SLURM_TASKNUM  = $SLURM_TASKNUM "
echo "SLURM_NUM_PPN  = $SLURM_NUM_PPN "
[ -z "$NOMP" ] && NOMP=0

# Enter directory from where the job was submitted
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_O_WORKDIR  "; exit 1; }

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

START="2000-01-01T00:00:00"
END="2000-01-20T00:00:00"
INITFLG="--init"
#INITFLG=""
echo "Start time in pbsjob.sh: $START"
echo "End   time in pbsjob.sh: $END"
# Generate atmospheric forcing :
atmo_synoptic.sh erai+all $START $END 

# Transfer data files to scratch - must be in "expt_XXX" dir for this script
expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

# Enter Scratch/run dir and Run model
cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
#srun -n $NMPI ./hycom_cice  > ../log/hycom.${SLURM_JOBID}.out 2>&1
srun -n $NMPI --cpu_bind=cores ./hycom_cice 

# Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
expt_postprocess.sh 

exit $?


