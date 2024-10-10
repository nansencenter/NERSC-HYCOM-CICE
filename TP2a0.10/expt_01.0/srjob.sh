#!/bin/bash -l

#SBATCH --account=nn9481k 
#SBATCH --job-name=TP2a010
#SBATCH --time="00:06:00"
#SBATCH --nodes=4 # number of nodes
#SBATCH --ntasks=504 # number of cores
#SBATCH --mail-type=END

#SBATCH -o log/HY_CICE.out
#SBATCH -e log/HY_CICE.err

#         
export NMPI=504
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
source $NHCROOT/betzy_env.sh || { echo "Could not source ~/betzy_env.sh "; exit 1; }
#source $NHCROOT/fram_env.sh || { echo "Could not source ~/fram_env.sh "; exit 1; }

echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

START="1990-01-01T00:00:00"
END="1990-01-05T00:00:00"
INITFLG="--init"
#INITFLG=""
echo "Start time in pbsjob.sh: $START"
echo "End   time in pbsjob.sh: $END"
# Generate atmospheric forcing :
#atmo_synoptic.sh erai+all $START $END 
../bin/atmo_synoptic.sh era5+lw $START $END

# Transfer data files to scratch - must be in "expt_XXX" dir for this script
# expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 2; }
../bin/expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

# Enter Scratch/run dir and Run model
cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
#srun -n $NMPI ./hycom_cice  > ../log/hycom.${SLURM_JOBID}.out 2>&1
srun -n $NMPI --cpu_bind=cores ./hycom_cice 

# Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
expt_postprocess.sh 

exit $?


