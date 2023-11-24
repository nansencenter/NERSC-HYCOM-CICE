#!/bin/bash -l

#SBATCH --account=nn9481k

#SBATCH -J  TP2a010

#SBATCH -N 10   # number of nodes
#SBATCH -n 22   # number of nodes
## Each compute node has 24 cores (See more details in section Hardware on Sisu User Guide).

##SBATCH -p small_long
## Choose a suitable queue <test,small,large>
## How to check queue limits: scontrol show part <queue name>
## for example: scontrol show part small

## System message output file
#SBATCH -o log/HYCICE.%J.out

## System error message file
#SBATCH -e log/HYCICE.%J.err

## How long job takes, wallclock time hh:mm:ss
#SBATCH -t 00:30:00

#SBATCH  --mail-type=END
#SBATCH --mail-user=jiping.xie@nersc.no

## option: -n (total number of mpi processes)
## option: -N (number of mpi processes per compute node)
## option: -S (number of mpi processes per NUMA node)
## option: -ss (allocate memory only from a local NUMA node)

## Calculate the total number of cores and store it in variable ncores
##(( ncores = SLURM_NNODES * 24 ))
##(( ncores = SLURM_NNODES * 32 ))  #on Fram

module restore system
module load FFTW/3.3.10-GCC-11.3.0
module load ESMF/8.3.0-intel-2022a
module load Python/3.10.4-GCCcore-11.3.0
module load UDUNITS/2.2.28-GCCcore-11.3.0
module load intel/2022a

export NMPI=220
export SLURM_SUBMIT_DIR=$(pwd)
# Enter directory from where the job was submitted
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }

# ------------------- Fetch Environment ------------------------------
# -------- these are needed in preprocess scripts---------------------
echo "SLURM_JOBID    = $SLURM_JOBID     "
echo "SLURM_SUBMIT_DIR= $SLURM_SUBMIT_DIR "
[ -z "$NOMP" ] && NOMP=0

# Enter directory from where the job was submitted
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_O_WORKDIR  "; exit 1; }


# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "


for ik in `seq 1 1` ; do
  if [ $ik -eq 0 ] ; then
    START="2001-12-25T00:00:00"
    END="2002-01-01T00:00:00"
    INITFLG="--init"
    cp -f ice_in.0 ice_in
  else
    INITFLG=""
    cp -f ice_in.1 ice_in
    if [ $ik -eq 1 ] ; then
      START="2001-12-25T00:00:00"
      END="2002-01-02T00:00:00"
    else 
      START="2010-01-02T00:00:00"
      END="2016-01-02T00:00:00"
    fi
  fi
  echo "Start time in pbsjob.sh: $START"
  echo "End   time in pbsjob.sh: $END"

  # Generate atmospheric forcing :
  atmo_synoptic.sh erai $START $END 

  # Transfer data files to scratch - must be in "expt_XXX" dir for this script
  expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

  # Enter Scratch/run dir and Run model
  cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
  ##aprun -n $NMPI -m 500M ./hycom_cice  > ../log/hycom.%J.out 2>&1
  srun --mpi=pmi2 -n $NMPI --cpu_bind=cores ./hycom_cice  

  # Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
  cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
  expt_postprocess.sh 

done

exit $?

