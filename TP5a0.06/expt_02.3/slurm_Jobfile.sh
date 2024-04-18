#!/bin/bash -l

#SBATCH --account=nn9481k

#SBATCH -J  TP5a010

#SBATCH --nodes 11   # number of nodes
#SBATCH -n 31  # number of tasks per note (1 task per cpu)
## Each compute node has 24 cores (See more details in section Hardware on Sisu User Guide).

##SBATCH -p test_large 
## Choose a suitable queue <test,small,large>
## How to check queue limits: scontrol show part <queue name>
## for example: scontrol show part small

## System message output file
#SBATCH -o log/HYCICE.%J.out

## System error message file
#SBATCH -e log/HYCICE.%J.err

## How long job takes, wallclock time hh:mm:ss
#SBATCH -t 00:10:00
##SBATCH --qos=devel

#SBATCH  --mail-type=END
#SBATCH --mail-user=alfatih.ali@nersc.no

## option: -n (total number of mpi processes)
## option: -N (number of mpi processes per compute node)
## option: -S (number of mpi processes per NUMA node)
## option: -ss (allocate memory only from a local NUMA node)

## Calculate the total number of cores and store it in variable ncores
##(( ncores = SLURM_NNODES * 24 ))
##(( ncores = SLURM_NNODES * 32 ))  #on Fram


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
export NMPI=336
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

for ik in `seq 3 3` ; do
  if [ $ik -eq 0 ] ; then
    START="2004-01-01T00:00:00"
    END="2006-01-01T00:00:00"
    INITFLG="--init"
  #  cp -f ice_in.0 ice_in
  else
    INITFLG=""
    cp -f ice_in.1 ice_in
    if [ $ik -eq 1 ] ; then
      START="1989-01-22T00:00:00"
      END="1989-02-02T00:00:00"
    else 
      START="2010-01-01T00:00:00"
      END="2013-01-01T00:00:00"
      #START="2006-08-02T00:00:00"
      #END="2010-01-01T00:00:00"
    fi
  fi
  echo "Start time in pbsjob.sh: $START"
  echo "End   time in pbsjob.sh: $END"

  # Generate atmospheric forcing :
  #atmo_synoptic.sh erai $START $END 

  # Transfer data files to scratch - must be in "expt_XXX" dir for this script
  expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

  # Enter Scratch/run dir and Run model
  cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
  ##aprun -n $NMPI -m 500M ./hycom_cice  > ../log/hycom.%J.out 2>&1
  #srun -n $NMPI --cpu_bind=cores ./hycom_cice  
  srun --mpi=pmi2 -n $NMPI --cpu_bind=cores ./hycom_cice ${e}
  # Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
  cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
  expt_postprocess.sh 

done

exit $?

