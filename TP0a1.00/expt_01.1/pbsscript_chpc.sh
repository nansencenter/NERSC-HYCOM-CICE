#!/bin/bash -l

#PBS -P ERTH0904
#PBS -N AGUa010

#PBS -l select=1:ncpus=4:mpiprocs=4

## System message output file
#PBS -o /mnt/lustre/users/asamuelsen/AGUa1.00/expt_01.0/log/HYCOM-CICE.out

## System error message file
#PBS -e /mnt/lustre/users/asamuelsen/AGUa1.00/expt_01.0/log/HYCOM-CICE.err

## How long job takes, wallclock time hh:mm:ss
#PBS -l walltime=0:30:00
#PBS -q test

##PBS -m abe
##PBS -M 

module load chpc/python/3.7.0
module load chpc/parallel_studio_xe/64/16.0.1/2016.1.150
module load chpc/compmech/fftw/3.3.10_intel2020u1
module load chpc/earth/netcdf/4.7.4/intel2020u1
module load chpc/compmech/netlib-lapack/3.10.1
module load chpc/earth/udunits2/udunits-2.2.28

export NMPI=$PBS_NP

# Enter directory from where the job was submitted
cd $PBS_O_WORKDIR       ||  { echo "Could not go to dir $PBS_O_WORKDIR  "; exit 1; }

# ------------------- Fetch Environment ------------------------------
# -------- these are needed in preprocess scripts---------------------
echo "PBS_JOBID    = $PBS_JOBID     "
echo "PBS_O_WORKDIR= $PBS_O_WORKDIR "
[ -z "$NOMP" ] && NOMP=0

# Enter directory from where the job was submitted
cd $PBS_O_WORKDIR       ||  { echo "Could not go to dir $PBS_O_WORKDIR  "; exit 1; }


# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

  START="2020-09-15T12:00:00"
  END="2020-09-19T12:00:00"
  INITFLG="--init"
#  INITFLG=""
#  ./bin/atmo_synoptic.sh era5+lw $START $END 
#  ./bin/expt_preprocess.sh $START $END $INITFLG

  echo "Start time in srjob.sh: $START"
  echo "End   time in srjob.sh: $END"

  # Generate atmospheric forcing :
#  atmo_synoptic.sh era5+lw $START $END 

  # Transfer data files to scratch - must be in "expt_XXX" dir for this script
  expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

  # Enter Scratch/run dir and Run model
  cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
  mpirun -np $NMPI ./hycom_alone  

  # Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
  cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
  expt_postprocess.sh 

done

exit $?



