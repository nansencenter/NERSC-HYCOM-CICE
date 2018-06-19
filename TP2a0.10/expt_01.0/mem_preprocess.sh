#!/bin/bash -l

module restore system
module load NCL/6.4.0-intel-2017a
module load FFTW/3.3.6-intel-2017a

export SLURM_SUBMIT_DIR=$(pwd)
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }

# ------------------- Fetch Environment ------------------------------
# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }


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
      START="2001-12-28T00:00:00"
      END="2002-01-04T00:00:00"
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


done

exit $?

