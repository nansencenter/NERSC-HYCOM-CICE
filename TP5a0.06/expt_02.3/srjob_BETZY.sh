#!/bin/bash -l

##SBATCH --account=nn9481k
#SBATCH --account=nn2993k
##SBATCH --account=ns2993k

#SBATCH -J  TP5a011

#SBATCH  -N 5
#SBATCH  -n 640
##SBATCH -N 1   # number of nodes
##SBATCH -n 128   # number of cores

#SBATCH  --mail-type=END
#SBATCH --mail-user=tar@dmi.dk
 
## System message output file
#SBATCH -o log/HYCICE.out

## System error message file
#SBATCH -e log/HYCICE.err

## small job. NOrmal requires use of 4 node3s
# smaller jobs than 1 node
##SBATCH --qos=devel
## How long job takes, wallclock time hh:mm:ss
##SBATCH -t 30:00:00
#SBATCH -t 00:30:00

#ml purge
#module load Python/3.8.2-GCCcore-9.3.0
#module load ESMF/8.0.1-intel-2020a
#module load ESMF/8.3.0-intel-2022a-debug
#module load FFTW/3.3.8-intel-2020a
#module load UDUNITS/2.2.26-GCCcore-9.3.0
#module load CMake/3.16.4-GCCcore-9.3.0
#module load intel/2021a

ml purge
module load Python/3.8.2-GCCcore-9.3.0
#module load intel/2021a
#module load FFTW/3.3.8-intel-2020a
module load UDUNITS/2.2.26-GCCcore-9.3.0
#module load CMake/3.16.4-GCCcore-9.3.0
module load ESMF/8.3.0-intel-2022a-debug

module list
ulimit -s 200000

export NMPI=636
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
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "
  START="2005-01-01T00:00:00"
#  END="2005-01-02T00:00:00"
  END="2009-01-01T00:00:00"

  INITFLG="--init"   ### Used for Clim Run
#  INITFLG=" "        ### Used for Nest Run

  echo "Start time in pbsjob.sh: $START"
  echo "End   time in pbsjob.sh: $END"

# Generate atmospheric forcing :
#if [ ${create_atm} ]
#then
   echo 'create atmospheric forcing'
#  ${NHCROOT}/bin/atmo_synoptic.sh era5+lw $START $END 
#else
#   echo 'do not create atmospheric forcing'
#fi
  # Transfer data files to scratch - must be in "expt_XXX" dir for this script
echo ${NHCROOT}
  ${NHCROOT}/bin/expt_preprocess.sh $START $END $INITFLG ||  { echo "Preprocess had fatal errors "; exit 1; }
  cp ${P}/SCRATCH/ice_in ${P}/SCRATCH/ice_in2 
             cat ${P}/SCRATCH/ice_in2                \
             | sed s%'.* ice_ic .*'%" ice_ic = 'read' "% \
             | sed s%'.* ocn_data_dir .*'%"   ocn_data_dir = './' "% \
             > ${P}/SCRATCH/ice_in
# Enter Scratch/run dir and Run model
  cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
  srun --mpi=pmi2 -n $NMPI --cpu_bind=cores ./hycom_cice_nersc

  # Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
  cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
#  /cluster/home/tar/NERSC-HOME-CICE/bin/expt_postprocess.sh 
  /cluster/home/tar/hycomv2_3/NERSC-HYCOM-CICE/bin/expt_postprocess.sh
exit $?

