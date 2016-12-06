#! /bin/bash
#
#  Make sure I use the correct shell.
#
#PBS -S /bin/bash
#
#  Give the job a name
#
#PBS -N TP0a1.00_X01.0
#
#  Specify the project the job belongs to
#
#PBS -A nn2993k
#
#  We want 24 hours on 32 cpu's:
#
##PBS -l walltime=00:40:00,mppwidth=4
#PBS -l walltime=00:15:00,mppwidth=4
#
#  The job needs 500mb  memory per cpu:
#PBS -l mppmem=500mb
#
#  Send me an email on  a=abort, b=begin, e=end
#
#PBS -m a
#
#  Use this email address (check that it is correct):
#PBS -M knut.lisaeter@nersc.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o log/mpijob.${PBS_JOBID}.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e log/mpijob.${PBS_JOBID}.err
#

# ------------------- Fetch Environment ------------------------------
# -------- these are needed in preprocess scripts---------------------
MYMPI=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
MYOMP=`qstat -f $PBS_JOBID | awk '/mppdepth/ {print $3}'`
echo "PBS_JOBID    = $PBS_JOBID     "
echo "PBS_JOBNAME  = $PBS_JOBNAME   "
echo "PBS_O_WORKDIR= $PBS_O_WORKDIR "
echo "PBS_TASKNUM  = $PBS_TASKNUM "
echo "PBS_NUM_PPN  = $PBS_NUM_PPN "
[ -z "$NOMP" ] && NOMP=0
echo "MYMPI=$MYMPI  (Number of system MPI tasks requested)"
echo "MYOMP=$MYOMP  (Number of systems OMP threads requested" # Not really used yet
# -----------------End Fetch Environment -----------------------------

# Enter directory from where the job was submitted
cd $PBS_O_WORKDIR       ||  { echo "Could not go to dir $PBS_O_WORKDIR  "; exit 1; }

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

START="2000-01-01T00:00:00"
END="2000-01-10T00:00:00"
INITFLG="--init"
#INITFLG=""
echo "Start time in pbsjob.sh: $START"
echo "End   time in pbsjob.sh: $END"

# Generate atmospheric forcing :
atmo_synoptic.sh erai $START $END 

# Transfer data files to scratch - must be in "expt_XXX" dir for this script
expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

# Enter Scratch/run dir and Run model
cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
aprun -n $NMPI -m 500M ./hycom_cice  > ../log/hycom.${PBS_JOBID}.out 2>&1

# Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
cd $P     ||  { echo "Could not go to dir $P  "; exit 1; }
expt_postprocess.sh 

exit $?

