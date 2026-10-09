#!/bin/bash -l

#SBATCH --account=nn9481k
#SBATCH --job-name=TP2a010
#SBATCH --time="48:00:00"
#SBATCH --nodes=4 # number of nodes
#SBATCH --ntasks=504 # number of cores

#SBATCH -o log/HY_CICE.%J.out
#SBATCH -e log/HY_CICE.%J.err

#
# Cycled spin-up: run the period CYCLE_START..CYCLE_END NCYCLES times with the same
# (GLORYS) boundary and atmospheric forcing. blkdat.input must have the GLORYS nesting
# settings. At the end of each cycle the restart files valid at CYCLE_END are wrapped
# to CYCLE_START (bin/cycle_wrap.py) and the next cycle starts from them.
#
# Cycle NN writes to its own data directory data/cycle_NN (D in EXPT.src must contain
# ${SPINUP_CYCLE:+/cycle_${SPINUP_CYCLE}}). The first cycle starts at SPINUP_START, which can lie
# inside the period (e.g. a September start from a GLORYS state). It needs the HYCOM
# restart for SPINUP_START in data/cycle_01; without a CICE restart at SPINUP_START,
# CICE is cold-started from ice_initial.nc (INITFLG="--init-ice").
#
# The job can be resubmitted at any time: it continues from the latest restart pair
# found in the cycle data directories. It runs up to CYCLES_PER_JOB cycles per job (an
# interrupted cycle is continued from its latest restart) and resubmits itself
# automatically until all cycles are done.
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
source $NHCROOT/environment/betzy_env.sh || { echo "Could not source betzy_env.sh "; exit 1; }
module load Miniforge3/24.1.2-0          || { echo "Could not load Miniforge3 module"; exit 1; }
source ${EBROOTMINIFORGE3}/bin/activate  || { echo "Could not activate Miniforge3 base"; exit 1; }
conda activate hycom-cice                || { echo "Could not activate hycom-cice conda environment"; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

# ------------------- Cycle settings ---------------------------------
SPINUP_START="1993-09-01T00:00:00"   # start of cycle 1
CYCLE_START="1993-01-01T00:00:00"    # start of cycles 2..NCYCLES
CYCLE_END="1998-01-01T00:00:00"      # end of every cycle (wrap date)
NCYCLES=5
CYCLES_PER_JOB=1                     # cycles per job; one cycle must fit in the wall-time limit
JOBSCRIPT=srjob_cycle.sh             # this script, for resubmission
ATMO_FORCING="era5+lw"               # forcing option for atmo_synoptic.sh; "" to never generate forcing
# --------------------------------------------------------------------

EXPTDIR=$P

# Seconds since epoch of an ISO date YYYY-mm-ddTHH:MM:SS
epoch() { date -u -d "${1/T/ } UTC" +%s ; }

# ISO date of HYCOM restart suffix YYYY_DDD_HH_SSSS
restart_iso() {
   local y=${1:0:4} d=${1:5:3} h=${1:9:2} s=${1:12:4}
   date -u -d "$y-01-01 00:00:00 UTC +$((10#$d-1)) days +$((10#$h)) hours +$((10#$s)) seconds" +%Y-%m-%dT%H:%M:%S
}

# CICE restart file (relative to data dir) valid at ISO date
cice_restart() {
   local dsec=$(( $(epoch $1) % 86400 ))
   printf "cice/iced.%s-%05d.nc" ${1:0:10} $dsec
}

# Latest date in [first,last] with a HYCOM restart in data dir $1 and a CICE restart
# (or no CICE restart needed: date == $4). Prints nothing if there is none.
latest_restart() {
   local dir=$1 first=$(epoch $2) last=$(epoch $3) coldice=$4 best="" bestsec=0 f iso sec
   for f in $dir/restart.[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9][0-9][0-9].b ; do
      [ -f "$f" ] || continue
      [ -f "${f%.b}.a" ] || continue
      iso=$(restart_iso $(basename $f .b | sed "s/^restart\.//"))
      sec=$(epoch $iso)
      [ $sec -ge $first -a $sec -le $last ] || continue
      [ -f $dir/$(cice_restart $iso) -o "$iso" == "$coldice" ] || continue
      if [ $sec -gt $bestsec ] ; then best=$iso ; bestsec=$sec ; fi
   done
   echo $best
}

[ $(epoch $CYCLE_START) -lt $(epoch $CYCLE_END) ] || { echo "CYCLE_START must be before CYCLE_END"; exit 1; }
[ $(epoch $SPINUP_START) -ge $(epoch $CYCLE_START) -a $(epoch $SPINUP_START) -lt $(epoch $CYCLE_END) ] || \
   { echo "SPINUP_START must lie in [CYCLE_START,CYCLE_END)"; exit 1; }

# Atmospheric forcing is the same for all cycles: generate it once for CYCLE_START..CYCLE_END
# unless the existing forcing in force/synoptic/$E already covers that period.
# Forcing times in the .b files are in HYCOM days (days since 1900-12-31).
hycom_day() { echo $(( ($(epoch $1) - $(epoch 1900-12-31T00:00:00)) / 86400 )) ; }
forcing_covers() {
   local b=$1 first last
   [ -s $b -a -s ${b%.b}.a ] || return 1
   first=$(head -n 6 $b | tail -n 1 | sed "s/.*=//" | awk '{print $1}')
   last=$(tail -n 1 $b | sed "s/.*=//" | awk '{print $1}')
   awk -v f=$first -v l=$last -v s=$2 -v e=$3 'BEGIN { exit !(f <= s && l >= e) }'
}
FORCDIR=$EXPTDIR/../force/synoptic/$E
if [ -n "$ATMO_FORCING" ] ; then
   cs=$(hycom_day $CYCLE_START) ; ce=$(hycom_day $CYCLE_END) ; ok=1
   for v in radflx shwflx vapmix airtmp precip mslprs wndewd wndnwd dewpt ; do
      forcing_covers $FORCDIR/$v.b $cs $ce || { ok=0 ; break ; }
   done
   if [ $ok -eq 1 ] ; then
      echo "Atmospheric forcing in $FORCDIR covers $CYCLE_START..$CYCLE_END"
   else
      echo "Generating atmospheric forcing ($ATMO_FORCING) for $CYCLE_START..$CYCLE_END in $FORCDIR"
      cd $EXPTDIR || { echo "Could not go to dir $EXPTDIR"; exit 1; }
      ../bin/atmo_synoptic.sh $ATMO_FORCING $CYCLE_START $CYCLE_END || { echo "atmo_synoptic.sh failed"; exit 1; }
   fi
fi

PREV_D=""
NDONE=0                              # cycles run in this job
for c in $(seq 1 $NCYCLES) ; do
   export SPINUP_CYCLE=$(printf "%02d" $c)
   cd $EXPTDIR || { echo "Could not go to dir $EXPTDIR"; exit 1; }
   source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
   [[ "$D" == */cycle_${SPINUP_CYCLE} ]] || { echo "D=$D does not end in cycle_${SPINUP_CYCLE}. Add \${SPINUP_CYCLE:+/cycle_\${SPINUP_CYCLE}} to D in EXPT.src"; exit 1; }

   # CICE restart only optional at SPINUP_START of cycle 1 (cold start of CICE)
   if [ $c -eq 1 ] ; then first=$SPINUP_START ; coldice=$SPINUP_START ; else first=$CYCLE_START ; coldice=none ; fi
   cur=$(latest_restart $D $first $CYCLE_END $coldice)
   if [ -z "$cur" ] ; then
      if [ $c -eq 1 ] ; then
         echo "Cycle $SPINUP_CYCLE: no HYCOM restart for $SPINUP_START in $D" ; exit 1
      fi
      echo "Cycle $SPINUP_CYCLE: wrapping restarts $CYCLE_END in $PREV_D to $CYCLE_START in $D"
      python ../bin/cycle_wrap.py $PREV_D $D $CYCLE_END $CYCLE_START --blkdat blkdat.input --ice-in ice_in \
         || { echo "cycle_wrap.py failed"; exit 1; }
      cur=$CYCLE_START
   fi
   echo "Cycle $SPINUP_CYCLE: continuing from $cur"

   if [ $(epoch $cur) -lt $(epoch $CYCLE_END) ] ; then
      START=$cur
      END=$CYCLE_END
      INITFLG=""
      if [ $c -eq 1 -a "$START" == "$SPINUP_START" -a ! -f $D/$(cice_restart $START) ] ; then
         INITFLG="--init-ice"
      fi
      echo "Cycle $SPINUP_CYCLE: Start time in pbsjob.sh: $START"
      echo "Cycle $SPINUP_CYCLE: End   time in pbsjob.sh: $END"
      echo "Cycle $SPINUP_CYCLE: INITFLG=$INITFLG"

      # Transfer data files to scratch - must be in "expt_XXX" dir for this script
      cd $EXPTDIR || { echo "Could not go to dir $EXPTDIR"; exit 1; }
      ../bin/expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

      # Enter Scratch/run dir and Run model
      cd $S  ||  { echo "Could not go to dir $S  "; exit 1; }
      srun -n $NMPI --cpu_bind=cores ./hycom_cice

      # Cleanup and move data files to data directory - must be in "expt_XXX" dir for this script
      cd $EXPTDIR     ||  { echo "Could not go to dir $EXPTDIR  "; exit 1; }
      ../bin/expt_postprocess.sh

      [ "$(latest_restart $D $END $END none)" == "$END" ] || \
         { echo "Cycle $SPINUP_CYCLE: no HYCOM/CICE restart for $END in $D after cycle, stopping"; exit 1; }
      NDONE=$((NDONE+1))
      if [ $c -lt $NCYCLES -a $NDONE -ge $CYCLES_PER_JOB ] ; then
         echo "Cycle $SPINUP_CYCLE finished, resubmitting $JOBSCRIPT"
         cd $EXPTDIR && sbatch $JOBSCRIPT
         exit 0
      fi
   fi
   echo "Cycle $SPINUP_CYCLE finished"
   PREV_D=$D
done
echo "All $NCYCLES cycles finished"
exit 0
