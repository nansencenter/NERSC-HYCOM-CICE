#!/bin/bash
#SBATCH --job-name=riverh_glofas
#SBATCH --account=nn2993k
#SBATCH --time=02:00:00
#SBATCH --qos=preproc
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -o log/riverh_glofas.%A_%a.out

# GloFAS river forcing riverh.[ab] for a long period, as a Slurm job array (one task per
# year) followed by a join job (concat_riverh.py). Gives the same file as river_glofas.sh.
#
# Usage, from the experiment directory. Run it directly, not with sbatch: it submits the
# Slurm jobs itself (the #SBATCH lines above are for those jobs):
#   $HOME/NERSC-HYCOM-CICE/bin/submit_river_glofas.sh [-s spreading] <START> <END> [dt]
#
# START and END are the first and last day (included) of the period, YYYY-MM-DD.
# Region, experiment number, topography and HYCOM version are read from REGION.src and
# EXPT.src. dt (fraction of a day) defaults to 1.0 for HYCOM 2.3, 0.25 otherwise.
# -s: spreading of each river mouth (bin/river_spreading.py), pointsource, coastal or diffuse
# (default diffuse).
# The result is written to force/rivers/<E>/riverh.[ab].

# The spreading reaches the Slurm jobs through the environment (sbatch exports it)
SPREADING=${SPREADING:-diffuse}
while getopts "s:" opt ; do
    case $opt in
        s) SPREADING=$OPTARG ;;
        *) exit 1 ;;
    esac
done
shift $((OPTIND-1))
export SPREADING

START=$1
END=$2
dt=$3
MODE=${4:-submit}   # submit (user), year or join (Slurm jobs submitted by this script)

if [[ -z "$START" || -z "$END" ]]; then
    echo "Usage (from the experiment directory): $(basename $0) [-s spreading] <START> <END> [dt]"
    echo "   START, END: first and last day (included), YYYY-MM-DD"
    exit 1
fi
if [ ! -f EXPT.src ]; then
    echo "Could not find EXPT.src. This script must be run in the experiment directory"
    exit 1
fi
set --    # conda's activate and EXPT.src would otherwise see the positional arguments
BASEDIR=$(cd .. && pwd)
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
if [ -z "$dt" ]; then
    if [ "$V" == "2.3" ]; then dt=1.0; else dt=0.25; fi
fi
RDIR=${BASEDIR}/force/rivers/${E}
CHUNKS=${RDIR}/SCRATCH/chunks_${START}_${END}

# --- Submit the job array and the join job
if [ "$MODE" == "submit" ]; then
    # Inside a Slurm job $0 is a copy of this script, so pass its real path on
    export GLOFAS_SCRIPT=$(cd $(dirname $0) && pwd)/$(basename $0)
    mkdir -p log
    jid=$(sbatch --parsable --array=${START:0:4}-${END:0:4}%12 \
          $GLOFAS_SCRIPT $START $END $dt year) || exit 1
    jid_join=$(sbatch --parsable --dependency=afterok:$jid --job-name=riverh_join \
          --time=00:30:00 --mem=2G --output=log/riverh_join.%j.out \
          $GLOFAS_SCRIPT $START $END $dt join) || exit 1
    echo "Region $R, experiment $E, topography $T, dt=$dt day, spreading $SPREADING"
    echo "Array job (one task per year): $jid"
    echo "Join job (starts when all years succeed): $jid_join"
    echo "Result: $RDIR/riverh.[ab]"
    exit 0
fi

source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh
module load Miniforge3/24.1.2-0
source ${EBROOTMINIFORGE3}/bin/activate
conda activate hycom-cice
SCRIPTDIR=$(dirname $GLOFAS_SCRIPT)

# --- Join the years into riverh.[ab]
if [ "$MODE" == "join" ]; then
    python -u ${SCRIPTDIR}/concat_riverh.py $RDIR/riverh ${CHUNKS}/riverh_*.a
    exit $?
fi

# --- One year (array task). Use the experiment's depth file, as river_glofas.sh does.
depthfile=${BASEDIR}/topo/depth_${R}_${T}.a
if [ ! -s "$depthfile" ] || [ ! -s "${depthfile%.a}.b" ]; then
    echo "Could not find depth files ${BASEDIR}/topo/depth_${R}_${T}.[ab]"
    exit 1
fi

# Run until 1 January of the next year, because a day's records are only written once
# the next day is processed; concat_riverh.py removes the overlap.
YEAR=$SLURM_ARRAY_TASK_ID
cs=${YEAR}-01-01
[[ "$cs" < "$START" ]] && cs=$START
ce=$((YEAR+1))-01-01
[[ "$ce" > "$END" ]] && ce=$END

mkdir -p $CHUNKS
cd ${SCRIPTDIR}/river-topaz/scripts/GloFAS
python -u hycom_river_hifre.py --dt $dt --spreading $SPREADING ${cs}T00:00:00 ${ce}T00:00:00 $CHUNKS/ $depthfile
