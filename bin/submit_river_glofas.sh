#!/bin/bash
#SBATCH --job-name=riverh_glofas
#SBATCH --account=nn2993k
#SBATCH --time=02:00:00
#SBATCH --qos=preproc
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -o log/riverh_glofas.%A_%a.out

# Slurm job array for GloFAS river forcing (one task per year).
#
# Experiment-specific config must be set in the calling environment (or in a
# thin wrapper script in the experiment directory) before submitting:
#
#   CONFIGNAME  e.g. TP2a0.10
#   IEXPT       e.g. 010
#   T           topography version, e.g. 04
#   DT          record interval in days: 1.0 for HYCOM 2.3, 0.25 for HYCOM 2.2.98
#
# Usage (from experiment directory):
#   sbatch --array=<FIRST_YEAR>-<LAST_YEAR>%12 \
#       $HOME/NERSC-HYCOM-CICE/bin/submit_river_glofas.sh <START> <END>
#
# START and END are dates in YYYY-MM-DD format; SLURM_ARRAY_TASK_ID is the year.
# The join step (concat_riverh.py) must be run separately once all years finish.

START=$1   # first day of the whole period, YYYY-MM-DD
END=$2     # last day of the whole period (included), YYYY-MM-DD

if [[ -z "$CONFIGNAME" || -z "$IEXPT" || -z "$T" || -z "$DT" ]]; then
    echo "Error: CONFIGNAME, IEXPT, T and DT must be set in the environment."
    exit 1
fi
if [[ -z "$START" || -z "$END" || -z "$SLURM_ARRAY_TASK_ID" ]]; then
    echo "Usage: sbatch --array=<FIRST_YEAR>-<LAST_YEAR>%12 submit_river_glofas.sh <START> <END>"
    exit 1
fi
set --    # conda's activate would otherwise see the positional arguments

source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh
module load Miniforge3/24.1.2-0
source ${EBROOTMINIFORGE3}/bin/activate
conda activate hycom-cice

DEPTH=$WORK/${CONFIGNAME}/topo/depth_${CONFIGNAME}_${T}.a

# This task's year. Run until 1 January of the next year, because a day's records
# are only written once the next day is processed; concat_riverh.py removes the overlap.
YEAR=$SLURM_ARRAY_TASK_ID
cs=${YEAR}-01-01
[[ "$cs" < "$START" ]] && cs=$START
ce=$((YEAR+1))-01-01
[[ "$ce" > "$END" ]] && ce=$END

OUT=$WORK/${CONFIGNAME}/force/rivers/${IEXPT}/SCRATCH/chunks_${START}_${END}/
mkdir -p $OUT
cd $HOME/NERSC-HYCOM-CICE/bin/river-topaz/scripts/GloFAS
python -u hycom_river_hifre.py --dt $DT ${cs}T00:00:00 ${ce}T00:00:00 $OUT $DEPTH
