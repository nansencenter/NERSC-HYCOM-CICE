#!/bin/bash
#SBATCH --job-name=riverh_ahype_v4b
#SBATCH --account=nn2993k
#SBATCH --time=04:00:00
#SBATCH --qos=preproc
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -o log/riverh_ahype_v4b.%A_%a.out

# AHYPE v4b river forcing — one calendar year per array task.
#
# Do not submit this script directly. Instead, create a thin wrapper in your
# experiment directory (see docs/forcing.md for a template) that sets the
# experiment-specific variables and calls this script:
#
#   sbatch --array=2002-2022%12 riverh_ahype_v4b_job.sh
#
# Required environment variables (set by the wrapper):
#   CONFIGNAME   e.g. TP2a0.10
#   IEXPT        e.g. 030
#   T            topography version, e.g. 04
#   SPREADING    pointsource | coastal | diffuse
#
# Spreading strategies:
#   pointsource  No spreading;                       ~5 min/year,  8 GB
#   coastal      50 km / 100 km, σ=1 cell;          ~30 min/year, 32 GB
#   diffuse      180 km / 360 km, σ=3 cells;        ~90 min/year, 32 GB
#
# Output: $WORK/<CONFIGNAME>/force/rivers/<IEXPT>/riverh_ahype_v4b_<spreading>_<YYYY>.[ab]
# Join with bin/concat_riverh.py after all tasks finish.

set -euo pipefail

for var in CONFIGNAME IEXPT T SPREADING; do
    [[ -z "${!var:-}" ]] && { echo "Error: $var not set. Run via the experiment wrapper."; exit 1; }
done

[[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] && { echo "Error: must run as a Slurm job array."; exit 1; }

YEAR=$SLURM_ARRAY_TASK_ID
BINDIR=$(cd "$(dirname "$0")" && pwd)
DEPTH=$WORK/${CONFIGNAME}/topo/depth_${CONFIGNAME}_${T}.a
GRID=$WORK/${CONFIGNAME}/topo/regional.grid.a
OUTDIR=$WORK/${CONFIGNAME}/force/rivers/${IEXPT}

[[ ! -s "$DEPTH" ]] && { echo "Depth file not found: $DEPTH"; exit 1; }
mkdir -p "$OUTDIR" log

set --   # conda activate reads $@ — clear before sourcing
source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh
module load Miniforge3/24.1.2-0
source ${EBROOTMINIFORGE3}/bin/activate
conda activate hycom-cice

echo "Year=${YEAR}  spreading=${SPREADING}  outdir=${OUTDIR}"

python -u "${BINDIR}/river_ahype_v4b.py" \
    --spreading "$SPREADING" \
    --year      "$YEAR"      \
    --outdir    "$OUTDIR"    \
    --grid      "$GRID"      \
    --depth     "$DEPTH"
