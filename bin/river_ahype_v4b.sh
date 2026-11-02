#!/bin/bash
#
# Prepare AHYPE v4b time-varying river forcing (riverh_ahype_v4b_<spreading>_<YYYY>.[ab]).
#
# Must be run from the experiment directory (where EXPT.src lives).
# Reads EXPT.src and REGION.src for R (CONFIGNAME), E (IEXPT), T (topography version).
#
# Usage:
#   river_ahype_v4b.sh <spreading> <start_year> <end_year>
#   river_ahype_v4b.sh <spreading> <year>
#
# spreading: pointsource | coastal | diffuse
#
# Examples:
#   $HOME/NERSC-HYCOM-CICE/bin/river_ahype_v4b.sh coastal 2002 2022
#   $HOME/NERSC-HYCOM-CICE/bin/river_ahype_v4b.sh diffuse 2010

if [ $# -lt 2 ] || [ $# -gt 3 ] ; then
   echo ""
   echo "Generates AHYPE v4b river forcing (riverh_ahype_v4b_<spreading>_<YYYY>.[ab])"
   echo ""
   echo "Usage: $(basename $0) <spreading> <start_year> [<end_year>]"
   echo "  spreading: pointsource | coastal | diffuse"
   echo ""
   echo "Examples:"
   echo "  $(basename $0) coastal 2002 2022   # full period"
   echo "  $(basename $0) diffuse 2010         # single year"
   echo ""
   exit 1
fi

SPREADING=$1
START_YEAR=$2
END_YEAR=${3:-$2}

if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in the experiment directory."
   exit 1
fi

BINDIR=$(cd $(dirname $0) && pwd)
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src             || { echo "Could not source ./EXPT.src" ; exit 1 ; }

# R = CONFIGNAME, E = IEXPT, T = topography version (from src files)
DEPTH="${BASEDIR}/topo/depth_${R}_${T}.a"
GRID="${BASEDIR}/topo/regional.grid.a"
OUTDIR="${BASEDIR}/force/rivers/${E}"

if [ ! -s "$DEPTH" ] || [ ! -s "${DEPTH%.a}.b" ] ; then
   echo "Could not find depth files ${BASEDIR}/topo/depth_${R}_${T}.[ab]"
   exit 1
fi

mkdir -p "$OUTDIR"

cmd="python -u ${BINDIR}/river_ahype_v4b.py \
   --spreading  ${SPREADING}  \
   --start-year ${START_YEAR} \
   --end-year   ${END_YEAR}   \
   --outdir     ${OUTDIR}     \
   --grid       ${GRID}       \
   --depth      ${DEPTH}"

echo "$cmd"
eval $cmd || { echo "Error running river_ahype_v4b.py" ; exit 1 ; }
