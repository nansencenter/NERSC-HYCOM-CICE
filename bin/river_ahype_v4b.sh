#!/bin/bash
#
# Generates AHYPE v4b river forcing (monthly, AHYPE v4b + E-HYPE + GEUS ice discharge) and
# writes it to force/rivers/<E>/riverh.[ab]. Must be run in the experiment directory.
# For a Slurm job instead of the login node, use submit_river_ahype_v4b.sh.

usage() {
   echo ""
   echo "Generates AHYPE v4b river forcing (riverh.[ab]) for whole years"
   echo ""
   echo " $(basename $0) [-s spreading] [first_year] [last_year]"
   echo "    -s: spreading (bin/river_spreading.py), pointsource | coastal | diffuse (default diffuse)"
   echo "    years default to the whole AHYPE v4b period, 2002 2022"
   echo ""
   echo "Example:"
   echo "    $(basename $0)"
   echo "    $(basename $0) -s coastal 2010 2012"
   echo ""
   exit 1
}
spreading=diffuse
while getopts "s:h" opt ; do
   case $opt in
      s) spreading=$OPTARG ;;
      *) usage ;;
   esac
done
shift $((OPTIND-1))
[ $# -gt 2 ] && usage
first=${1:-2002}
last=${2:-2022}

if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
set --    # EXPT.src would otherwise see the positional arguments
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }

D=$BASEDIR/force/rivers/$E
S=$D/SCRATCH
mkdir -p $S

# Use the experiment's depth file (topography version T in EXPT.src), as river_glofas.sh does
depthfile="${BASEDIR}/topo/depth_${R}_${T}.a"
if [ ! -s "$depthfile" ] || [ ! -s "${depthfile%.a}.b" ]; then
   echo "Could not find depth files ${BASEDIR}/topo/depth_${R}_${T}.[ab]"
   exit 1
fi

out=$S/riverh_ahype_v4b_${spreading}_${first}_${last}
cmd="python -u ${BINDIR}river_ahype_v4b.py --spreading $spreading --start-year $first --end-year $last \
   --grid ${BASEDIR}/topo/regional.grid.a --depth $depthfile --out $out"
echo "$cmd"
eval $cmd || { echo "Error running $cmd" ; exit 1 ; }

for ext in a b ; do
   mv ${out}.${ext} $D/riverh.${ext}
   echo "Created  $D/riverh.${ext}"
done
