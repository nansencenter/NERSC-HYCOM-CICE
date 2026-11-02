#!/bin/bash
#

# Experiment number needed
if [ $# -lt 2 -o $# -gt 3 ] ; then
   echo ""
   echo "Generates GloFAS river forcing (riverh.[ab]) over the time span"
   echo ""
   echo " $(basename $0) needs start time, end time and optionally the record interval dt"
   echo " (fraction of a day). If dt is omitted it is set from the HYCOM version V in EXPT.src:"
   echo "    V=2.3    : dt=1.0  (daily records; priver=2 interpolates linearly between records)"
   echo "    otherwise: dt=0.25 (6-hourly records; highfq_river in 2.2.98 must match the"
   echo "                        atmospheric forcing times)"
   echo ""
   echo "Example:"
   echo "    $(basename $0) 2013-01-01T00:00:00  2013-01-05T00:00:00"
   echo "    $(basename $0) 2013-01-01T00:00:00  2013-01-05T00:00:00 0.25"
   echo ""
   exit
fi
export start=$1
export stop=$2
dt=$3

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
# --- S is scratch directory,
# --- D is permanent directory,
# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }

if [ -z "$dt" ] ; then
   if [ "$V" == "2.3" ] ; then
      dt=1.0
   else
      dt=0.25
   fi
fi
echo "River forcing record interval dt=$dt day (HYCOM version V=$V)"

D=$BASEDIR/force/rivers/$E
S=$D/SCRATCH
[ ! -d $D ] && mkdir -p $D
[ ! -d $S ] && mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input. Function in common_functions.sh
#
#copy_setup_files $S

# Use the experiment's depth file (topography version T in EXPT.src), the same file the
# model run and the climatological rivers (copy_topo_files) use, so the land masks match.
# topo/regional.depth.[ab] may point to another version and is not used.
depthfile="${BASEDIR}/topo/depth_${R}_${T}.a"
if [ ! -s "$depthfile" ] || [ ! -s "${depthfile%.a}.b" ]; then
   echo "Could not find depth files ${BASEDIR}/topo/depth_${R}_${T}.[ab]"
   exit 1
fi

# Needs the hycom-cice Python environment (xarray etc.). Do not load the matplotlib module here:
# it puts the module's Python 3.10 numpy on PYTHONPATH, which breaks the environment's python.
cd ${BINDIR}/river-topaz/scripts/GloFAS

cmd="python ${BINDIR}/river-topaz/scripts/GloFAS/hycom_river_hifre.py --dt $dt $start $stop $S/ $depthfile"
echo "$cmd"
eval $cmd   ||  { echo "Error running $cmd " ; exit 1 ; }

# The nersc river forcing is region-independent 
cd $S
echo ""
for i in riverh_${start:0:10}_${stop:0:10}_*.[ab] ; do
   echo ""
   new=riverh.${i##*.}
   [ -s $D/$new ] && rm $D/$new
   mv $i $D/$new
   echo "Created  $D/$new"
done




