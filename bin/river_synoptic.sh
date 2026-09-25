#!/bin/bash
#

# Experiment number needed
if [ $# -ne 2 ] ; then
   echo ""
   echo "Generates river forcing over the time span"
   echo ""
   echo " $(basename $0) needs synoptic river forcing option, start time and end time"
   echo ""
   echo "Example:"
   echo "    $(basename $0) 2013-01-01T00:00:00  2013-01-05T00:00:00"
   echo ""
   exit
fi
export start=$1
export stop=$2

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

D=$BASEDIR/force/rivers/$E
S=$D/SCRATCH
[ ! -d $D ] && mkdir -p $D
[ ! -d $S ] && mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input. Function in common_functions.sh
#
#copy_setup_files $S

# check the right depth files:
if [ -s "${BASEDIR}/topo/regional.depth.a" ] && [ -s "${BASEDIR}/topo/regional.depth.b" ]; then
   depthfile="${BASEDIR}/topo/regional.depth.a"
elif [ -s "${BASEDIR}/topo/depth_${R}_${T}.a" ] && [ -s "${BASEDIR}/topo/depth_${R}_${T}.b" ]; then
   depthfile="${BASEDIR}/topo/depth_${R}_${T}.a"
else
   echo "Could not find depth files in ${BASEDIR}/topo/: regional.depth.[ab] or depth_${R}_${T}.[ab]"
   exit 1
fi

ml load matplotlib/3.5.2-foss-2022a
cd ${BINDIR}/river-topaz/scripts/GloFAS

echo "python ${BINDIR}/river-topaz/scripts/GloFAS/hycom_river_hifre.py $start $stop $S/ $depthfile"

cmd="python ${BINDIR}/river-topaz/scripts/GloFAS/hycom_river_hifre.py $start $stop $S/ $depthfile"
eval $cmd   ||  { echo "Error running $cmd " ; exit 1 ; }

# The nersc river forcing is region-independent 
cd $S
echo ""
for i in riverh_${start:0:10}_${stop:0:10}_*.[ab] ; do
   echo ""
   new=riverh${i#*6.0h}
   [ -s $D/$new ] && rm $D/$new
   mv $i $D/$new
   echo "Created  $D/$new"
done




