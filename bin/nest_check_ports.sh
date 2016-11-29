#!/bin/bash
# Set up 2d tiling for region

# Input args
if [ $# -ne 2 ]  ; then
   echo "Does a consistency check on the ports file set up for this experiment"
   echo "Check is done using HYCOM_ALL routine topo/topo_ports"
   echo "Need topo version and ports file"
   echo
   echo "Checks a port for consistency and writes to stdout in case of errors in port config"
   exit 1
fi
#iqr=$1
#jqr=$2
#sfudge=$3
T=$1
PF=$2



# Must be in expt dir to run this script
if [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
elif [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find REGION.src. This script must be run in expt dir or in region dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }

# Check that pointer to HYCOM_ALL is set (from EXPT.src)
if [ -z ${HYCOM_ALL} ] ; then
   echo "Environment not set "
   exit
else
   if [ ! -d ${HYCOM_ALL} ] ; then
      echo "HYCOM_ALL not properly set up"
      echo "HYCOM_ALL not a directory at ${HYCOM_ALL}"
      exit
   fi
fi

# Pointers to programs
export TOPO_PORTS=${HYCOM_ALL}/topo/src/topo_ports


# Create work dir, and copy files to it
SUBRDIR=$BASEDIR/subregion/
SCRATCH=$BASEDIR/subregion/SCRATCH/
[ ! -d $SUBRDIR ] && mkdir -p $SUBRDIR
[ ! -d $SCRATCH ] && mkdir -p $SCRATCH
cp  $PF $SCRATCH/ports.input
cd $SCRATCH || { echo "Could not descend dir  $SCRATCH" ; exit 1 ;}
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b fort.21 fort.31
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b fort.21 fort.31
cp  $BASEDIR/topo/regional.grid.a . || { echo "Could not copy regional.grid.a " ; exit 1 ;}
cp  $BASEDIR/topo/regional.grid.b . || { echo "Could not copy regional.grid.b " ; exit 1 ;}
cp  $BASEDIR/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not copy depth_${R}_${T}.a " ; exit 1 ;}
cp  $BASEDIR/topo/depth_${R}_${T}.b regional.depth.b || { echo "Could not copy depth_${R}_${T}.b " ; exit 1 ;}


# The tiler program needs these
export FOR051=regional.depth.b
export FOR051A=regional.depth.a
export FOR021=fort.21

${TOPO_PORTS}
exit 0
