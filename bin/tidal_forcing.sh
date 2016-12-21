#!/bin/bash

#Set environment
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
export STARTDIR=$(pwd)
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }



# Create SCRATCH  dir
SUBDIR=$BASEDIR/nest/${E}_with_tides/
SCRATCH=$SUBDIR/SCRATCH/
[ ! -d $SUBDIR ] && mkdir -p $SUBDIR
[ ! -d $SCRATCH ] && mkdir -p $SCRATCH
cd $SCRATCH || { echo "Could not descend dir  $SCRATCH" ; exit 1 ;}

# remove old archove files lying around
rm $SUBDIR/archv*.[ab]

# Copy in  new files
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b
cp $BASEDIR/topo/regional.grid.a . || { echo "Could not copy regional.grid.a " ; exit 1 ;}
cp $BASEDIR/topo/regional.grid.b . || { echo "Could not copy regional.grid.b " ; exit 1 ;}
cp $BASEDIR/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not copy depth_${R}_${T}.a " ; exit 1 ;}
cp $BASEDIR/topo/depth_${R}_${T}.b regional.depth.b || { echo "Could not copy depth_${R}_${T}.b " ; exit 1 ;}
cp $P/blkdat.input . || { echo "Could not copy blkdat.input " ; exit 1 ;}


# Let python take over remaining logic. archv_basedir is set in case any archive files provided use relative paths
echo "Running command in $SCRATCH:  tidal_forcing.py $* "
tidal_forcing.py --archv-basedir=$P $*
rc=$?
if [ $rc -ne 0 ] ; then
   echo "tidal_forcing.py failed - see above"
   exit 1
fi

