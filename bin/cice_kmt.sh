#!/bin/bash
# Set up 2d tiling for region

# Input args
if [ $# -ne 1 ]  ; then
   echo "Create cice mask file from hycom bathymetry file"
    echo
    echo "Usage:" 
    echo "   $(basename $0) topo_version"
    exit 1
fi
T=$1


# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
elif [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
else
   echo "Could not find EXPT.src or REGION.src . This script must be run in expt dir or region dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
source ${BINDIR}/common_functions.sh  || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
D=$BASEDIR/topo/
S=$D/SCRATCH 
mkdir -p $S
cd $S || { echo "Could not descend dir  $S" ; exit 1 ;}

infile=$(topo_file $R ${T})
kmtfile=$(kmt_file $R $T)
cp ${BASEDIR}/topo/$infile.a $S  || { echo "Could not get $infile.a file " ; exit 1 ; }
cp ${BASEDIR}/topo/$infile.b $S  || { echo "Could not get $infile.b file " ; exit 1 ; }
cp ${BASEDIR}/topo/regional.grid.b $S  || { echo "Could not get regional.grid.b file " ; exit 1 ; }
cp ${BASEDIR}/topo/regional.grid.a $S  || { echo "Could not get regional.grid.a file " ; exit 1 ; }

echo "Running routine cice_kmt.py in directory $S"
[ -f cice_kmt.nc ] && rm cice_kmt.nc
#$BASEDIR/../python/cice_kmt.py $infile
cice_kmt.py $infile
cp cice_kmt.nc $D/$kmtfile  || { echo "Could not get cice_kmt.nc file " ; exit 1 ; }
echo "Normal exit. created kmt file $D/$kmtfile"
