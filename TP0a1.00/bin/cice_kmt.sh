#!/bin/bash
# Set up 2d tiling for region

# Input args
if [ $# -ne 1 ]  ; then
    echo
    echo "Usage:" 
    echo "   $(basename $0) topo_version"
    exit 1
fi
T=$1



# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd $(dirname $0)/.. && pwd)/  
source ${BASEDIR}/bin/common_functions.sh  || { echo "Could not source ${BASEDIR}/bin/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
D=$BASEDIR/topo/
S=$D/SCRATCH 
mkdir -p $S
cd $S || { echo "Could not descend dir  $S" ; exit 1 ;}

infile=$(topo_file $R ${T})
kmtfile=$(kmt_file $R $T)
cp ${BASEDIR}/topo/$infile.a $S  || { echo "Could not get $infile.a file " ; exit 1 ; }
cp ${BASEDIR}/topo/$infile.b $S  || { echo "Could not get $infile.b file " ; exit 1 ; }

[ -f cice_kmt.nc ] && rm cice_kmt.nc
cice_kmt.py $infile
cp cice_kmt.nc $D/$kmtfile  || { echo "Could not get cice_kmt.nc file " ; exit 1 ; }
echo "Normal exit. created kmt file $D/$kmtfile"
