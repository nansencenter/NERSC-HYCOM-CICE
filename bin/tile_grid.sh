#!/bin/bash
# Set up 2d tiling for region
sfudge="9.5"
usage="This script will set up the partition files needed when running HYCOM with
MPI parallelization. The input is the number of partitions along the 1st
and 2nd dimensions, and the topography version to apply the partitioning to

The fudge factor can be specified as optional argument, default value is $sfudge 
See also text about sfudge below

Usage:
   $(basename $0) [-s sfudge] 1st_tile_dimension 2nd_tile_dimension topo_version
Example:
   $(basename $0) -2 -2 01 
  tip  - two negative tile dimensions gives a uniform grid
from partit routine: 
 c ---   sfudge:  size fudge factor (0.5 to 9.9, larger for more variation)
 c ---              < 1.0 to keep all  constant-sized tiles
 c ---              > 9.0 to shift     constant-sized tiles
 c ---              > 9.8 to double-up constant-sized tiles "

if [ "$1" == "-s" ]   ; then
   sfudge=$2
   shift 2
fi
echo "size fudge factor: $sfudge"


# Input args
if [ $# -ne 3 ]  ; then
    echo "$usage"
    exit 1
fi
iqr=$1
jqr=$2
T=$3

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
elif [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
else
   echo "Could not find EXPT.src or REGION.src. This script must be run in expt dir or region dir"
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
export PARTIT=${HYCOM_ALL}/topo/src/partit
export TOPO_PPM=${HYCOM_ALL}/topo/src/topo_ppm


# Create work dir, and copy files to it
mkdir -p $BASEDIR/topo/partit/SCRATCH 
cd $BASEDIR/topo/partit/SCRATCH || \
   { echo "Could not descend dir  $BASEDIR/topo/partit/SCRATCH" ; exit 1 ;}
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b fort.21 fort.31
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b fort.21 fort.31
cp  $BASEDIR/topo/regional.grid.a . || { echo "Could not copy regional.grid.a " ; exit 1 ;}
cp  $BASEDIR/topo/regional.grid.b . || { echo "Could not copy regional.grid.b " ; exit 1 ;}
cp  $BASEDIR/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not copy depth_${R}_${T}.a " ; exit 1 ;}
cp  $BASEDIR/topo/depth_${R}_${T}.b regional.depth.b || { echo "Could not copy depth_${R}_${T}.b " ; exit 1 ;}
cp  $INPUTDIR/xbathy.pal  . ||  { echo "Could not copy xbathy.pal " ; exit 1 ;}

# The tiler program needs these
export FOR051=regional.depth.b
export FOR051A=regional.depth.a
export FOR021=fort.21

# Calculate "load-balancing" patches
echo "$iqr $jqr $sfudge" | $PARTIT
nm=$(head -n2 fort.21 | tail -n1 | tr -s  " " | cut -d " " -f 2)
#nm=$(printf "%4.4d" $nm)
nm=`echo 000$nm | tail -5c`
echo "partitioning finished"

cp fort.21 depth_${R}_${T}.$nm

# Create graphics
#echo 1 | $TOPO_PPM
$TOPO_PPM
mv fort.31 depth_${R}_${T}.$nm.ppm
rm fort.*

# Show it if ImageMagick (display) is present...
echo "Output partition file is  depth_${R}_${T}.$nm, with $nm tiles"
#which display &> /dev/null && display  depth_${R}_${T}.$nm.ppm 

# Copy partitions to parent directory - overwrite if present
mv depth_${R}_${T}.$nm.ppm depth_${R}_${T}.$nm $BASEDIR/topo/partit

echo
echo "Created partition file    :  $BASEDIR/topo/partit/depth_${R}_${T}.$nm"
echo "Created ppm graphics file :  $BASEDIR/topo/partit/depth_${R}_${T}.$nm.ppm"
