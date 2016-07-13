#!/bin/bash
# Set up 2d tiling for region

# Input args
if [ $# -ne 2 ]  ; then
#    echo "This script will set up the partition files needed when running HYCOM with"
#    echo "MPI parallelization. The input is the number of partitions along the 1st"
#    echo "and 2nd dimensions, and the topography version to apply the partitioning to"
#    echo
#    echo "Usage:" 
#    echo "   $(basename $0) [1st tile dimension] [2nd tile dimension] [sfudge] [topo version]" 
#    echo "Example:" 
#    echo "   $(basename $0) -2 -2 01 " 
#    echo "  tip  - two negative tile dimensions gives a uniform grid "
#    echo "sfudge : "
#    echo " c ---   sfudge:  size fudge factor (0.5 to 9.9, larger for more variation)"
#    echo " c ---              < 1.0 to keep all  constant-sized tiles"
#    echo " c ---              > 9.0 to shift     constant-sized tiles"
#    echo " c ---              > 9.8 to double-up constant-sized tiles"
   echo "Need topo version and ports file"
    exit 1
fi
#iqr=$1
#jqr=$2
#sfudge=$3
T=$1
PF=$2



# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd $(dirname $0)/.. && pwd)/  
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
export TOPO_PORTS=/home/nersc/knutali/Models/hycom/HYCOM_ALL_2.2.72/ALL/topo/src/topo_ports


# Create work dir, and copy files to it
SUBRDIR=$BASEDIR/subregion/
SCRATCH=$BASEDIR/subregion/SCRATCH/
[ ! -d $SUBRDIR ] && mkdir -p $SUBRDIR
[ ! -d $CSCRATCH ] && mkdir -p $CSCRATCH
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
exit

## Calculate "load-balancing" patches
#echo "$iqr $jqr $sfudge" | $PARTIT
#nm=$(head -n2 fort.21 | tail -n1 | tr -s  " " | cut -d " " -f 2)
##nm=$(printf "%4.4d" $nm)
#nm=`echo 000$nm | tail -5c`
#echo "partitioning finished"

cp fort.21 depth_${R}_${T}.$nm

## Create graphics
##echo 1 | $TOPO_PPM
#$TOPO_PPM
#mv fort.31 depth_${R}_${T}.$nm.ppm
#rm fort.*
#
## Show it if ImageMagick (display) is present...
#echo "Output partition file is  depth_${R}_${T}.$nm, with $nm tiles"
#which display &> /dev/null && display  depth_${R}_${T}.$nm.ppm
#
## Copy partitions to parent directory - overwrite if present
#mv depth_${R}_${T}.$nm.ppm depth_${R}_${T}.$nm $BASEDIR/topo/partit
#
#echo
#echo "partition file depth_${R}_${T}.$nm placed in $BASEDIR/topo/partit"
