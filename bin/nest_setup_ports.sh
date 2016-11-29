#!/bin/bash
# Set up 2d tiling for region

# Input args
if [ $# -ne 2 ]  ; then
   echo " Script calculates ports and relaxation zones for use when nesting"
    echo "Need  width of relax zone, efold time in days"
    exit 1
fi
width=$1
efold=$2



# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
export TARGETDIR=$BASEDIR/nest/$E
export SCRATCH=$BASEDIR/nest/SCRATCH/
[ ! -d $TARGETDIR ] && mkdir $TARGETDIR
[ ! -d $SCRATCH ] && mkdir $SCRATCH

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
#export TOPO_FIND=/home/nersc/knutali/Models/hycom/HYCOM_ALL_2.2.72/ALL/topo/src/topo_ports_find
export PORT_ROUTINE=hycom_bathy_ports.py

# Create work dir, and copy files to it
echo "running routine $PORT_ROUTINE"
echo "working in  directory $SCRATCH"
cd $SCRATCH || { echo "Could not descend dir  $SCRATCH" ; exit 1 ;}
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b fort.21 fort.31
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b fort.21 fort.31
cp  $BASEDIR/topo/regional.grid.a . || { echo "Could not copy regional.grid.a " ; exit 1 ;}
cp  $BASEDIR/topo/regional.grid.b . || { echo "Could not copy regional.grid.b " ; exit 1 ;}
cp  $BASEDIR/topo/depth_${R}_${T}.a . || { echo "Could not copy depth_${R}_${T}.a " ; exit 1 ;}
cp  $BASEDIR/topo/depth_${R}_${T}.b . || { echo "Could not copy depth_${R}_${T}.b " ; exit 1 ;}
cp  $INPUTDIR/xbathy.pal  . ||  { echo "Could not copy xbathy.pal " ; exit 1 ;}


touch  ports.input.tmp && rm ports.input.tmp
echo ${PORT_ROUTINE}  depth_${R}_${T} $width $efold
${PORT_ROUTINE}  depth_${R}_${T} $width $efold
if [ $? -eq 0 ] ; then
   cp  ports.input.tmp $TARGETDIR/ports.nest
   echo "Created ports file $TARGETDIR/ports.nest"
   cp  rmu.a $TARGETDIR/rmu.a
   cp  rmu.b $TARGETDIR/rmu.b
else 
   echo "port creation routine failed "
   exit 1
fi


## The tiler program needs these
#export FOR051=regional.depth.b
#export FOR051A=regional.depth.a
#export FOR021=fort.21

#${TOPO_FIND}




## Calculate "load-balancing" patches
#echo "$iqr $jqr $sfudge" | $PARTIT
#nm=$(head -n2 fort.21 | tail -n1 | tr -s  " " | cut -d " " -f 2)
##nm=$(printf "%4.4d" $nm)
#nm=`echo 000$nm | tail -5c`
#echo "partitioning finished"
#
#cp fort.21 depth_${R}_${T}.$nm
#
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
