#!/bin/bash
#
#set -x
pget=cp
pput=cp

# Experiment number needed
if [ $# -ne 3 ] ; then
   echo " $(basename $0) needs experiment version as input (ex 01.0) "
   echo " in addition to synoptic and climate forcing options"
   echo "Example:"
   echo "$(basename $0) 01.0 era40 era40"
   exit
fi
export X=$1

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd  $(dirname $0)/../ && pwd)/
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }
# --- S is scratch directory,
# --- D is permanent directory,
D=$BASEDIR/force/synoptic/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


# Check that pointer to MSCPROGS is set
if [ -z ${MSCPROGS} ] ; then
   echo "MSCPROGS Environment not set "
   exit
else
   if [ ! -d ${MSCPROGS} ] ; then
      echo "MSCPROGS not properly set up"
      echo "MSCPROGS not a directory at ${MSCPROGS}"
      exit
   fi
fi



#
# --- Input.
#
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b infile.in blkdat.input
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b infile.in blkdat.input
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b     || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a     || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }
${pget} ${BASEDIR}/expt_${X}/blkdat.input blkdat.input      || { echo "Could not get file blkdat.input " ; exit 1 ; }
${pget} ${BASEDIR}/expt_${X}/infile.in infile.in            || { echo "Could not get file infile.in " ; exit 1 ; }


$MSCPROGS/bin_setup/forfun_nersc-2.2.37 $2 $3 ||  { echo "Error running forfun_nersc-2.2.37 " ; exit 1 ; }

# The nersc era40 forcing is region-independent 
mkdir -p $D/$E  || { echo " Could not create data  destination dir $D/$E" ; exit 1;}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done




