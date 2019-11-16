#!/bin/bash
#
pget=cp
pput=cp
#set -x

# Experiment  needs experiment number
if [ $# -ne 3 ] ; then
 if [ $# -ne 4 ] ; then  
   echo 
   echo "Usage:"
   echo "   $(basename $0) normal_radius alongshore_radius rivers.dat-file"
   echo "   $(basename $0) normal_radius alongshore_radius rivers.dat-file biorivers.dat-file"
   echo 
   echo "This script will set up river forcing files used by HYCOM. You will"
   echo "need to have a \"rivers.dat\" file. There should be one present in GITROOT/input."
   echo "Input to this script is normal and alongshore radius and rivers.dat file."
   echo "The two first are used to place river discharge, along and normal to the coase. "
   echo 
   echo "Example:"
   echo "   $(basename $0)h 100 300 ../../inputdir/rivers.dat"
   echo "If biology river forcing needs to be activated:"
   echo "   $(basename $0)h 100 300 ../../inputdir/rivers.dat ../../inputdir/biorivers.dat"
   echo "Will place rivers at most 100 km from land, and at most 300 km along the land "
   echo "contour."
   echo "biorivers.dat is created by glodapBIO_to_hycom.py" 
   exit 1
 else
   echo "river biology forcing is active"
 fi
fi

# Full path 
riverfile=$(readlink -f $3)
if [ "$#" = 4 ]; then
   bio=true
   bioriverfile=$(readlink -f $4)
fi
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
export BASEDIR=$(cd .. && pwd)/  
source ${BASEDIR}/REGION.src         || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
source ./EXPT.src                    || { echo "Could not source ./EXPT.src" ; exit 1 ; }
D=$BASEDIR/force/rivers/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data

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
#
# --- Create HYCOM atmospheric forcing (era40 climatology)
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input.
#
touch blkdat.input infile.in
rm    blkdat.input infile.in 
copy_topo_files $S
copy_grid_files $S
${pget} ${BASEDIR}/expt_${X}/blkdat.input blkdat.input  || { echo "Could not get file blkdat.input " ; exit 1 ; }
${pget} $riverfile Data/rivers.dat || { echo "Could not get file rivers.dat " ; exit 1 ; }
if [ "$bio" = true ]; then
   ${pget} $bioriverfile Data/biorivers.dat || { echo "Could not get file biorivers.dat " ; exit 1 ; }
fi
${pget} ${BASEDIR}/topo/grid.info grid.info  || { echo "Could not get file grid.info " ; exit 1 ; }



$MSCPROGS/bin_setup/rivers $1 $2 || { echo "Error when running rivers (see errors above)" ; exit 1 ; }
if [ "$bio" = true ]; then
    $MSCPROGS/bin_setup/biorivers $1 $2 || { echo "Error when running biorivers (see errors above)" ; exit 1 ; }
fi

# For now the river forcing is  experiment-dependent
mkdir -p $D/${E}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done

echo "Diagnostics can be found in directory $S"
echo "HYCOM river forcing files can be found in $D"


