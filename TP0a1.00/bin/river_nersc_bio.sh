#!/bin/bash
#
pget=cp
pput=cp
#set -x

# Experiment  needs experiment number
if [ $# -ne 3 ] ; then
   echo "This script will set up river nutrient forcing files used by HYCOM. You will"
   echo "need to have a \"biorivers.dat\" file present in the experiment directory."
   echo "Input to this script is experiment number, and normal and alongshore radius."
   echo "The two latter are used to place river nutrient discharge, with the majority of "
   echo "the discharge concentrated along land."
   echo 
   echo "Example:"
   echo "   nersc_biorivers.sh 01.0 100 300"
   echo "Will place rivers at most 100 km from land, and at most 300 km along the land "
   echo "contour. This will be done for experiment 01.0."
   exit
fi
export X=$1


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd $(dirname $0)/../ && pwd)
source ${BASEDIR}/bin/common_functions.sh || { echo "Could not source ${BASEDIR}/bin/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }

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
D=$BASEDIR/force/rivers/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input.
#
touch regional.grid.a regional.grid.b blkdat.input infile.in regional.depth.a regional.depth.b
rm    regional.grid.a regional.grid.b blkdat.input infile.in regional.depth.a regional.depth.b
copy_grid_files $S
copy_topo_files $S
${pget} ${BASEDIR}/expt_${X}/blkdat.input blkdat.input  || { echo "Could not get file blkdat.input " ; exit 1 ; }
${pget} ${BASEDIR}/topo/grid.info grid.info  || { echo "Could not get file grid.info " ; exit 1 ; }
${pget} ${BASEDIR}/expt_${X}/biorivers.dat Data/biorivers.dat  || { echo "Could not get file biorivers.dat " ; exit 1 ; }


# Run program
$MSCPROGS/bin_setup/biorivers $2 $3 || { echo "Error when running biorivers (see errors above)" ; exit 1 ; }

# For now the river forcing is  experiment-dependent
mkdir -p $D/${E}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done

echo "Diagnostics can be found in directory $S"
echo "HYCOM river forcing files can be found in $D"


