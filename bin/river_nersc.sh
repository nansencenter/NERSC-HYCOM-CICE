#!/bin/bash
#
pget=cp
pput=cp
#set -x

# Experiment  needs experiment number
if [ $# -ne 2 ] ; then
   echo "This script will set up river forcing files used by HYCOM. You will"
   echo "need to have a \"rivers.dat\" file present in the experiment directory."
   echo "Input to this script is experiment number, and normal and alongshore radius."
   echo "The two latter are used to place river discharge, with the majority of "
   echo "the discharge concentrated along land."
   echo 
   echo "Example:"
   echo "   nersc_rivers.sh 100 300"
   echo "Will place rivers at most 100 km from land, and at most 300 km along the land "
   echo "contour."
   exit 1
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
${pget} ${BASEDIR}/expt_${X}/rivers.dat Data/rivers.dat || { echo "Could not get file rivers.dat " ; exit 1 ; }
${pget} ${BASEDIR}/topo/grid.info grid.info  || { echo "Could not get file rivers.dat " ; exit 1 ; }



$MSCPROGS/bin_setup/rivers $1 $2 || { echo "Error when running rivers (see errors above)" ; exit 1 ; }

# For now the river forcing is  experiment-dependent
mkdir -p $D/${E}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done

echo "Diagnostics can be found in directory $S"
echo "HYCOM river forcing files can be found in $D"


