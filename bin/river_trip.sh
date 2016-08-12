#!/bin/bash
#
pget=cp
pput=cp
#set -x

# Experiment  needs experiment number
echo "This script will set up river forcing files used by HYCOM. The forcing"
echo "is based on the ERA40 runoff data set, coupled with the TRIP data base "
echo "for directing the runoff towards the ocean. "
echo 
echo "Example:"
echo "   nersc_trip"
echo 
echo "NB :River discharge radius set to 200 km - modify trip_to_hycom to change this"
echo "NB2:This script might take some time to finish since it calculates river flow "
echo "    for ~50 years. However, the riverflow need only be calculated once, and "
echo "    the resulting trip_era40_clim.nc file will be valid for ALL model grids.  "
echo "    If you store trip_era40_clim.nc somewhere safe you might modify this script "
echo "    so that only the routine trip_to_hycom is called. That routine will use  "
echo "    an existing trip_era40_clim.nc to calculate a hycom river climatology.     "
echo "NB3:This script will overwrite any other river forcing files in force/rivers/E !!"
sleep 6



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
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src            || { echo "Could not source ./EXPT.src" ; exit 1 ; }

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
# --- Create HYCOM river forcing
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
#NB: overwrites any 
D=$BASEDIR/force/rivers/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input.
#
touch regional.grid.a regional.grid.b blkdat.input grid.info regional.depth.a regional.depth
rm    regional.grid.a regional.grid.b blkdat.input grid.info regional.depth.a regional.depth
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a  || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/grid.info grid.info  || { echo "Could not get file rivers.dat " ; exit 1 ; }


# TODO: These two need only be called once if you keep the resulting netcdf files from trip_flow
$MSCPROGS/bin_setup/trip_weights   || { echo "Error when running riverweight (see errors above)" ; exit 1 ; }
$MSCPROGS/bin_setup/trip_flow      || { echo "Error when running riverflow  (see errors above)" ; exit 1 ; }

# This uses trip_era40_clim.nc (from trip_flow) to interpolate to model grid. Uses conformal mapping
$MSCPROGS/bin_setup/trip_tohycom      || { echo "Error when running trip_tohycom  (see errors above)" ; exit 1 ; }


# For now the river forcing is  experiment-dependent
mkdir -p $D/${E}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done
#
echo "Diagnostics can be found in directory $S"
echo "HYCOM river forcing files can be found in $D"


