#!/bin/bash
#
pget=cp
pput=cp
#set -x

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

# Need TRIP_PATH to be set
if [ -z "${TRIP_PATH}" ] ; then
   echo "TRIP_PATH is not set. Set that first (pref in REGION.src)"
   exit 1
fi

# Experiment  needs experiment number
along=200
across=60
create=0
usage="
This script will set up river forcing files used by HYCOM. The forcing
is based on the ERA40 runoff data set, coupled with the TRIP data base 
for directing the runoff towards the ocean. 

In addition to create a river climatology (forcing.rivers.[ab], this routine
will also create precipitation fields that represents realtime river flow as described
bvy the TRIP routines.

Note that the precip fields will have to be further manipulated if they are to be used by 
hycom (they will need to include precipitation fields as well, and they need to 
be copied in by the preprocessing script).


Example:
   $(basename $0) [ -t along_shore_discharge_radius ] [ -n across_shore_discharge_radius] [-c] runoff_source time1 time2 deltatime

Argument
   runoff_source                     : specify runoff source, erai or era40 
   time1                             : First time to process (format YYYY-mm-ddTHH:MM:ss)
   time1                             : Last  time to process (format YYYY-mm-ddTHH:MM:ss)
   deltat                            : Time step (in days)

Optional arguments:
   -t along_shore_discharge_radius  : rivers will be spread along shore using this radius (default $along)
   -n across_shore_discharge_radius : rivers will be spread normal to the shore using this radius (default $across)
   -c                               : Create TRIP weights and flow in local scratch dir. Otherwise use existing data. Defult: use existing

   
NB:This script will overwrite any other river forcing files in force/rivers/E !!
"


# This will process optional arguments
options=$(getopt -o ct:n: -- "$@")
[ $? -eq 0 ] || {
    echo "$usage"
    echo "Error: Incorrect options provided"
    exit 1
}
eval set -- "$options"
while true; do
    case "$1" in
    -t)
       shift;
       along=$1
        ;;
    -n)
       shift;
       across=$1
        ;;
    -c)
       create=1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done

if [ $# -ne 4 ] ; then
   echo -e "$usage"
   echo "No source specified!"
   exit 1;
fi
src=$1
starttime=$2
endtime=$3
deltat=$4


# Convert times to hycom dtime
dtime1=$(hycom_date.py datetime dtime $starttime)
if [ $? -ne 0 ] ; then
   echo "$usage"
   echo "An error occured when converting from datetime=$starttime to hycom time"
   exit 1
fi

dtime2=$(hycom_date.py datetime dtime $endtime)
if [ $? -ne 0 ] ; then
   echo "$usage"
   echo "An error occured when converting from datetime=$endtime to hycom time"
   exit 1
fi
echo $dtime1
echo $dtime2


flowpath="${TRIP_PATH}/$src"
flowfile="trip_${src}_clim.nc"
echo "Along-shore  radius               :$along"
echo "Across-shore radius               :$across"
echo "Runoff source                     :${src}"
echo "Create riverweights and riverflow :${create} (0=False)"
if [ $create -eq 0 ] ; then
   echo "Location of existing weight and flow data :${flowpath}"
   echo "Flow file                                 :${flowpath}/${flowfile}"
fi
rwcellinfo="rw_cellinfo.uf"
rwmaxncell="rw_maxncells.asc"




#
#
# --- Create HYCOM river forcing
#
#
# --- S is scratch directory,
# --- D is target directory,
#
#NB: overwrites any 
D=$BASEDIR/force/rivers/$E/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}
#
# --- Input.
#
touch regional.grid.a regional.grid.b blkdat.input grid.info regional.depth.a regional.depth $rwcellinfo $rwmaxncell
rm    regional.grid.a regional.grid.b blkdat.input grid.info regional.depth.a regional.depth $rwcellinfo $rwmaxncell
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a  || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/grid.info grid.info  || { echo "Could not get file rivers.dat " ; exit 1 ; }


# Create new files
if [ $create -ne 0 ] ; then
   echo
   echo "**Calculating TRIP weights"
   $MSCPROGS/bin_setup/trip_weights $src  || { echo "Error when running riverweight (see errors above)" ; exit 1 ; }
   echo
   echo "**Calculating TRIP river flow"
   $MSCPROGS/bin_setup/trip_flow    $src  || { echo "Error when running riverflow  (see errors above)" ; exit 1 ; }
# Copy from previously stored files 
else 
   cp $flowpath/$flowfile   .  || { echo "Could not get file $flowpath/$flowfile"   ; exit 1 ; }
   cp $flowpath/$rwcellinfo .  || { echo "Could not get file $flowpath/$rwcellinfo" ; exit 1 ; }
   cp $flowpath/$rwmaxncell .  || { echo "Could not get file $flowpath/$rwmaxncell" ; exit 1 ; }
fi

# This uses trip_era40_clim.nc (from trip_flow) to interpolate to model grid. Uses conformal mapping
echo
echo "**Interpolating TRIP river flow to hycom"
#$MSCPROGS/bin_setup/trip_tohycom $src $along $across    || { echo "Error when running trip_tohycom  (see errors above)" ; exit 1 ; }
$MSCPROGS/bin_setup/trip_tohycomrt $src $along $across $dtime1 $dtime2 $deltat   || { echo "Error when running trip_tohycomrt  (see errors above)" ; exit 1 ; }

# TODO: Its possible to use the river flow data and create fake precipitation fields. This can be used as realtime river forcing 


for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D
done
#
echo "Diagnostics can be found in directory $S"
echo "HYCOM river forcing files can be found in $D"


