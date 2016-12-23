#!/bin/bash

usage="
This routine will use FES to calculate the didal forcing on the boundary of your
model domain. The tidal forcing will then be used to calculate hycom archive
files suitable for use as nesting files.

There are two ways to call this routine:

1)Input to this routine is a bunch of archive files. The script will extract the
times of the archive files, provide tidal forcing for those, and add to the
archive files.

2) Input to this routine is a start time, an end time, and a time step in days.
Tidal forcing will be generated from start time to end time with the provided
time step. Any archive files given on the command line will be modified to take
the tidal forcing into account IF the times match those of the tidal data



Usage 1:
   $(basename $0) [--include-uv] [--tidal-database dbname] 1 archive_file1 [archive_file2 ...]

Usage 2:
   $(basename $0) [--include-uv] [--tidal-database dbname] 2 start_time end_time delta_t [archive_file1 archive_file2 ..]


Arguments:
   archive_file : a HYCOM archv file
   tidal_dataset: FES2004 or FES2014
   start_time   : (Only usage 2) Start of tide generation. Format YYYY-mm-ddTHH:mm:ss
   end_time     : (Only usage 2) End of tide generation.   Format YYYY-mm-ddTHH:mm:ss
   delta_t      : (Only usage 2) Time step in hours

Optional arguments
   --tidal-database dbase  : Tidal database to use; FES2004 or FES2014(default)
   --include-uv            : Include tidal currents (only FES2014)


Examples:
   $(basename $0) 1 2013-01-01T00:00:00 2013-01-10T00:00:00 1.0 
   $(basename $0) 1 2013-01-01T00:00:00 2013-01-10T00:00:00 1.0  archv.2013_003_12.a archv.2013_004.a 
   $(basename $0) 2 archv.2013_003_12.a archv.2013_004.a 
"


#Set environment
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
export STARTDIR=$(pwd)
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }


if [ $#  -lt 2 ] ; then
   echo "$usage"
   exit 1
fi


# Create SCRATCH  dir
SUBDIR=$BASEDIR/nest/${E}_with_tides/
SCRATCH=$SUBDIR/SCRATCH/
[ ! -d $SUBDIR ] && mkdir -p $SUBDIR
[ ! -d $SCRATCH ] && mkdir -p $SCRATCH
cd $SCRATCH || { echo "Could not descend dir  $SCRATCH" ; exit 1 ;}

# remove old archove files lying around
rm $SUBDIR/archv*.[ab]

# Copy in  new files
touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b
cp $BASEDIR/topo/regional.grid.a . || { echo "Could not copy regional.grid.a " ; exit 1 ;}
cp $BASEDIR/topo/regional.grid.b . || { echo "Could not copy regional.grid.b " ; exit 1 ;}
cp $BASEDIR/topo/depth_${R}_${T}.a regional.depth.a || { echo "Could not copy depth_${R}_${T}.a " ; exit 1 ;}
cp $BASEDIR/topo/depth_${R}_${T}.b regional.depth.b || { echo "Could not copy depth_${R}_${T}.b " ; exit 1 ;}
cp $P/blkdat.input . || { echo "Could not copy blkdat.input " ; exit 1 ;}



# Let python take over remaining logic. archv_basedir is set in case any archive files provided use relative paths
echo "Running command in $SCRATCH:  tidal_forcing.py $* "
tidal_forcing.py --archv-basedir=$P $*
rc=$?
if [ $rc -ne 0 ] ; then
   echo "tidal_forcing.py failed - see above"
   exit 1
fi

