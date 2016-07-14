#!/bin/bash
#

# Experiment number needed
if [ $# -ne 4 ] ; then
   echo " $(basename $0) needs experiment version as input (ex 01.0) "
   echo " in addition to synoptic forcing option, start time and end time"
   echo ""
   echo "Example:"
   echo "    $(basename $0) 01.0 erai 2013-01-01T00:00:00  2013-01-05T00:00:00 "
   echo "Generates forcing from erai over the time span"
   exit
fi
export X=$1
export forcing=$2
export start=$3
export stop=$4

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
# --- S is scratch directory,
# --- D is permanent directory,
export BASEDIR=$(cd  $(dirname $0)/../ && pwd)/
source ${BASEDIR}/bin/common_functions.sh || { echo "Could not source ${BASEDIR}/bin/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }
D=$BASEDIR/force/synoptic/
S=$D/SCRATCH
[ ! -d $D ] && mkdir -p $D
[ ! -d $S ] && mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


#
# --- Sanity check on forcing option
#
if [ $forcing == "erai" ] ; then
   xmlfile=$BASEDIR/force/input/era-interim.xml
   xmlident="era-interim+lw"
else 
   tellerror "Forcing option is erai only..."
   exit 1
fi


#
# --- Input. Function in common_functions.sh
#
copy_setup_files $S





cmd="$BASEDIR/../python/hycom_atmfor.py $start $stop $xmlfile $xmlident"
eval $cmd   ||  { echo "Error running $cmd " ; exit 1 ; }

# The nersc era40 forcing is region-independent 
mkdir -p $D/$E  || { echo " Could not create data  destination dir $D/$E" ; exit 1;}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
   echo "Created  $D/$E/$new"
done




