#!/bin/bash
#

# Experiment number needed
if [ $# -ne 3 ] ; then
   echo " $(basename $0) needs experiment version as input (ex 01.0) "
   echo " in addition to synoptic forcing option, start time and end time"
   echo ""
   echo "Example:"
   echo "    $(basename $0) erai 2013-01-01T00:00:00  2013-01-05T00:00:00 "
   echo "Generates forcing from erai over the time span"
   exit
fi
export forcing=$1
export start=$2
export stop=$3

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
# --- S is scratch directory,
# --- D is permanent directory,
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
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
D=$BASEDIR/force/synoptic/
S=$D/SCRATCH
[ ! -d $D ] && mkdir -p $D
[ ! -d $S ] && mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


#
# --- Sanity check on forcing option
#
if [ $forcing == "erai" ] ; then
   xmlfile=$BASEDIR/../input/era-interim.xml
   xmlident="era-interim+lw"
elif [ $forcing == "erai-lw" ] ; then
   xmlfile=$BASEDIR/../input/era-interim.xml
   xmlident="era-interim"
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




