#!/bin/bash
#
#
# --- Create HYCOM atmospheric forcing (era40 climatology)
#
#
#
#set -x
pget=cp
pput=cp

# Experiment number needed
if [ $# -ne 4 ] ; then
   echo " $(basename $0) needs experiment version as input (ex 01.0) "
   echo " in addition to synoptic and climate forcing options"
   echo "Example:"
   echo "$(basename $0) 01.0 erai starttime endtime"
   exit
fi
export X=$1
starttime=$3
frcopt=$2
endtime=$4

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd $(dirname $0)/../ && pwd)/
source ${BASEDIR}/bin/common_functions.sh || { echo "Could not source ${BASEDIR}/bin/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }
# --- S is scratch directory,
# --- D is permanent directory,
D=$BASEDIR/force/synoptic/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input. See bin/common_functions.sh
#
copy_setup_files $S


# Source necessary setup for python libs
. ~/Repos/Git/knutalnersc/modeltools/set_env.sh
if [ "$frcopt" == "erai" ] ; then
   ~/Repos/Git/knutalnersc/modeltools/scripts/hycom_atmfor.py  $starttime $endtime /home/nersc/knutali/Repos/Git/knutalnersc/modeltools/input/era-interim.xml era-interim 

   [ $? -ne 0 ] && { echo "hycom_atmfor.py failed" ; exit 1;}
else 
   echo "Could not get forcing for option $frcopt"
   exit 1
fi


# The nersc era40 forcing is region-independent 
mkdir -p $D/$E  || { echo " Could not create data  destination dir $D/$E" ; exit 1;}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done




