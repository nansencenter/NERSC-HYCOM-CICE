#!/bin/bash
# Must be in expt dir to run this script
if [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
elif [ -f EXPT.src ] ; then
   export BASEDIR=$(pwd)
else
   echo "Could not find EXPT.src or REGION.src This script must be run in expt dir or region dir"
   exit 1
fi
export BASEDIR=$BASEDIR/
echo BASEDIR=$BASEDIR

doit=$1
scratchdirs=$(find $BASEDIR -type d -name SCRATCH)
for SCRATCH in $scratchdirs ; do

   if [[ ! "$1" == "doit" ]] ; then
      echo "Candidate for removal $SCRATCH (not removed until 'doit' is  provided as argument to this routine)"
   else 
      echo "Removing files in $SCRATCH "
      find $SCRATCH -type f -exec rm {} \;
   fi
done
