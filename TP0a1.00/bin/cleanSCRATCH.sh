#!/bin/bash

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
BASEDIR=`cd $(dirname $0)/.. && pwd`
export BASEDIR=$BASEDIR/
echo BASEDIR=$BASEDIR


#find $BASEDIR -type d -name SCRATCH -exec rm -r {} \;
scratchdirs=$(find $BASEDIR -type d -name SCRATCH)
for SCRATCH in $scratchdirs ; do
   echo "Removing in $SCRATCH"
   find $SCRATCH -type f -exec rm {} \;
done
