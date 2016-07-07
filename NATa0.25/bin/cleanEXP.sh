#!/bin/bash

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
cd $(dirname $0)/../
export BASEDIR=$(pwd)/  
cd -


if [ $# -ne 1 ] ; then
   echo "Inut is experiment number E"
   exit 1
fi


if ! echo $1 | egrep ^[[:digit:]]{3}$ &> /dev/null ; then
   echo "Experiment E must be a text string consisting of three numbers"
   exit 1
fi



find $BASEDIR -type d -name $1 -exec rm -r {} \;
