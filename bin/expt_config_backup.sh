#!/bin/bash

msg="
Usage: $(basename $0) path_to_storage
 
Example: $(basename $0) $HOME/ModelInput/TOPAZ5/TP5a0.06/
 
Purpose: Copies config files for a experiment to a backup location
"



set -e
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src "; exit 1; }
source $BINDIR/common_functions.sh  || { echo "Could not source common_functions.sh "; exit 1; }


if [ $# -ne 1 ] ; then
   echo -e "$msg"
   tellerror "Need target dir as input"
   exit 1
fi

targetdir=$1
if [ ! -d $targetdir ] ; then
   echo -e $msg
   tellerror "target dir must be an existing directory"
   exit 1
fi
exptdir=$(basename $(pwd))

[ ! -d $targetdir/$exptdir ] && mkdir $targetdir/$exptdir

echo "Copying config files in $(pwd) to  $targetdir/$exptdir"
cp * $targetdir/$exptdir

echo "Normal exit"
echo "PS1 :Dont forget to backup region files as well (region_backup.sh)"
echo "PS2 :Dont forget to backup data files as well"
exit 0

