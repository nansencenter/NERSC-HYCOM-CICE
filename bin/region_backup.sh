#!/bin/bash

msg="
Usage: $(basename $0) path_to_storage
 
Example: $(basename $0) $HOME/ModelInput/TOPAZ5/TP5a0.06/
 
Purpose: Copies config files for a region (Any files in region dir + topo files) to a backup location
"



set -e
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
elif [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
else
   echo "Could not find REGION.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
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

[ ! -d $targetdir/topo ] && mkdir $targetdir/topo

echo "Copying topo files in $BASEDIR/topo to  $targetdir/topo"
cp -r $BASEDIR/topo/* $targetdir/topo
echo "Copying files in $BASEDIR/ to  $targetdir/"
cp $BASEDIR/* $targetdir/

echo "Normal exit"
echo "PS1 :Dont forget to backup config files as well (expt_backup_config.sh)"
echo "PS2 :Dont forget to backup data files as well"
exit 0

