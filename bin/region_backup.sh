#!/bin/bash

msg="
Usage: $(basename $0) [-n] path_to_storage
 
Example: $(basename $0) $HOME/ModelInput/TOPAZ5/TP5a0.06/
 
Purpose: Copies config files for a region (Any files in region dir + topo files) to a backup location
"


options=$(getopt -o n -- "$@")
norstore=0
eval set -- "$options"
while true; do
    case "$1" in
    -n)
       norstore=1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done


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
source $BASEDIR/REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source $BINDIR/common_functions.sh  || { echo "Could not source common_functions.sh "; exit 1; }


if [ $# -ne 1 ] ; then
   echo -e "$msg"
   tellerror "Need target dir as input"
   exit 1
fi
targetdir=$1


if [[ $norstore -eq 0 ]] ; then 
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
else

   # Else do a transfer to norstore: /projects/NS2993K/HYCOM-CICE/Setup/$USER
   NORSTORE=norstore.uio.no
   NSUSER=$USER
   NSBASEDIR=/projects/NS2993K/HYCOM-CICE/Setup/$NSUSER/
   if [[ $norstore -eq 1 ]] ; then
      echo "Transferring to backup on $NORSTORE  $USER@$NORSTORE:$NSBASEDIR/$R/topo/"
      echo "ssh $NSUSER@$NORSTORE mkdir -p $NSBASEDIR/$R/topo/ "
      ssh $NSUSER@$NORSTORE mkdir -p $NSBASEDIR/$R/topo/ 
      scp * $NSUSER@$NORSTORE:$NSBASEDIR/$R/topo/ 
   fi
fi

echo "Normal exit"
echo "PS1 :Dont forget to backup config files as well (expt_backup_config.sh)"
echo "PS2 :Dont forget to backup data files as well"
exit 0

