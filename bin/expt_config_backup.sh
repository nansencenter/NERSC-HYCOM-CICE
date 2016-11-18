#!/bin/bash
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

msg="
Usage: $(basename $0) [-n] path_to_storage
 
Example: $(basename $0) $HOME/ModelInput/TOPAZ5/TP5a0.06/
 
Purpose: Copies config files for a experiment to a backup location
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


if [[ $norstore -eq 0 ]] ; then
   [ ! -d $targetdir/$exptdir ] && mkdir $targetdir/$exptdir
   echo "Copying config files in $(pwd) to  $targetdir/$exptdir"
   cp * $targetdir/$exptdir
   echo "..."

else

   # Also do a transfer to norstore: /projects/NS2993K/HYCOM-CICE/Setup/
   NORSTORE=norstore.uio.no
   NSUSER=$USER
   NSBASEDIR=/projects/NS2993K/HYCOM-CICE/Setup/$NSUSER/
   echo "hh" $norstore
   if [[ $norstore -eq 1 ]] ; then
      echo "Transferring to backup on $NORSTORE  $USER@$NORSTORE:$NSBASEDIR/$R/expt_$X/"
      echo "ssh $NSUSER@$NORSTORE mkdir -p $NSBASEDIR/$R/expt_$X/ "
      ssh $NSUSER@$NORSTORE mkdir -p $NSBASEDIR/$R/expt_$X/ 
      scp * $NSUSER@$NORSTORE:$NSBASEDIR/$R/expt_$X/ 
   fi
fi

echo "Normal exit"
echo "PS1 :Dont forget to backup region files as well (region_backup.sh)"
echo "PS2 :Dont forget to backup data files as well"
exit 0

