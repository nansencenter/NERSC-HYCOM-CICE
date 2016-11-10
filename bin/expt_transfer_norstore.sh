#!/bin/bash
set -e
NSUSER=$USER
NSBASEDIR=/scratch/$NSUSER/
#NSBASEDIR=/projects/NS2993K/$NSUSER/
#NSBASEDIR=/norstore_osl/home/$NSUSER/

msg="
Usage: $(basename $0) relative_path_to_storage
 
Example: $(basename $0) 
 
Purpose: Copies expt config and data files to norstore, placement is relative to $NSBASEDIR


Example:  
   $(basename $0) TP4a0.12 
puts experiment data in  $NSBASEDIR/TP4a0.12
"

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
   echo "$msg"
   tellerror "Need target dir as input"
   exit 1
fi
target=$1

# Set dryrun for testing
#dryrun="--dry-run"
dryrun=""

# 
exptdir=$(basename $(pwd))
cd $BASEDIR

# username assumed the same on local system and norstore
ssh $NSUSER@norstore  mkdir -p $NSBASEDIR/$target/
rsync -avh --exclude="SCRATCH" ./$exptdir/         $NSUSER@norstore:$NSBASEDIR/$target/$exptdir/ $dryrun
rsync -avh --exclude="SCRATCH" ./topo/             $NSUSER@norstore:$NSBASEDIR/$target/topo/     $dryrun
rsync -avh -R --exclude="SCRATCH" ./relax/$E/      $NSUSER@norstore:$NSBASEDIR/$target/relax/    $dryrun
rsync -avh -R --exclude="SCRATCH" ./force/**/$E/   $NSUSER@norstore:$NSBASEDIR/$target/force/    $dryrun

# Other dirs? subregion, nest?
