#!/bin/bash
#######################################################################
# Start post-processing
set -e
trap 'abort' ERR

# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
if [ -z ${BASEDIR} ] ; then
   tellerror "BASEDIR Environment not set "
   exit
fi

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ ) - must be in
# experiment dir for this to work
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src "; exit 1; }
source $BINDIR/common_functions.sh  || { echo "Could not source common_functions.sh "; exit 1; }


[ ! -d $D ] && mkdir $D
[ ! -d $D/cice ] && mkdir $D/cice

# Print error if we cannot descend scratch dir $S
cd       $S || { cd $P ; echo "BADRUN" > log/hycom.stop ;  exit 1 ;}

touch   PIPE_DEBUG
/bin/rm PIPE_DEBUG

# Get names of archive files, restart files, etc
restarto=$(blkdat_get_string blkdat.input nmrsto "restart_out")
restarti=$(blkdat_get_string blkdat.input nmrsti "restart_in")
nmarcv=$(blkdat_get_string blkdat.input nmarcv "archv.")
nmarcs=$(blkdat_get_string blkdat.input nmarcs "archs.")
nmarcm=$(blkdat_get_string blkdat.input nmarcm "archm.")
#echo $restarto
#echo $restarti
#echo "test"
#exit 1


# Move restart and daily files to data location $D
#Delete link first to avoid recursive links
for i in $(ls -- ${restarto}*.[ab]) 
do
  if [ -L "${i}" ]
    then
       rm $i
    else 
       echo "Moving $i to $D"
       mv $i $D
  fi
done

# Move archive files
for i in $(ls ${nmarcv}*.[ab] ${nmarcs}*.[ab] ${nmarcm}*.[ab]) 
do
  if [ -L "${i}" ]
    then
       echo "$i is a sympolic link, I am removing it"
       rm $i
    else
       echo "Moving $i to  $D"
       mv $i $D
   fi
done

# Move cice archive files 
for i in $(ls -- cice/)  ; do
   if [ -L "cice/${i}" ]
     then
        echo "cice/$i is a sympolic link, I am removing it"
        rm cice/$i
     else
        echo "Moving cice/$i  to $D/cice/"
        mv cice/$i $D/cice/
   fi
done
[ -f summary.out ] && cp summary_out $D/


# Copy some files useful for analysis
for i in $(ls regional.grid.* regional.depth.* blkdat.input ice_in cice_*.nc) ; do
   cp $i $D
done

#
# --- HYCOM error stop is implied by the absence of a normal stop.
#
if  [ `tail -1 summary_out | grep -c "^normal stop"` == 0 ] ; then
  cd $P
  echo "BADRUN"  > log/hycom.stop
else
  cd $P
  echo "GOODRUN" > log/hycom.stop
fi
