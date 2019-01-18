#!/bin/bash
#
# Script will set up a new experiment, copying routines from the present experiment. 
#Input is the number of the new experiment. The new experiment has no data files, only 
#the experiment directory expt_X

if [ $# -ne 2 ] ; then
   echo "This script will create a new experiment subdirectory based on"
   echo "an already existing experiment directory. Some basic variables"
   echo "will be set in files EXPT.src and blkdat.input - the rest must be "
   echo "set by user. Script will stop if the new experiment directory "
   echo "exists."
   echo "Usage:"
   echo "   $(basename $0) old_experiment_number new_experiment_number"
   echo
   echo "Example:"
   echo "   $(basename $0) 01.0 01.1"
   echo
   exit 1
fi

# Must be in expt or region dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
elif [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
else
   echo "Could not find EXPT.src REGION.src . This script must be run in expt dir or region dir"
   exit 1
fi
exptold=$1
exptnew=$2

# Dont create expt if one already exists
if [ -a $BASEDIR/expt_$exptnew ] ; then
   echo "Directory or file $T$BASEDIR/expt_$exptnew already exists"
   exit 1
fi

# Dont create expt if one regional.grid files are missing
if [ ! -f $BASEDIR/topo/regional.grid.b ] ; then
   echo "Region does not seem to be set up.. "
   echo "Can not find the file $BASEDIR/topo/regional.grid.b "
   exit 1
fi

mkdir $BASEDIR/expt_$exptnew
rsync -av --exclude 'data' --exclude 'log' --exclude 'SCRATCH' $BASEDIR/expt_$exptold/* $BASEDIR/expt_$exptnew/

#cp -r $BASEDIR/expt_$exptold/subprogs $BASEDIR/expt_$exptnew/subprogs
mkdir $BASEDIR/expt_$exptnew/log
mkdir $BASEDIR/expt_$exptnew/data

# Set up new EXPT.src file
cd $BASEDIR/expt_$exptnew
mv EXPT.src EXPT.src.tmp
enew=$(echo $exptnew | sed "s/\.//")
cat EXPT.src.tmp | \
   sed "s/^X=.*/X=\"$exptnew\"                # X based on dir name (expt_$exptnew)/" | \
   sed "s/^E=.*/E=\"$enew\"                 # E is X without \".\"/  " \
   > EXPT.src
rm EXPT.src.tmp

source ../REGION.src
source EXPT.src

# Get grid size from regional file
export IDM=`grep "'idm   '" $BASEDIR/topo/regional.grid.b | awk '{printf("%1d", $1)}'`
export JDM=`grep "'jdm   '" $BASEDIR/topo/regional.grid.b | awk '{printf("%1d", $1)}'`


# Set up new blkdat.input
mv blkdat.input blkdat.input.tmp
cat blkdat.input.tmp | \
   sed "s/^[ ]*[0-9]*[\t ]*'iexpt '/ $enew\t  'iexpt '/" |\
   sed "s/^[ ]*[0-9]*[\t ]*'idm   '/ $IDM\t  'idm   '/"  |\
   sed "s/^[ ]*[0-9]*[\t ]*'jdm   '/ $JDM\t  'jdm   '/" \
   > blkdat.input
rm blkdat.input.tmp

# Modify pbs job scripts with descriptive id
for i in pbsjob*.sh ; do
   mv $i $i.tmp
   cat $i.tmp | sed "s/PBS[ ]*-N.*$/PBS -N ${R}_X${X}/"  > $i
done

echo
echo
echo
echo "setup  done, now check 'iexpt' number (and any other parameter)  in:"
echo "   $BASEDIR/expt_$exptnew/blkdat.input "
echo "Also check $BASEDIR/EXPT.src"
echo "Then its time to set topo, relax, forcing etc etc"
echo "Need to edit the flag in $BASEDIR/Build_V${V}_X$exptnew and compile as well"
