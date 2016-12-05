#!/bin/bash
#
# Script will set up a new region, copying routines from the present region. 
#Input is name of new region. The new region will be emptied of data files

if [ $# -ne 2 ] ; then
   echo "This script will set up the directory structure for a new region"
   echo "The directory structure is retrieved by following the path of this script"
   echo "Input to this script is a descriptive name of the region, and "
   echo "optionally a path to where the new region directory structure will be "
   echo "placed. The new region will not have any data files, but will have the scripts"
   echo "needed to import/create them."
   echo
   echo "Example:"
   echo "   $(basename $0) TP3a0.12 /work/knutali/hycom/"
   echo "Will create a new region directory structure under /work/knutali/hycom/"
   echo
   exit 1
else
   regname=$1
   if [ $# -eq 2 ] ; then
      TARGETDIR=$2
   else
      TARGETDIR=../
   fi
fi

# Dont create region if one already exists
# Must be in expt or region dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
elif [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
else
   echo "Could not find EXPT.src or REGION.src. This script must be run in expt  opr region dir"
   exit 1
fi

TARGETDIR=$(cd $TARGETDIR && pwd)
if [ -a $TARGETDIR/$regname ] ; then
   echo "Directory or file $TARGETDIR/$regname already exists"
   exit 1
fi

# Safety check - see that basedir does not contain ";"
if echo $BASEDIR | grep ";" ; then
   echo "$BASEDIR contains the letter \";\" - this script will not work"
   exit 1
fi

# Safety check - see that regname does not contain "_"
if echo $regname | grep "_" ; then
   echo "$regname contains the letter \"_\" - this script will not work"
   exit 1
fi

mkdir -p $TARGETDIR/$regname

# Copy files in top-level directory (only files)
cp $BASEDIR/* $TARGETDIR/$regname/

# Copy all experiment dirs, but not their data subdirectory
for i in $BASEDIR/expt_* ; do
   #echo "test2: $i"
   newdir=$TARGETDIR/$regname/$(basename $i)/
   [ ! -d $newdir ] && mkdir $newdir
   [ ! -d $newdir/log ] && mkdir $newdir/log
   [ ! -d $newdir/data ] && mkdir $newdir/data
   #cp -r $i/subprogs $TARGETDIR/$regname/$(basename $i)/
   cp $i/* $TARGETDIR/$regname/$(basename $i)/

   # Set up pbsjob.sh with suitable job name
   source $i/EXPT.src
   cat $i/pbsjob.sh | sed -s "s/^#PBS[ ]*-N[ ]*.*/#PBS -N ${regname}_${X}/" >  $TARGETDIR/$regname/$(basename $i)/pbsjob.sh

done

# Set up EXPT.src so that it is correct
cat  $TARGETDIR/$regname/REGION.src | sed "s/^[ ]*export[ ]\+R=.*/export R=$regname/" >  $TARGETDIR/$regname/REGION.src.new
mv $TARGETDIR/$regname/REGION.src.new $TARGETDIR/$regname/REGION.src

echo "Base structure set up, note that topo subdir is empty, so you may want to start with that ..."

