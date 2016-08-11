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
if [ -a $TARGETDIR/$regname ] ; then
   echo "Directory or file $TARGETDIR/$regname already exists"
   exit 1
fi

BASEDIR=$(dirname $0)/..
cd $BASEDIR
BASEDIR=$(pwd)/
cd - > /dev/null


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

# start copying bin 
bindirs=$(find $BASEDIR -type d -name bin )
for i in $bindirs ; do
   relpath=$(echo $i | sed "s;${BASEDIR};;")

   # Check relative path
   cd $BASEDIR
   if [ ! -d $relpath ] ; then
      echo "Error in relative path generation - you will have to copy by hand"
      exit 1
   fi
   cd - > /dev/null

   #echo "test0: $i $relpath"
   mkdir -p $TARGETDIR/$regname/$relpath
   cp -r  $i/* $TARGETDIR/$regname/$relpath
done

# start copying READMEs
bindirs=$(find $BASEDIR -type f -name "README*" )
for i in $bindirs ; do
   relpath=$(echo $i | sed "s;${BASEDIR};;")

   # Check relative path
   cd $BASEDIR
   if [ ! -f $relpath ] ; then
      echo "Error in relative path generation - you will have to copy by hand"
      exit 1
   fi
   cd - > /dev/null

   #echo "test1: $i $relpath"
   mkdir -p $(dirname $TARGETDIR/$regname/$relpath)
   cp -r  $i $TARGETDIR/$regname/$relpath
done

# Copy files in top-level directory (only files)
cp $BASEDIR/rmu.in $TARGETDIR/$regname/
cp $BASEDIR/REGION.src $TARGETDIR/$regname/
#cp $BASEDIR/README.KAL $TARGETDIR/$regname/


# Copy all experiment dirs, nut not their data subdirectory
for i in $BASEDIR/expt_* ; do
   #echo "test2: $i"
   newdir=$TARGETDIR/$regname/$(basename $i)/
   [ ! -d $newdir ] && mkdir $newdir
   [ ! -d $newdir/log ] && mkdir $newdir/log
   [ ! -d $newdir/data ] && mkdir $newdir/data
   cp -r $i/subprogs $TARGETDIR/$regname/$(basename $i)/
   cp $i/* $TARGETDIR/$regname/$(basename $i)/
   #cp    $i/blkdat.input $TARGETDIR/$regname/$(basename $i)/
   #cp    $i/EXPT.src $TARGETDIR/$regname/$(basename $i)/
   #cp    $i/infile* $TARGETDIR/$regname/$(basename $i)/
   #cp    $i/*.sh $TARGETDIR/$regname/$(basename $i)/
   #cp    $i/ports.input $TARGETDIR/$regname/$(basename $i)/
   #cp    $i/README $TARGETDIR/$regname/$(basename $i)/
   #cp    $i/rivers.dat $TARGETDIR/$regname/$(basename $i)/
done

# Set up EXPT.src so that it is correct
cat  $TARGETDIR/$regname/REGION.src | sed "s/R=.*/R=$regname/" >  $TARGETDIR/$regname/REGION.src.new
mv $TARGETDIR/$regname/REGION.src.new $TARGETDIR/$regname/REGION.src

