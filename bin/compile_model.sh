#!/bin/bash
#
# Script for quickly setting up model source code and compiling it. 
update=$1


# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
export EDIR=$(pwd)
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }


# Get some useful info from blkdat.input
THFLAG=$(blkdat_get blkdat.input thflag)
IDM=$(blkdat_get blkdat.input idm)
JDM=$(blkdat_get blkdat.input jdm)

#Get statement function include file from SIGVER (in EXPT.src)
if [ $SIGVER -eq 1 ] ; then
   TERMS=7
   MYTHFLAG=0
elif [ $SIGVER -eq 2 ] ; then
   TERMS=7
   MYTHFLAG=2
elif [ $SIGVER -eq 3 ] ; then
   TERMS=9
   MYTHFLAG=0
elif [ $SIGVER -eq 4 ] ; then
   TERMS=9
   MYTHFLAG=2
elif [ $SIGVER -eq 5 ] ; then
   TERMS=17
   MYTHFLAG=0
elif [ $SIGVER -eq 6 ] ; then
   TERMS=17
   MYTHFLAG=2
else
   echo "SIGVER = $SIGVER"
   echo "So far only 7 term eq of state is supported (SIGVER=1 or 2) is supported"
   exit 1
fi
TERMS2=$(echo 0$TERMS | tail -c3)
echo "SIGVER      = $SIGVER .There are $TERMS terms in equation of state"

# Sanity check
if [ $THFLAG -ne $MYTHFLAG ] ; then
   echo "thflag in blkdat.input ($THFLAG) does not match thflag from eq of state ($MYTHFLAG)"
   exit 1
fi

# Set up rel path and statment function name fnc
stmt=stmt_fns_SIGMA${MYTHFLAG}_${TERMS}term.h
targetdir=$(source_dir $V $TERMS $THFLAG)
targetdir=$EDIR/build/$targetdir
targetconfdir=$EDIR/build/config/

# Create hycom executable. Copy code to expt dir
sourcedir=$BASEDIR/../hycom/RELO/src_${V}ZA-07Tsig0-i-sm-sse_relo_mpi/ # Yes, we always use this version
sourceconfdir=$BASEDIR/../hycom/RELO/config/
#if [ ! -d $EDIR/build/ ] ; then 
#   mkdir $EDIR/build
#   cp -r -L $sourcedir $targetdir
#   cp -r -L $sourceconfdir $targetconfdir
#   echo "build dir $EDIR/build not found. Setting it up with repo code from $sourcedir"
#else 
#   echo "build dir $EDIR/build found. Using code in that subdirectory"
#fi
if [ ! -d $EDIR/build/ ] ; then 
   mkdir $EDIR/build
   echo "build dir $EDIR/build not found. Setting it up with repo code from $sourcedir"
   rsync -avhL $sourcedir/ $targetdir/
   rsync -avhL $sourceconfdir/ $targetconfdir/
else 
   if [ "$update" == "update" ] ; then
      echo "build dir $EDIR/build found. Updating code in that subdirectory"
      rsync -avhL $sourcedir/ $targetdir/
      rsync -avhL $sourceconfdir/ $targetconfdir/
   else 
      echo "build dir $EDIR/build found. Using code in that subdirectory [not updated with code in $sourcedir]"
   fi
fi

# Create hycom executable. Set up correct eq of state
cd $targetdir
echo "Now setting up stmt_fns.h in $targetdir"
source set_env.sh
rm stmt_fns.h
ln -s ALT_CODE/$stmt stmt_fns.h

# Create hycom executable. Set up correct domain size in CICE code, and pass this through to the compile script
# NB: hycom is domain-independent, CICE is not. So the grid size info needs to be there. It is passed on to the CICE
# compilation script by hycoms makefile
echo "Now compiling hycom_cice in $targetdir." 
env RES=gx3 GRID=${IDM}x${JDM} csh Make_cice.csh 
exit $?
