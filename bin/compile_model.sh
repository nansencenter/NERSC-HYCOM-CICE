#!/bin/bash
#

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


# Get some useful info
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
else
   echo "SIGVER = $SIGVER"
   echo "So far only 7 term eq of state is supported (SIGVER=1 or 2) is supported"
   exit 1
fi
TERMS2=$(echo 0$TERMS | tail -c3)
echo "SIGVER      = $SIGVER .There are $TERMS terms in equation of state"

# Set up rel path and stmt fnc
stmt=stmt_fns_SIGMA${MYTHFLAG}_${TERMS}term.h
targetdir=$EDIR/build/src_${V}ZA-${TERMS2}Tsig${MYTHFLAG}-i-sm-sse_relo_mpi/
targetconfdir=$EDIR/build/config/

# Sanity check
if [ $THFLAG -ne $MYTHFLAG ] ; then
   echo "thflag in blkdat.input ($THFLAG) does not match thflag from eq of state ($MYTHFLAG)"
   exit 1
fi

# Create hycom executable. Copy code to expt dir
sourcedir=$BASEDIR/../hycom/RELO/src_${V}ZA-07Tsig0-i-sm-sse_relo_mpi/ # Yes, we always use this version
sourceconfdir=$BASEDIR/../hycom/RELO/config/
[ ! -d $EDIR/build/ ] && mkdir $EDIR/build
cp -r -L $sourcedir $targetdir
cp -r -L $sourceconfdir $targetconfdir

# Create hycom executable. Set up correct eq of state
cd $targetdir
source set_env.sh
rm stmt_fns.h
ln -s ALT_CODE/$stmt stmt_fns.h

# Create hycom executable. Set up correct domain size in CICE code, and pass this through to the compile script
# NB: hycom is domain-independent, CICE is not. So the grid size info needs to be there. It is passed on to the CICE
# compilation script by hycoms makefile
echo "compiling hycom_cice in $targetdir"
env RES=gx3 GRID=${IDM}x${JDM} csh Make_cice.csh > $EDIR/log/ref_hycom.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."
