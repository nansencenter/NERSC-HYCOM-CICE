#!/bin/bash
EDIR=$(cd $(dirname $0) && pwd)/         # Location of this script
BASEDIR=$(cd $(dirname $0)/.. && pwd)/   # Location of basedir
source $BASEDIR/bin/common_functions.sh
source $BASEDIR/REGION.src
source $EDIR/EXPT.src
echo "Logs can be found in $EDIR/log"
echo ".."

#Various settings
myclim="woa2013" # Climatology to use

# Get some useful info
THFLAG=$(blkdat_get blkdat.input thflag)
IDM=$(blkdat_get blkdat.input idm)
JDM=$(blkdat_get blkdat.input jdm)

#Get statement function include file from SIGVER
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
echo "climatology = $myclim"

# Set up rel path and stmt fnc
stmt=stmt_fns_SIGMA${MYTHFLAG}_${TERMS}term.h
targetdir=$EDIR/build/src_${V}ZA-${TERMS2}Tsig${MYTHFLAG}-i-sm-sse_relo_mpi/
targetconfdir=$EDIR/build/config/

# Sanity check
if [ $THFLAG -ne $MYTHFLAG ] ; then
   echo "thflag in blkdat.input ($THFLAG) does not match thflag from eq of state ($MYTHFLAG)"
   exit 1
fi


# Create hycom_all executables
cd $EDIR
dir="$BASEDIR/../hycom_ALL/hycom_2.2.72_ALL/"
cd $dir
echo "compiling hycom_all in $dir"
csh Make_all.com ARCH=amd64 > $EDIR/log/ref_hycom_all.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

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
echo "compiling hycom_cice in $dir"
env RES=gx3 GRID=${IDM}x${JDM} csh Make_cice.csh > $EDIR/log/ref_hycom.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."



# Create z-level relaxation files
cd $EDIR
echo "z climatology"
../bin/z_generic.sh ${X} $myclim > $EDIR/log/ref_z_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create relaxation files on hybrid coordinates
cd $EDIR
echo "hybrid climatology"
../bin/relaxi.sh ${X} $myclim 2    > $EDIR/log/ref_hybrid_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create relaxation mask
echo "relaxation mask"
cat <<EOF | ../bin/relax_rmu.sh ${X}  > $EDIR/log/ref_rmu_mask.out 2>&1
F
F
F
T
F
20
20
EOF
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create simple river forcing
cd $EDIR
echo "river forcing"
../bin/river_nersc.sh ${X} 100 300 > $EDIR/log/ref_river_nersc.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create kpar file
echo "kpar forcing"
../bin/seawifs_mon_kpar.sh > $EDIR/log/ref_seawifs.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create tiling 
echo "grid tiling"
../bin/tile_grid.sh -2 -2 9.5 01 > $EDIR/log/ref_tiling.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

echo "If things went fine, you can now generate test forcing like this: "
echo "    ../bin/atmo_synoptic_new.sh ${X} erai 2015-01-02T00:00:00  2015-01-05T00:00:00"
echo "Then edit the job script pbsjob.sh to read and make sure expt_preprocess.sh is called like this"
echo "    ../bin/expt_preprocess_new.sh 2015-01-02T00:00:00 2015-01-05T00:00:00 --init "
echo
echo "TODO: Current limitations on synoptic forcing routines mean that only erai forcing from 2015 can be used."

