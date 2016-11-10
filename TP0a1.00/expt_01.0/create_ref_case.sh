#!/bin/bash



myclim="woa2013" # Climatology to use

# What ARCH value to pass to hycom_all compilation script
# NB: On some machines (cray/unicos) you may have to modify module settings as well


#hycom_all_arch=gfortran   # Most linux machines have gnu fortran
hycom_all_arch=amd64       # Uses pgi compiler
hycom_arch=amd64





# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
BINDIR="../../bin"
EDIR=$(pwd)/                           # Location of this script
BASEDIR=$(cd $(dirname $0)/.. && pwd)/ # Location of basedir
source $BASEDIR/REGION.src
source $EDIR/EXPT.src
source ../../bin//common_functions.sh
echo "Logs can be found in $EDIR/log"
[ ! -d $EDIR/log ] && mkdir $EDIR/log/
echo ".."

#Various settings
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
THFLAG=$(blkdat_get blkdat.input thflag)
IDM=$(blkdat_get blkdat.input idm)
JDM=$(blkdat_get blkdat.input jdm)

# Create hycom_all executables
cd $EDIR
dir="$HYCOM_ALL"
cd $dir
echo "compiling hycom_all in $dir"
csh Make_all.com ARCH=amd64 > $EDIR/log/ref_hycom_all.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

cd $EDIR
echo "Compiling hycom_cice"
$BINDIR/compile_model.sh > $EDIR/log/ref_hycom.out 2>&1
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create z-level relaxation files
cd $EDIR
echo "z climatology"
$BINDIR/z_generic.sh $myclim > $EDIR/log/ref_z_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create relaxation files on hybrid coordinates
cd $EDIR
echo "hybrid climatology"
$BINDIR/relaxi.sh $myclim   > $EDIR/log/ref_hybrid_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create relaxation mask
echo "relaxation mask"
cat <<EOF | $BINDIR/relax_rmu.sh ${X}  > $EDIR/log/ref_rmu_mask.out 2>&1
F
F
F
T
T
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
$BINDIR/river_nersc.sh 100 300 ../../input/rivers.dat > $EDIR/log/ref_river_nersc.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create kpar file
echo "kpar forcing"
$BINDIR/seawifs_mon_kpar.sh > $EDIR/log/ref_seawifs.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create tiling. 
echo "grid tiling"
$BINDIR/tile_grid.sh -2 -2 9.5 ${T} > $EDIR/log/ref_tiling.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

echo "If things went fine, you can now generate test forcing like this: "
echo "    $BINDIR/atmo_synoptic.sh erai 2015-01-02T00:00:00  2015-01-05T00:00:00"
echo "Then edit the job script pbsjob.sh to read and make sure expt_preprocess.sh is called like this"
echo "    $BINDIR/expt_preprocess.sh 2015-01-02T00:00:00 2015-01-05T00:00:00 --init "
echo
echo "TODO: Current limitations on synoptic forcing routines mean that only erai forcing from 2015 can be used."

