#!/bin/bash
EDIR=$(cd $(dirname $0) && pwd)/
BASEDIR=$(cd $(dirname $0)/.. && pwd)/
source $BASEDIR/bin/common_functions.sh
source $BASEDIR/REGION.src
source $EDIR/EXPT.src
THFLAG=$(blkdat_get blkdat.input thflag)
echo "Logs can be found in $EDIR/log"
echo ".."

# Create hycom_all exacutables
cd $EDIR
dir="$BASEDIR/../hycom_ALL/hycom_2.2.72_ALL/"
cd $dir
echo "compiling hycom_all in $dir"
csh Make_all.com ARCH=amd64 > $EDIR/log/ref_hycom_all.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create hycom exacutable
cd $EDIR
dir=$BASEDIR/../hycom/RELO/src_${V}ZA-07Tsig${THFLAG}-i-sm-sse_relo_mpi/
cd $dir
source set_env.sh
echo "compiling hycom_cice in $dir"
csh Make_cice.csh > $EDIR/log/ref_hycom.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create z-level relaxation files
cd $EDIR
echo "z climatology"
../bin/z_generic.sh 01.0 woa2013 > $EDIR/log/ref_z_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create relaxation files on hybrid coordinates
cd $EDIR
echo "hybrid climatology"
../bin/relaxi.sh 01.0 woa2013 2    > $EDIR/log/ref_hybrid_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create relaxation mask
echo "relaxation mask"
cat <<EOF | ../bin/relax_rmu.sh 01.0  > $EDIR/log/ref_rmu_mask.out 2>&1
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
../bin/river_nersc.sh 01.0 100 300 > $EDIR/log/ref_river_nersc.out 2>&1
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
echo "    ../bin/atmo_synoptic_new.sh 01.0 erai 2015-01-02T00:00:00  2015-01-05T00:00:00"
echo "Then edit the job script pbsjob.sh to read and make sure expt_preprocess.sh is called like this"
echo "    ../bin/expt_preprocess_new.sh 2015-01-02T00:00:00 2015-01-05T00:00:00 --init "
echo
echo "TODO: Current limitations on synoptic forcing routines mean that only erai forcing from 2015 can be used."

