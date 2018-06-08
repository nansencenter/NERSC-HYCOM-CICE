#!/bin/bash

myclim="woa2013" # Climatology to use

# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
EDIR=$(pwd)/                           # Location of this script
BASEDIR=$(cd $(dirname $0)/.. && pwd)/ # Location of basedir
source $BASEDIR/REGION.src
source $EDIR/EXPT.src

# Create z-level relaxation files
cd $EDIR
echo "z climatology"
z_generic.sh $myclim > $EDIR/log/ref_z_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_z_relax.out"
echo ".."

# Create relaxation files on hybrid coordinates
cd $EDIR
echo "hybrid climatology"
relaxi.sh $myclim   > $EDIR/log/ref_hybrid_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_hybrid_relax.out"
echo ".."

# Create relaxation mask
echo "relaxation mask"
cat <<EOF | relax_rmu.sh ${X}  > $EDIR/log/ref_rmu_mask.out 2>&1
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
#river_nersc.sh 100 300 $INPUTDIR/rivers.dat > $EDIR/log/ref_river_nersc.out 2>&1
river_trip_bio.sh erai > $EDIR/log/ref_river_nersc.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_river_nersc.out"
echo ".."

# Create kpar file
echo "kpar forcing"
seawifs_mon_kpar.sh > $EDIR/log/ref_seawifs.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_seawifs.out"
echo ".."

# Create tiling. 
echo "grid tiling"
tile_grid.sh -2 -2 ${T} > $EDIR/log/ref_tiling.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_tiling.out" 
echo ".."

echo "If things went fine, you can now generate test forcing like this: "
echo "    atmo_synoptic.sh erai 2001-01-01T00:00:00  2001-01-05T00:00:00"
echo "Then edit the job script pbsjob.sh to read and make sure expt_preprocess.sh is called like this"
echo "    expt_preprocess.sh 2001-01-01T00:00:00 2001-01-05T00:00:00 --init "
echo
echo "TODO: Current limitations on synoptic forcing routines mean that only erai forcing from 2015 can be used."

