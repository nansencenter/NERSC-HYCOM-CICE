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
BASEDIR=$(cd .. && pwd)/ # Location of basedir
source $BASEDIR/REGION.src
source $EDIR/EXPT.src



MERCID="MRCb0.25"
MERCROOT="../../"
MERCDIR="$MERCROOT/$MERCID"

NEMO_MESH_FILE="/work/shared/nersc/msc/mercator/GLORYS2V3_mesh_mask.nc"
NEMO_2D_FILES="/work/shared/nersc/msc/mercator/GLORYS2V4_1dAV_*grid2D*.nc"

# Create a region directory that will keep the Mercacator files. Set up in "HYCOM-fashion"
# to make it easier to use existing routines
region_new.sh $MERCID $MERCROOT
if [ $? -ne 0 ]  ; then
   echo " Could not create $MERCDIR, see above"
   #exit 1
fi

# Enter exp dir.
cd $MERCDIR/expt_01.0

# Get first temp file. Extract idm, jdm, kdm
files=(${NEMO_2D_FILES})
ffile=${files[0]}
ftfile=${ffile/grid2D/gridT}
idm=$(ncdump -h $ftfile | egrep "^[ 	]*x[ ]*="| sed "s/.*x[ ]*=[ ]*\([0-9]\+\).*/\1/")
jdm=$(ncdump -h $ftfile | egrep "^[ 	]*y[ ]*="| sed "s/.*y[ ]*=[ ]*\([0-9]\+\).*/\1/")
kdm=$(ncdump -h $ftfile | egrep "^[ 	]*deptht[ ]*="| sed "s/.*deptht[ ]*=[ ]*\([0-9]\+\).*/\1/") 
echo "idm=$idm jdm=$jdm kdm=$kdm (from $ftfile)"


# Edit blkdat.input and use new dimensions
cat blkdat.input  | sed  "s/^.*'idm   '.*/ ${idm}      'idm   ' = longitudinal array size/"  \
                  | sed  "s/^.*'jdm   '.*/ ${jdm}      'jdm   ' = latitudinal array size/"  \
                  | sed  "s/^.*'kdm   '.*/ ${kdm}      'kdm   ' = number of layers/" > blkdat.input.new

# Now run main script, that converts nemo netcdf files to hycom acrhv files
nemo_to_hycom.py ${NEMO_MESH_FILE} ${NEMO_2D_FILES}

exit 
# In the MERCID directory, run the remapping routine from mercator grid to this region/experiment. 
# this=the location of this script ...
isuba_gmapi.sh $BASEDIR





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
river_nersc.sh 100 300 ../../input/rivers.dat > $EDIR/log/ref_river_nersc.out 2>&1
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
echo "    expt_preprocess.sh 2015-01-02T00:00:00 2015-01-05T00:00:00 --init "

