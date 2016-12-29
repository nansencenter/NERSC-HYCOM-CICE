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
#NEMO_2D_FILES="/work/shared/nersc/msc/mercator/GLORYS2V4_1dAV_*grid2D*.nc"
NEMO_2D_FILES="/work/shared/nersc/msc/mercator/GLORYS2V4_1dAV_20130101_20130102_grid2D_R20130102.nc"

# Create a region directory that will keep the Mercacator files. Set up in "HYCOM-fashion"
# to make it easier to use existing routines
region_new.sh $MERCID $MERCROOT
if [ $? -ne 0 ]  ; then
   echo " Could not create $MERCDIR, see above"
   #exit 1
fi

# Enter exp dir.
MERCDIR=$(cd $MERCDIR && pwd)
cd $MERCDIR/expt_01.0

# Get first temp file. Extract idm, jdm, kdm
files=(${NEMO_2D_FILES})
ffile=${files[0]}
ftfile=${ffile/grid2D/gridT}
idm=$(ncdump -h $ftfile | egrep "^[ 	]*x[ ]*="| sed "s/.*x[ ]*=[ ]*\([0-9]\+\).*/\1/")
jdm=$(ncdump -h $ftfile | egrep "^[ 	]*y[ ]*="| sed "s/.*y[ ]*=[ ]*\([0-9]\+\).*/\1/")
kdm=$(ncdump -h $ftfile | egrep "^[ 	]*deptht[ ]*="| sed "s/.*deptht[ ]*=[ ]*\([0-9]\+\).*/\1/") 
echo "**idm=$idm jdm=$jdm kdm=$kdm (from $ftfile)"


# Edit blkdat.input and use new dimensions
echo "**Setting up blkdat.input for Mercator grid"
cat blkdat.input  | sed  "s/^.*'idm   '.*/ ${idm}      'idm   ' = longitudinal array size/"  \
                  | sed  "s/^.*'jdm   '.*/ ${jdm}      'jdm   ' = latitudinal array size/"  \
                  | sed  "s/^.*'kdm   '.*/ ${kdm}      'kdm   ' = number of layers/" > blkdat.input.new

# Now run main script, that converts nemo netcdf files to hycom acrhv files
mkdir data
cd data
echo "**Converting from NEMO netcdf to hycom archv"
nemo_to_hycom.py ${NEMO_MESH_FILE} ${NEMO_2D_FILES}
cd ..

# Copy grid and depth files to topo dir  in the MERCATOR dir..
mv data/regional.grid.* $MERCDIR/topo/
mv data/depth_bathy.a $MERCDIR/topo/depth_${MERCID}_01.a
mv data/depth_bathy.b $MERCDIR/topo/depth_${MERCID}_01.b

# In the MERCID directory, run the remapping routine from mercator grid to this region/experiment. 
# this=the location of this script ...
echo "**Calculating remap from NEMO to hycom grid"
isuba_gmapi.sh $BASEDIR

# Now we need to create a bathymetry that merges Mercator bathy and this bathy








