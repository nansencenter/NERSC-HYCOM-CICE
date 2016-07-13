#!/bin/bash
#
#set -x
pget=cp
pput=cp

# Experiment number needed
if [ $# -ne 2 ] ; then
   echo "This script will set up climatology forcing files for hycom. The"
   echo "climatology will be based on either NCEP, ERA40 or old raw climatology"
   echo "files, and are interpolated onto the model grid. Input to this routine"
   echo "is the experiment number (exampel 01.0), and the climatology option"
   echo "(ncepr, era40 or old)."
   echo
   echo "Example:"
   echo "   $(basename $0) 01.0 era40"
   echo "Will create climatology files for experiment 01.0 basexd on ERA40 data"
   exit
fi
export X=$1
export clim=$2

if [ "$clim" != "era40" -a "$clim" != "ncepr" -a "$clim" != "old"  ] ; then
   echo "Unknown climatology $clim"
   echo " Climatology can be era40, ncepr and old"
   exit 1
fi

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd $(dirname $0)/.. && pwd)/  
echo $BASEDIR
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }
D=$BASEDIR/force/nersc_${clim}/
S=$D/SCRATCH
mkdir -p $S
mkdir -p $S/Data


# Check that pointer to MSCPROGS is set
if [ -z ${MSCPROGS} ] ; then
   echo "MSCPROGS Environment not set "
   exit
else
   if [ ! -d ${MSCPROGS} ] ; then
      echo "MSCPROGS not properly set up"
      echo "MSCPROGS not a directory at ${MSCPROGS}"
      exit
   fi
fi

# Check that environment variables are set
if [ "$clim" == "era40" ] ; then
   [ -z ${ERA40_CLIM_PATH} ] && { echo "ERA40_CLIM_PATH is not set " ; exit 1 ; }
   [ ! -d ${ERA40_CLIM_PATH} ] && { echo "Directory ${ERA40_CLIM_PATH} does not exist " ; exit 1 ; }
elif [ "$clim" == "ncepr" ] ; then
   [ -z ${NCEP_CLIM_PATH} ] && { echo "NCEP_CLIM_PATH is not set " ; exit 1 ; }
   [ ! -d ${NCEP_CLIM_PATH} ] && { echo "Directory ${NCEP_CLIM_PATH} does not exist " ; exit 1 ; }
elif [ "$clim" != "old" ] ; then
   # Data for "old" is under MSCPROGS directory
   echo " Unknown climatology $clim"
   exit 1
fi


#
#
# --- Create HYCOM atmospheric forcing (era40 climatology)
#
#
# --- S is scratch directory,
# --- D is permanent directory,
#
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Input.
#
if [ "$clim" == "old" ] ; then
   touch Data/sss_nodc.dat Data/sst_nodc.dat
   rm    Data/sss_nodc.dat Data/sst_nodc.dat
   cp ${BASEDIR}/relax/${E}/sss_nodc.dat Data/  || \
     { echo "Could not get file sss_nodc.dat (need to run old_lev_nersc.sh in relax first)" ; exit 1 ; }
   cp ${BASEDIR}/relax/${E}/sst_nodc.dat Data/  || \
      { echo "Could not get file sst_nodc.dat (need to rub old_lev_nersc.sh in relax first)" ; exit 1 ; }
fi

touch regional.grid.a regional.grid.b regional.depth.a regional.depth.b 
rm    regional.grid.a regional.grid.b regional.depth.a regional.depth.b 
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b     || { echo "Could not get regional.grid.b file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a     || { echo "Could not get regional.grid.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Could not get regional.depth.a file " ; exit 1 ; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }
#${pget} ${BASEDIR}/expt_${X}/blkdat.input blkdat.input      || { echo "Could not get file blkdat.input " ; exit 1 ; }
#${pget} ${BASEDIR}/expt_${X}/infile.in infile.in            || { echo "Could not get file infile.in " ; exit 1 ; }


if [ "$clim" == "old" ] ; then
   echo "Copying old Forcing to $S"
   cp -r $MSCPROGS/data/Data $S
   cp -r $MSCPROGS/data/Data/fields.inp* $S
   $MSCPROGS/bin_setup/old_newf;
fi


$MSCPROGS/bin_setup/forfun_nersc-2.2.37 month $clim ||  { echo "Error running forfun_nersc-2.2.37 " ; exit 1 ; }

# The nersc era40 forcing is region-independent 
mkdir -p $D/$E  || { echo " Could not create data  destination dir $D/$E" ; exit 1;}
for i in forcing.*.[ab] ; do
   new=$(echo $i | sed "s/^forcing\.//")
   mv $i $D/$E/$new
done

echo
echo "Done. $clim climatology files in $D/$E"

