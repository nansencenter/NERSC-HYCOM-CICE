#!/bin/bash

# KAL - get X from input
if [ $# -ne 3 ] ; then
   echo "This script will set up relaxation files from the NORCPM monthly files. The "
   echo "fileds will be given on z-levels, and interpolated in "
   echo "the horizontal - also onto land points. The resulting files"
   echo "are needed for setting up the relaxation files to be used in"
   echo "HYCOM."
   echo
   echo "The input is ksigma (0 or 2) and the start and end year"
   echo "netcdf-filed have been regridded to a regular grid from NORCPM"
   exit
fi
export KSIGMA=$1
export syear=$2
export eyear=$3

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly

# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
   echo $BASEDIR
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi

#cd $(dirname $0)/../../ 
#export BASEDIR=$(cd $(dirname $0)/.. && pwd)/
#cd -
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }

# Check that pointer to HYCOM_ALL is set (from REGION.src)
if [ -z ${HYCOM_ALL} ] ; then
   echo "Environment not set "
   exit
else
   if [ ! -d ${HYCOM_ALL} ] ; then
      echo "HYCOM_ALL not properly set up"
      echo "HYCOM_ALL not a directory at ${HYCOM_ALL}"
      exit
   fi
fi




#set -x
pget=cp
pput=cp

#
# --- Interpolate NORCPM filed to a HYCOM region.
# --- KSIGMA=0 for _sig0_ fields
#
#
# --- P is primary path,
# --- S is scratch   directory,
# --- D is permanent directory,
# --- L is Levitus   directory (LEVITUS_PATH)
#
D=$BASEDIR/relax/norcpm/
S=$D/SCRATCH
L=$NORCPM_CLIM_PATH/
[  "${L}" == ""  ] && { echo "NORCPM_CLIM_PATH is not set " ; exit 1; }
[ ! -d $L ] && { echo "Cant find NORCPM_CLIM_PATH directory $L" ; exit 1; }

mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


touch   z_norcpm
/bin/rm z_norcpm
touch regional.grid.a regional.grid.b
rm    regional.grid.a regional.grid.b
#
# --- 12 months
#
for ((year=syear; year<=eyear; year+=1)); do
for ((MM=1; MM<=12; MM+=1)); do
CMM=`echo -n 0${MM} | tail -2c` 
echo MONTH=$CMM
#for MM  in  01 02 03 04 05 06 07 08 09 10 11 12 ; do
   #
   # --- Input.
   #
   # --- Download Levitus from:
   # --- ftp://obelix.rsmas.miami.edu/awall/hycom/levitus_for_hycom.tar.gz
   #
   
   export CDF_TEMP=${NORCPM_CLIM_PATH}thetao_Omon_NorCPM1_historical_r20i1p1f1_gr_${year}s${CMM}_rgrid.nc
   export CDF_SALN=${NORCPM_CLIM_PATH}so_Omon_NorCPM1_historical_r20i1p1f1_gr_${year}s${CMM}_rgrid.nc

   if [ ! -s regional.grid.b ] ; then
     ${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Cant copy regional.grid.b "; exit 1 ; }
   fi
   if [ ! -s regional.grid.a ] ; then
     ${pget} ${BASEDIR}//topo/regional.grid.a regional.grid.a || { echo "Cant copy regional.grid.a "; exit 1 ; }
   fi
   touch z_norcpm
   if [ ! -s  z_woa2013_bio ] ; then
     ${pget} ${MSCPROGS}/src/Relax/z_norcpm .  || { echo "Cant copy z_woa2013_bio "; exit 1 ; }
   fi
   wait

   chmod a+rx z_norcpm
   /bin/rm -f core
   touch core

   /bin/rm -f fort.10 fort.010a fort.11 fort.011a fort.12 fort.012a

   ./z_norcpm <<E-o-D
    &AFTITL
     CTITLE = '1234567890123456789012345678901234567890',
     CTITLE = 'NORCPM     monthly',
    /
    &AFFLAG
     ICTYPE =   3,
     SIGVER =   1,
     MONTH  =   ${MM}
     INTERP =   0,
     ITEST  =  40,
     JTEST  =  17,
    /
E-o-D
   #
   # --- Required Output, potential density and temperature.
   #
   ${pput} fort.10  ${D}/temp_sig${KSIGMA}_y${year}_m${CMM}.b
   ${pput} fort.010a ${D}/temp_sig${KSIGMA}_y${year}_m${CMM}.a
   ${pput} fort.11  ${D}/saln_sig${KSIGMA}_y${year}_m${CMM}.b
   ${pput} fort.011a ${D}/saln_sig${KSIGMA}_y${year}_m${CMM}.a
   ${pput} fort.12  ${D}/dens_sig${KSIGMA}_y${year}_m${CMM}.b
   ${pput} fort.012a ${D}/dens_sig${KSIGMA}_y${year}_m${CMM}.a
   #
   # --- Optional Output.
   #
   #
   # --- end of month foreach loop
   #
done
done
#
# --- Delete all files.
#
#/bin/rm -f *
