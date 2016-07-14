#!/bin/bash

# KAL - get X from input
if [ $# -ne 2 ] ; then
   echo "This script will set up files with the  climatology. The "
   echo "climatology will be given on z-levels, and interpolated in "
   echo "the horizontal - also onto land points. The resulting files"
   echo "are needed for setting up the relaxation files to be used in"
   echo "HYCOM."
   echo
   echo "The only input is experiment number and climatology"
   echo "climatology can be phc, levitus or woa2013"
   echo
   echo "Example:"
   echo "$(basename $0) 01.0 phc"
   echo
   exit
fi
export CLIM_CHOICE=$2
export X=$1


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd $(dirname $0)/../ && pwd)/
source ${BASEDIR}/bin/common_functions.sh || { echo "Could not source ${BASEDIR}/bin/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }
D=$BASEDIR/relax/${CLIM_CHOICE}/
S=$D/SCRATCH
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

# Check that pointer to HYCOM_ALL is set
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
# --- Interpolate Levitus climtology to a HYCOM region.
# --- KSIGMA=0 for _sig0_ fields
#
#
# --- P is primary path,
# --- S is scratch   directory,
# --- D is permanent directory,
# --- PHC_PATH is PHC directory

if [ ${CLIM_CHOICE} == "phc" ] ; then
   CLIM_PATH=${PHC_PATH}
   CLIM_TITLE="PHC Climatology (monthly)"
elif [ ${CLIM_CHOICE} == "levitus" ] ; then
   CLIM_PATH=${LEVITUS_PATH}
   CLIM_TITLE="Levitus Climatology (monthly)"
elif [ ${CLIM_CHOICE} == "woa2013" ] ; then
   CLIM_PATH=${WOA2013_PATH}
   CLIM_TITLE="WOA2013 Climatology (monthly)"
else 
   echo "Unknown climatology ${CLIM_CHOICE}"
   exit 2
fi




[  "${CLIM_PATH}" == ""  ] && { echo "CLIM_PATH is not set " ; exit 1; }
[ ! -d ${CLIM_PATH} ] && { echo "Cant find CLIM_PATH directory at ${CLIM_PATH}" ; exit 1; }

copy_grid_files $S
cp $P/blkdat.input $S || tellerror "No $P/blkdat.input file" 
export KSIGMA=$(blkdat_get blkdat.input thflag)

touch   z_woa_dynamic && /bin/rm z_woa_dynamic
#
# --- 12 months
#
for MM  in  01 02 03 04 05 06 07 08 09 10 11 12 ; do

   if [ -f ${D}/temp_sig${KSIGMA}_m${MM}.b -a -f ${D}/temp_sig${KSIGMA}_m${MM}.a -a -f ${D}/saln_sig${KSIGMA}_m${MM}.b -a -f ${D}/saln_sig${KSIGMA}_m${MM}.a ]
   then
      echo "Found existing z-climatology files  ${D}/[temp|saln]_sig${KSIGMA}_m${MM}.b  - delete these first if you want to regenerate them"
      continue
   fi

   #
   # --- Input.
   #
   # --- Download Levitus from:
   # --- ftp://obelix.rsmas.miami.edu/awall/hycom/levitus_for_hycom.tar.gz
   #
   touch      fort.71 fort.73
   /bin/rm -f fort.71 fort.73
   ${pget} ${CLIM_PATH}/r_m${MM}.d fort.71 || { echo "Cant copy ${CLIM_PATH}/r_m${MM}.d "; exit 1 ; }
   ${pget} ${CLIM_PATH}/s_m${MM}.d fort.73 || { echo "Cant copy ${CLIM_PATH}/r_m${MM}.d "; exit 1 ; }

   touch z_woa_dynamic
   if [ ! -s  z_woa_dynamic ] ; then
     ${pget} ${HYCOM_ALL}/relax/src/z_woa_dynamic .  || { echo "Cant copy z_woa_dynamic "; exit 1 ; }
   fi
   wait

   chmod a+rx z_woa_dynamic
   /bin/rm -f core
   touch core
   export FOR010A=fort.10A
   export FOR011A=fort.11A
   export FOR012A=fort.12A
   /bin/rm -f fort.10 fort.10A fort.11 fort.11A fort.12 fort.12A
   ./z_woa_dynamic <<E-o-D
    &AFTITL
     CTITLE = '1234567890123456789012345678901234567890',
     CTITLE = '${CLIM_TITLE}',
    /
    &AFFLAG
     ICTYPE =   2,
     KSIGMA =   ${KSIGMA},
     INTERP =   1,
     ITEST  =  40,
     JTEST  =  17,
     MONTH  = $MM,
    /
E-o-D
   #
   # --- Required Output, potential density and temperature.
   #
   ${pput} fort.10  ${D}/temp_sig${KSIGMA}_m${MM}.b
   ${pput} fort.10A ${D}/temp_sig${KSIGMA}_m${MM}.a
   ${pput} fort.12  ${D}/dens_sig${KSIGMA}_m${MM}.b
   ${pput} fort.12A ${D}/dens_sig${KSIGMA}_m${MM}.a
   #
   # --- Optional Output.
   #
   ${pput} fort.11  ${D}/saln_sig${KSIGMA}_m${MM}.b
   ${pput} fort.11A ${D}/saln_sig${KSIGMA}_m${MM}.a
   #
   # --- end of month foreach loop
   #
   /bin/rm fort.7[123]
done
#
# --- Delete all files.
#
#/bin/rm -f *
