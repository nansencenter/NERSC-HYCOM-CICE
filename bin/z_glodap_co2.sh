#!/bin/bash

# KAL - get X from input
if [ $# -ne 1 ] ; then
   echo "This script will set up files with the WOA2013 climatology. The "
   echo "climatology will be given on z-levels, and interpolated in "
   echo "the horizontal - also onto land points. The resulting files"
   echo "are needed for setting up the relaxation files to be used in"
   echo "HYCOM."
   echo
   echo "The only input is ksigma (0 or 2)"
   echo ".d files to be read are created by bin/hycom_glodap_zfiles.py"
   exit
fi
export KSIGMA=$1

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
# --- Interpolate WOA2013 climtology to a HYCOM region.
# --- KSIGMA=0 for _sig0_ fields
#
#
# --- P is primary path,
# --- S is scratch   directory,
# --- D is permanent directory,
# --- L is Levitus   directory (LEVITUS_PATH)
#
D=$BASEDIR/relax/woa2013/
S=$D/SCRATCH
L=$GLODAPCO2_PATH/
[  "${L}" == ""  ] && { echo "GLODAPCO2_PATH is not set " ; exit 1; }
[ ! -d $L ] && { echo "Cant find GLODAPCO2_PATH directory $L" ; exit 1; }

mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


touch   z_glodap_co2
/bin/rm z_glodap_co2
touch regional.grid.a regional.grid.b
rm    regional.grid.a regional.grid.b
#
# --- 12 months
#
for MM  in  01 02 03 04 05 06 07 08 09 10 11 12 ; do
   #
   # --- Input.
   #
   # --- Download Levitus from:
   # --- ftp://obelix.rsmas.miami.edu/awall/hycom/levitus_for_hycom.tar.gz
   #
   touch      fort.71 fort.72
   /bin/rm -f fort.71 fort.72
   ${pget} $L/TCO2${MM}.d fort.71 || { echo "Cant copy $L/TCO2${MM}.d "; exit 1 ; }
   ${pget} $L/TAlk${MM}.d fort.72 || { echo "Cant copy $L/TAlk${MM}.d "; exit 1 ; }

   if [ ! -s regional.grid.b ] ; then
     ${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Cant copy regional.grid.b "; exit 1 ; }
   fi
   if [ ! -s regional.grid.a ] ; then
     ${pget} ${BASEDIR}//topo/regional.grid.a regional.grid.a || { echo "Cant copy regional.grid.a "; exit 1 ; }
   fi
   touch z_glodap_co2
   if [ ! -s  z_glodap_co2 ] ; then
       ${pget} ${MSCPROGS}/bin_setup/z_glodap_co2 .  || { echo "Cant copy z_glodap_co2 "; exit 1 ; }
   fi
   wait

   chmod a+rx z_glodap_co2
   /bin/rm -f core
   touch core
   export FOR010A=fort.10A
   export FOR011A=fort.11A
   export FOR012A=fort.12A
   export FOR013A=fort.13A
   /bin/rm -f fort.10 fort.10A fort.11 fort.11A fort.12 fort.12A fort.13 fort.13A
   ./z_glodap_co2 <<E-o-D
    &AFTITL
     CTITLE = '1234567890123456789012345678901234567890',
     CTITLE = 'PHC     monthly',
    /
    &AFFLAG
     ICTYPE =   1,
     KSIGMA =   ${KSIGMA},
     INTERP =   0,
     ITEST  =  40,
     JTEST  =  17,
     MONTH  = $MM,
    /
E-o-D
   #
   # --- Required Output, nitrate, phosphate, silicate, and oxygen
   #
   ${pput} fort.10  ${D}/dicr_sig${KSIGMA}_m${MM}.b
   ${pput} fort.10A ${D}/dicr_sig${KSIGMA}_m${MM}.a
   ${pput} fort.11  ${D}/alks_sig${KSIGMA}_m${MM}.b
   ${pput} fort.11A ${D}/alks_sig${KSIGMA}_m${MM}.a
   #
   # --- Optional Output.
   #
   #
   # --- end of month foreach loop
   #
   /bin/rm fort.7[123]
done
#
# --- Delete all files.
#
#/bin/rm -f *
