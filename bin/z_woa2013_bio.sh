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
L=$WOA2013BIO_PATH/
[  "${L}" == ""  ] && { echo "WOA2013_PATH is not set " ; exit 1; }
[ ! -d $L ] && { echo "Cant find WOA2013_PATH directory $L" ; exit 1; }

mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}


touch   z_woa2013_bio
/bin/rm z_woa2013_bio
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
   touch      fort.71 fort.72 fort.73 fort.74
   /bin/rm -f fort.71 fort.72 fort.73 fort.74
   ${pget} $L/n${MM}.d fort.71 || { echo "Cant copy $L/n${MM}.d "; exit 1 ; }
   ${pget} $L/p${MM}.d fort.72 || { echo "Cant copy $L/p${MM}.d "; exit 1 ; }
   ${pget} $L/i${MM}.d fort.73 || { echo "Cant copy $L/i${MM}.d "; exit 1 ; }
   ${pget} $L/o${MM}.d fort.74 || { echo "Cant copy $L/o${MM}.d "; exit 1 ; }

   if [ ! -s regional.grid.b ] ; then
     ${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Cant copy regional.grid.b "; exit 1 ; }
   fi
   if [ ! -s regional.grid.a ] ; then
     ${pget} ${BASEDIR}//topo/regional.grid.a regional.grid.a || { echo "Cant copy regional.grid.a "; exit 1 ; }
   fi
   touch z_woa2013_bio
   if [ ! -s  z_woa2013_bio ] ; then
     ${pget} ${MSCPROGS}/src/Relax/z_woa2013_bio .  || { echo "Cant copy z_woa2013_bio "; exit 1 ; }
   fi
   wait

   chmod a+rx z_woa2013_bio
   /bin/rm -f core
   touch core
   export FOR010A=fort.10A
   export FOR011A=fort.11A
   export FOR012A=fort.12A
   export FOR013A=fort.13A
   /bin/rm -f fort.10 fort.10A fort.11 fort.11A fort.12 fort.12A fort.13 fort.13A
   ./z_woa2013_bio <<E-o-D
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
   ${pput} fort.10  ${D}/nitr_sig${KSIGMA}_m${MM}.b
   ${pput} fort.10A ${D}/nitr_sig${KSIGMA}_m${MM}.a
   ${pput} fort.11  ${D}/phos_sig${KSIGMA}_m${MM}.b
   ${pput} fort.11A ${D}/phos_sig${KSIGMA}_m${MM}.a
   ${pput} fort.12  ${D}/sili_sig${KSIGMA}_m${MM}.b
   ${pput} fort.12A ${D}/sili_sig${KSIGMA}_m${MM}.a
   ${pput} fort.13  ${D}/oxyg_sig${KSIGMA}_m${MM}.b
   ${pput} fort.13A ${D}/oxyg_sig${KSIGMA}_m${MM}.a
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
