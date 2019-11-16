#!/bin/bash
#
# KAL - get X from input
if [ $# -ne 2 ] ; then
   echo "This script will set up the final relaxation files to be used by HYCOM."
   echo "Before this you should have set up the z-level files based on either PHC or "
   echo "Levitus climatologies."
   echo
   echo "You must input experiment name (example: 01.0) and climatology (phc or levitus) "
   echo "when running this script"
   exit
fi
export X=$1
export CLIM=$2



# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
#cd $(dirname $0)/../../ 
#export BASEDIR=$(pwd)/  
#cd -
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
   echo $BASEDIR
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }

# Check that pointer to HYCOM_ALL is set (from EXPT.src)
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
#
# --- Convert z-level climatology to HYCOM layers.
#
pget=cp
pput=cp

#
# --- E is experiment number, from EXPT.src
# --- R is region identifier, from EXPT.src
# --- T is topog. identifier, from EXPT.src
#
# --- P is primary path,
# --- S is scratch directory,
# --- D is permanent directory,
#
#
D=${BASEDIR}/relax/$E/   # Where data ends up
S=$D/SCRATCH              # SCRATCH dir 

echo $D
echo $S

mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}
KSIGMA=$(egrep "'thflag'"  ${BASEDIR}/expt_${X}/blkdat.input  | sed "s/.thflag.*$//" | tr -d "[:blank:]")
IVERSN=$(egrep "'iversn'"  ${BASEDIR}/expt_${X}/blkdat.input  | sed "s/.iversn.*$//" | tr -d "[:blank:]")
IVERSN=$(echo $IVERSN | sed "s/^0*//")
echo "IVERSN = $IVERSN"
echo "KSIGMA = $KSIGMA"
echo "SIGVER = $SIGVER"
tmp=$(expr $SIGVER % 2)
if [ $KSIGMA -eq 0 -a $tmp -ne 1 ] ; then
   echo "Recommend SIGVER=1,3,5 or 7 when thflag=0"
   exit 1
elif [ $KSIGMA -eq 2 -a $tmp -ne 0 ] ; then
   echo "Recommend SIGVER=2,4,6 or 8 when thflag=2"
   exit 1
fi

# Retrieve blkdat.input and create subset
touch blkdat.input
rm    blkdat.input
echo "World Ocean Atlas 2013 Climatology"              > blkdat.input
echo "  00        'month ' = month of climatology (01 to 12)"            >> blkdat.input
egrep "'iversn'|'iexpt '|'yrflag'" ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'idm   '|'jdm   '"          ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'itest '|'jtest '|'kdm   '" ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'nhybrd'|'nsigma'"          ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'dp00. '|'ds00. '"   ${BASEDIR}/expt_${X}/blkdat.input | egrep -v "'dp00i '"    >> blkdat.input
egrep "'thflag'|'thbase'|'sigma '" ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
egrep "'thkmin'"                   ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.input
if [ ! -s  blkdat.input ] ; then
   echo "Couldnt get blkdat.input " ; exit 1 ;
fi

touch   relax_tracer fort.51 fort.51A
/bin/rm relax_tracer fort.51 fort.51A

#
# --- 12 months
#
for MM in   01 02 03 04 05 06 07 08 09 10 11 12 ; do
   #
   # --- Input.
   #

   touch      fort.71 fort.71A fort.12 fort.12A 
   /bin/rm -f fort.71 fort.71A fort.12 fort.12A 
   ${pget} ${BASEDIR}/relax/$CLIM/dicr_sig${KSIGMA}_m${MM}.b fort.71  || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$CLIM/dicr_sig${KSIGMA}_m${MM}.a fort.71A || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$E/relax_int.b fort.12  || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$E/relax_int.a fort.12A || { echo "Couldnt get z-climatology" ; exit 1 ;}

   touch fort.51 fort.51A
   rm fort.51 fort.51A
   if [ ! -s fort.51 ] ; then
     ${pget} ${BASEDIR}/topo/depth_${R}_${T}.b fort.51 || { echo "Couldnt get depth_${R}_${T}.b " ; exit 1 ;}
   fi
   if [ ! -s  fort.51A ] ; then
     ${pget} ${BASEDIR}/topo/depth_${R}_${T}.a fort.51A  || { echo "Couldnt get depth_${R}_${T}.a " ; exit 1 ;}
   fi

   touch regional.grid.a regional.grid.b
   rm regional.grid.a regional.grid.b
   if [ ! -s regional.grid.b ]; then
     ${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b || { echo "Couldnt get regional.grid.b " ; exit 1 ;}
   fi
   if [ ! -s regional.grid.a ] ; then
     ${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a || { echo "Couldnt get regional.grid.b " ; exit 1 ;}
   fi
 
   
   touch relax_tracer
   if [ ! -s  relax_tracer ] ; then
    ${pget} ${HYCOM_ALL}/relax/src/relax_tracer .  || { echo "Couldnt get relax_tracer " ; exit 1 ;}

   fi
   wait
   chmod a+rx relax_tracer
   
   sed -e "s/^[ 	0-9]*'month ' =/  ${MM}	  'month ' =/" blkdat.input > fort.99
   
   /bin/rm -f core
   touch core
   export FOR010A=fort.10A
   export FOR012A=fort.12A
   export FOR021A=fort.21A
   export FOR051A=fort.51A
   export FOR071A=fort.71A
   /bin/rm -f fort.10  fort.21
   /bin/rm -f fort.10A fort.21A
   
   ./relax_tracer
   #
   # --- Output.
   #
   mv fort.10  relax_dic_m${MM}.b
   mv fort.10A relax_dic_m${MM}.a
   
   export DAYM=`echo ${MM} | awk '{printf("0000_%3.3d_00\n",30*($1-1)+16)}'`
   #export DAYM=`echo ${MM} | awk '{printf("0000_%3.3d_00\n",30.5*($1-1)+16)}'`
   #export DAYM=`echo ${MM} | awk '{printf("0000_%3.3d_00\n",30.5*($1-1)+1)}'`
   ${pput} fort.21  ${D}/relax.${DAYM}.b
   ${pput} fort.21A ${D}/relax.${DAYM}.a
   
   # --- end of month foreach loop
   /bin/rm fort.7[12]
done
#
# --- Merge monthly climatologies into one file.
#
cp relax_dic_m01.b relax.CO2_dic.b
#
for  MM  in  02 03 04 05 06 07 08 09 10 11 12 ; do
  tail -n +6 relax_dic_m${MM}.b >> relax.CO2_dic.b
done
#
cp relax_dic_m01.a relax.CO2_dic.a
#
for MM in  02 03 04 05 06 07 08 09 10 11 12 ; do
  cat relax_dic_m${MM}.a >> relax.CO2_dic.a
done
${pput} relax.CO2_dic.b ${D}/relax.CO2_dic.b
${pput} relax.CO2_dic.a ${D}/relax.CO2_dic.a
#
# --- delete the monthly files
#
/bin/rm relax_dic_m??.[ab]
#C
#C --- Delete all scratch directory files.
#C
#/bin/rm -f *

echo
echo "Finito... Final relaxation files should now be in directory"
echo "$D"

