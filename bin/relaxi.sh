#!/bin/bash
#
# KAL - get X from input
if [ $# -ne 1 ] ; then
   echo "This script will set up the final relaxation files to be used by HYCOM."
   echo "Before this you should have set up the z-level files based on either PHC or "
   echo "Levitus climatologies."
   echo
   echo "You must input climatology (phc or levitus)"
   echo "when running this script"
   echo ""
   echo "Example:"
   echo "   $0 woa2013"
   exit 1
fi
export CLIM=$1


# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }

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
mkdir -p $S
cd       $S || { echo " Could not descend scratch dir $S" ; exit 1;}

#
# --- Existing files must be deleted first
#
proceed=1
if [ -f $D/relax_sal.a -a -f $D/relax_sal.b ] ; 
then
   echo "Salinity   relaxation files already exist : $D/relax_sal.[ab]"
   proceed=0
fi
if [ -f $D/relax_tem.a -a -f $D/relax_tem.b ] ; 
then
   echo "Temperature relaxation files already exist: $D/relax_tem.[ab]"
   proceed=0
fi
if [ -f $D/relax_int.a -a -f $D/relax_int.b ] ; 
then
   echo "Interface   relaxation files already exist: $D/relax_int.[ab]"
   proceed=0
fi
if [ $proceed -eq 0 ] ; then
   echo "one or more relaxation files already present in $D "
   echo "Delete these first if you want to regenerate the relax files"
   exit 0
fi
   
   


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
touch blkdat.subset
rm    blkdat.subset
echo "$CLIM Climatology"                              > blkdat.subset
echo "  00        'month ' = month of climatology (01 to 12)"                        >> blkdat.subset
echo "  $SIGVER        'sigver ' = Version of eqn of state  "                        >> blkdat.subset
echo "  1        'levtop ' = top level of input clim. to use (optional, default 1)"  >> blkdat.subset
egrep "'iversn'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'iexpt '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'mapflg'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'yrflag'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'idm   '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'jdm   '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
echo "  0        'jdw    ' = width of zonal average (optional, default 0)"  >> blkdat.subset
echo "  -1       'itest  ' = grid point where detailed diagnostics are desired"  >> blkdat.subset
echo "  -1       'jtest  ' = grid point where detailed diagnostics are desired"  >> blkdat.subset
egrep "'kdm   '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'nhybrd'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'nsigma'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
if [ $IVERSN -lt 20 ] ; then
   egrep "'dp00s '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp00  '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp00x '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp00f '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
else 
   egrep "'isotop'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp00  '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp00x '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp00f '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'ds00  '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'ds00x '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'ds00f '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'ds0k  '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp0k  '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
   egrep "'dp00i '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
fi
egrep "'thflag'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'thbase'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'vsigma'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'sigma '"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
egrep "'thkmin'"  ${BASEDIR}/expt_${X}/blkdat.input >> blkdat.subset
if [ ! -s  blkdat.subset ] ; then
   echo "Couldnt get blkdat.input " ; exit 1 ;
fi

#get  thflag from blkdat.input
touch   relaxi fort.51 fort.51A
/bin/rm relaxi fort.51 fort.51A

#
# --- 12 months
#
for MM in   01 02 03 04 05 06 07 08 09 10 11 12 ; do
   #
   # --- Input.
   #

   touch      fort.71 fort.71A fort.72 fort.72A
   /bin/rm -f fort.71 fort.71A fort.72 fort.72A
   ${pget} ${BASEDIR}/relax/$CLIM/dens_sig${KSIGMA}_m${MM}.b fort.71  || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$CLIM/dens_sig${KSIGMA}_m${MM}.a fort.71A || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$CLIM/temp_sig${KSIGMA}_m${MM}.b fort.72  || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$CLIM/temp_sig${KSIGMA}_m${MM}.a fort.72A || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$CLIM/saln_sig${KSIGMA}_m${MM}.b fort.73  || { echo "Couldnt get z-climatology" ; exit 1 ;}
   ${pget} ${BASEDIR}/relax/$CLIM/saln_sig${KSIGMA}_m${MM}.a fort.73A || { echo "Couldnt get z-climatology" ; exit 1 ;}

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
 
   touch relaxi
   if [ ! -s  relaxi ] ; then
    ${pget} ${HYCOM_ALL}/relax/src/relaxi .  || { echo "Couldnt get relaxi " ; exit 1 ;}

   fi
   wait
   chmod a+rx relaxi

   # Create subset of blkdat.input
   sed -e "s/^[ 	0-9]*'month ' =/  ${MM}	  'month ' =/" blkdat.subset > fort.99
   
   /bin/rm -f core
   touch core
   export FOR010A=fort.10A
   export FOR011A=fort.11A
   export FOR012A=fort.12A
   export FOR021A=fort.21A
   export FOR051A=fort.51A
   export FOR071A=fort.71A
   export FOR072A=fort.72A
   export FOR073A=fort.73A
   /bin/rm -f fort.10  fort.11  fort.12  fort.21
   /bin/rm -f fort.10A fort.11A fort.12A fort.21A
 
   logfile=$S/log_$MM.out
   echo "Running relaxi, log in $logfile"
   ./relaxi > $logfile 2>&1
   #
   # --- Output.
   #
   mv fort.10  relax_tem_m${MM}.b
   mv fort.10A relax_tem_m${MM}.a
   mv fort.11  relax_sal_m${MM}.b
   mv fort.11A relax_sal_m${MM}.a
   mv fort.12  relax_int_m${MM}.b
   mv fort.12A relax_int_m${MM}.a
   
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
cp relax_int_m01.b relax_int.b
cp relax_sal_m01.b relax_sal.b
cp relax_tem_m01.b relax_tem.b
#
for  MM  in  02 03 04 05 06 07 08 09 10 11 12 ; do
  tail +6 relax_int_m${MM}.b >> relax_int.b
  tail +6 relax_sal_m${MM}.b >> relax_sal.b
  tail +6 relax_tem_m${MM}.b >> relax_tem.b
done
#
cp relax_int_m01.a relax_int.a
cp relax_sal_m01.a relax_sal.a
cp relax_tem_m01.a relax_tem.a
#
for MM in  02 03 04 05 06 07 08 09 10 11 12 ; do
  cat relax_int_m${MM}.a >> relax_int.a
  cat relax_sal_m${MM}.a >> relax_sal.a
  cat relax_tem_m${MM}.a >> relax_tem.a
done
${pput} relax_int.b ${D}/relax_int.b
${pput} relax_int.a ${D}/relax_int.a
${pput} relax_sal.b ${D}/relax_sal.b
${pput} relax_sal.a ${D}/relax_sal.a
${pput} relax_tem.b ${D}/relax_tem.b
${pput} relax_tem.a ${D}/relax_tem.a
#
# --- delete the monthly files
#
/bin/rm relax_int_m??.[ab]
/bin/rm relax_sal_m??.[ab]
/bin/rm relax_tem_m??.[ab]
echo "Final files in $D"




#KAL - move to separate routine
##C
##C --- Delete all scratch directory files.
##C
##/bin/rm -f *
## Check that pointer to MSCPROGS is set (from EXPT.src)
#if [ -z ${MSCPROGS} ] ; then
#   echo "MSCPROGS Environment not set "
#   exit
#else
#   if [ ! -d ${MSCPROGS} ] ; then
#      echo "MSCPROGS not properly set up"
#      echo "MSCPROGS not a directory at ${MSCPROGS}"
#      exit
#   fi
#fi
## Create the spatially varying threshold sor salinity above which we do not relax (currently 0.5psu)
#cd ${D}
#${MSCPROGS}/src/Model_input-2.2.37/sssrmx-2.2.37
#
#echo
#echo "Finito... Final relaxation files should now be in directory"
#echo "$D"
#
