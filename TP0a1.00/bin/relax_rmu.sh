#!/bin/bash
#
#
# --- Generate a HYCOM relaxation mask.
#
# KAL - get X from input
if [ $# -ne 1 ] ; then
   echo "You must input experiment name (ex 01.0)"
   exit
fi
export X=$1

# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd $(dirname $0)/../ && pwd)/
source ${BASEDIR}/bin/common_functions.sh || { echo "Could not source ${BASEDIR}/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$X/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }

D=$BASEDIR/relax/$E/
C=$D/SCRATCH
mkdir -p $C
cd       $C || { echo " Could not descend scratch dir $C" ; exit 1;}



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



#set -x
pget=cp
pput=cp
#
# --- E is experiment number, from EXPT.src
# --- R is region identifier, from EXPT.src
# --- T is topog. identifier, from EXPT.src
#
# --- P is primary path,
# --- C is scratch directory,
# --- D is permanent directory,
#
#

touch   rmu fort.51 fort.51A regional.grid.a regional.grid.b regional.depth.a regional.depth.b
/bin/rm rmu fort.51 fort.51A regional.grid.a regional.grid.b regional.depth.a regional.depth.b
#
# --- Input.
#
touch fort.51 fort.51A
rm    fort.51 fort.51A
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b fort.51  || { echo "Cant retrieve depth_${R}_${T}.b" ; exit 1; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a fort.51A || { echo "Cant retrieve depth_${R}_${T}.b" ; exit 1; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Cant retrieve depth_${R}_${T}.a" ; exit 1; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Cant retrieve depth_${R}_${T}.b" ; exit 1; }

#
touch regional.grid.a regional.grid.b
rm   regional.grid.a regional.grid.b
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b  || { echo "Cant retrieve regional_grid.b" ; exit 1; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a  || { echo "Cant retrieve regional_grid.a" ; exit 1; }

touch rmu 
rm    rmu 
cp      ${HYCOM_ALL}/relax/src/rmu . || { echo "Cant retrieve rmu from HYCOM_ALL dir" ; exit 1; }

#
/bin/rm -f core
touch core
/bin/rm -f fort.21 fort.21A
#
export FOR021A=fort.21A
export FOR051A=fort.51A
#
# KAL - must place this in 010 ?
#./rmu <<E-o-D
# &MASK
#  CTITLE = 'N and S boundary: 5 grid points with 20 to 120 day e-folding time',
#  IF     =   1,     1,     1,     1,     1,
#             1,     1,     1,     1,     1,
#  IL     =  56,    56,    56,    56,    56,
#            56,    56,    56,    56,    56,
#  JF     =   1,     2,     3,     4,     5,
#            51,    50,    49,    48,    47,
#  JL     =   1,     2,     3,     4,     5,
#            51,    50,    49,    48,    47,
#  EFOLD  =  20.0,  45.0,  70.0,  95.0, 120.0,
#            20.0,  45.0,  70.0,  95.0, 120.0,
# /
#E-o-D

# KAL use routine "rmunew" in stead - Option 1: et input from user 
# - this only activates the "simple" option for relaxation along
# walls



echo "Do you want to include port in your model that is not on a relaxed model boundary (T/F)"
read port
if [ $port == "T" ]
 then
      rm rmu.in
      cp  ${BASEDIR}/expt_$X/rmu.in .
      [ ! -f "rmu.in" ] && echo "I cant find file \"rmu.in\" ... This file is needed if you " && \
      echo " want to use relaxation in the port... A sample file can be found in $INFILEDIR" && exit
else
   echo 
   echo "Now input where you want relaxation along sidewalls"
   echo "There will be no sidewall relaxation activated if you "
   echo "have land covering 2 and more grid cells from the edges"
   echo "Input T (true) or F(false) below, and a relaxation width"
   echo
   echo "Do you want relaxation along Eastern wall ? (T/F)"
   read least
   echo "Do you want relaxation along Western wall ? (T/F)"
   read lwest
   echo "Do you want relaxation along Southern wall? (T/F)"
   read lsouth
   echo "Do you want relaxation along Northern wall? (T/F)"
   read lnorth
   echo "Input relaxation width along walls (typically 20)"
   read rwidth
   echo "Input relaxation time scale maximum  in days (typically 20)"
   echo "This time scale will be reduced linearly towards zero away from the wall"
   read rscale

   echo $least  >  rmu.in
   echo $lwest  >> rmu.in
   echo $lsouth >> rmu.in
   echo $lnorth >> rmu.in
   echo $rwidth >> rmu.in
   echo $rwidth >> rmu.in
   echo $rscale >> rmu.in
fi

$MSCPROGS/bin_setup/rmunew
#knutali/Progs/HYCOM_UTILITY_ROUTINES/HYCOM_2.1.03/ALL/relax/src_2.1.00/rmunew


#
#
# --- Output.
#
${pput} fort.21  ${D}/relax_rmu.b
${pput} fort.21A ${D}/relax_rmu.a

echo
echo "relax_rmu.[ab] copied to $D"
echo "You should check that its ok by looking at it with matlab or with plotfp.sh"
echo "The values in that file are 1./(timescale*86400)"
