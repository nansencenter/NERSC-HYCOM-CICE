#!/bin/bash
#
#set -x
pget=cp
pput=cp
#
# --- Create HYCOM kpar files.
#
# --- P is primary path,
# --- S is scratch directory,
# --- D is permanent directory,
#


# Set basedir based on relative paths of script
# Can be troublesome, but should be less prone to errors
# than setting basedir directly
export BASEDIR=$(cd  $(dirname $0)/../ && pwd )/
source ${BASEDIR}/bin/common_functions.sh || { echo "Could not source ${BASEDIR}/bin/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
D=$BASEDIR/force/seawifs/
S=$D/SCRATCH
mkdir -p $S
cd       $S

# Check that pointer to HYCOM_ALL is set
if [ -z ${HYCOM_ALL} ] ; then
   echo "HYCOM_ALL Environment not set "
   exit
else
   if [ ! -d ${HYCOM_ALL} ] ; then
      echo "HYCOM_ALL not properly set up"
      echo "HYCOM_ALL not a directory at ${HYCOM_ALL}"
      exit
   fi
fi
echo "jau"


#
# --- Input.
#
copy_grid_files $S
touch      fort.71  && /bin/rm -f fort.71
${pget}  ${SEAWIFS}/seawifs_ann_kpar.d fort.71 || { echo "Cant get seawifs data at   ${SEAWIFS}" ; exit 1 ; }
cp ${HYCOM_ALL}/force/src/kp . ||  { echo "Cant get kp routine" ; exit 1 ; }
wait
chmod a+rx kp

/bin/rm -f fort.1[01234]*
export FOR010A=fort.10A
/bin/rm -f core
touch core
#
./kp <<E-o-D
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'SeaWifs 1997:2001 monthly clim, 1/m',
 /
 &AFTIME
  PARMIN =    0.04, !minimum kpar (default 0.04)
  PARMAX =    0.20  !maximum kpar (default 0.20)
 /
 &AFFLAG
  INTERP =    0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =    0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 /
E-o-D
#
#
# --- Output.
#
mv fort.10  kpar.b
mv fort.10A kpar.a
${pput} kpar.b ${D}/kpar.b
${pput} kpar.a ${D}/kpar.a
#
# --- Delete all files.
#
#/bin/rm -f *
echo "kpar.[ab] placed in $D"
