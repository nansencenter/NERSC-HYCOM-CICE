#!/bin/bash
#
#
# --- Generate a thkdf4 file (horizaontal diffusion) suitable for hycom
#

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BASEDIR=$(cd .. && pwd)/ 
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src            || { echo "Could not source ${BASEDIR}/expt_$X/EXPT.src" ; exit 1 ; }



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
D=$BASEDIR/relax/$E/
C=$D/SCRATCH
mkdir -p $C
cd       $C || { echo " Could not descend scratch dir $C" ; exit 1;}

touch   regional.grid.a regional.grid.b regional.depth.a regional.depth.b
/bin/rm regional.grid.a regional.grid.b regional.depth.a regional.depth.b
#
# --- Input.
#
touch fort.51 fort.51A
rm    fort.51 fort.51A
${pget} ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Cant retrieve depth_${R}_${T}.a" ; exit 1; }
${pget} ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Cant retrieve depth_${R}_${T}.b" ; exit 1; }

#
touch regional.grid.a regional.grid.b
rm   regional.grid.a regional.grid.b
${pget} ${BASEDIR}/topo/regional.grid.b regional.grid.b  || { echo "Cant retrieve regional_grid.b" ; exit 1; }
${pget} ${BASEDIR}/topo/regional.grid.a regional.grid.a  || { echo "Cant retrieve regional_grid.a" ; exit 1; }



$MSCPROGS/bin/thkdf4-2.2.37


#
#
# --- Output.
#
${pput} thkdf4.a ${D}/
${pput} thkdf4.b ${D}/

echo
echo "thkdf4.[ab] copied to $D"
echo "You should check that its ok by looking at it with matlab or with plotfp.sh"
