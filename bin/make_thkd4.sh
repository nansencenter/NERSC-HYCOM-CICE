#! /bin/bash
# Must be in expt dir to run this script
# MBP
#
if [ -f EXPT.src ] ; then
    export BASEDIR=$(cd .. && pwd)
else
    echo "Could not find EXPT.src. This script must be run in expt dir"
    exit 1
fi
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
#
expt_path=${BASEDIR}/expt_$X
topo_path=${BASEDIR}/topo
data_path=${expt_path}/../relax/${E}/
#
export SCRATCH=$expt_path/SCRATCH


prog=$MSCPROGS/bin/thkdf4-2.2.37


cp $topo_path/depth_${R}_$T.a regional.depth.a || { echo "Cannot access to depth_${R}_${T}.a" ; exit 1; }
cp $topo_path/depth_${R}_$T.b regional.depth.b || { echo "Cannot access to depth_${R}_${T}.b" ; exit 1; }
cp $topo_path/regional.grid.* .                || { echo "Cannot access to regional.grid.[ab] files" ; exit 1; }

chmod a+x $prog
logfile=$SCRATCH/thkd4.log

$prog

mv -f thkdf4.* ${data_path}
