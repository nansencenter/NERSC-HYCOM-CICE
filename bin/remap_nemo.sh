#!/bin/bash

# Name: nemo2archvz.py
# Purpose: Convert MMERCATOR data to HYCOM archive files in z coordinates
# Author: Mostafa Bakhoday-Paskyabi (Mostafa.Bakhoday@nersc.no)
# Created: Sep 2017
# Copyright: (c) NERSC Norway 2017
# Licence:
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# Refrences:

#
#
# --- This code carries out 2D horizontal mapping and vertical re-gridding between outer
# --- (i.e. NEMO archive [ab] file) and inner regions.
# --- for a nesting scenario. Generally, the script contains following steps (covering both HYCOM2HYCOM and NEMO2HYCOM nesting):
# --- (1) 2D horizontal interpolation from the outer region to the inner region
#     (2) Vertical grid remapping from the original density layers to reference density layer listed in the "blkdat.input".
#
# --- Author: Mostafa Bakhoday-Paskyabi, Ocean Modeling group, NERSC, Bergen
# --- Mostafa.Bakhoday@nersc.no
# --- 8 December 2017.
#


iscan=15
usage="
   In order to use this routine, you need to first do mapping between outter and inner regions,
   using ../bin/isuba_gmapi.sh new region. Note that your current directory is in outer region
   experiment folder.

   Usage:
      $(basename $0) [-s iscan] NEMO mesh file  new_experiment_path MERCATOR netcdf files
   NEMO mesh file is located in a folder containing another NEMO-coordinate netcdf file with very similar 
   naming convention.

   Example:
      $(basename $0) /work/$USER/hycom/TP4a0.12/expt_03.1 archv.2013_003_12.a archv.2013_004.a
      $(basename $0) -s $iscan /work/$USER/hycom/TP4a0.12/expt_03.1 archv.2013_003_12.a

   Optional argument 'iscan' has default value of 15. This is the distance that will be scanned on
   this region grid to find  a sea point for the new region grid points.
"
# Must be in expt dir to run this script
#
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
export STARTDIR=$(pwd)
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
#
# Explore provided input path to get target region 
#
SIGVERSN=$SIGVER
thisexpt=$X
RR=$R
RT=$T
export RX=`echo ${E} | awk '{printf("%04.1f", $1*0.1)}'`

newexptpath=$(cd $1 && pwd)
newregionpath=$(cd  $newexptpath/.. && pwd)

echo "new experiment path $newexptpath"
echo "new region path $newregionpath"
source ${newregionpath}/REGION.src || { echo "Could not source ${newregionpath}/REGION.src" ; exit 1 ; }
source ${newexptpath}/EXPT.src || { echo "Could not source ${newexptpath}/EXPT.src" ; exit 1 ; }
echo "This region name    :$RR"

export U=$R
export E=$E
export X=`echo ${E} | awk '{printf("%04.1f", $1*0.1)}'`
#
export D=${BASEDIR}/expt_${RX}/data
export N=${BASEDIR}/subregion/${E}
export NEST=${newregionpath}/nest/${E}
export EXPTNEW=${newregionpath}/expt_${X}
#
mkdir -p $N
mkdir -p $NEST
#
# Change current directory to the outer region data directory
#
cd $D
echo ${BASEDIR}/topo/depth_${RR}_${RT}.b
echo ${newregionpath}/topo/depth_${U}_${T}.b
echo "################"
#
#   From outer region
#
cp  ${BASEDIR}/topo/regional.grid.a           $D/regional.grid.a
cp  ${BASEDIR}/topo/regional.grid.b           $D/regional.grid.b
cp  ${BASEDIR}/topo/depth_${RR}_${RT}.a       $D/regional.depth.a
cp  ${BASEDIR}/topo/depth_${RR}_${RT}.b       $D/regional.depth.b
#
#   From inner region
#
cp  ${BASEDIR}/topo/ZL50.txt                  $N/ZL50.txt
cp  ${BASEDIR}/subregion/${U}.gmap.a          $N/regional.gmap.a
cp  ${BASEDIR}/subregion/${U}.gmap.b          $N/regional.gmap.b
cp  ${EXPTNEW}/blkdat.input                   $N/blkdat.input
cp  ${newregionpath}/topo/regional.grid.a     $N/regional.grid.a
cp  ${newregionpath}/topo/regional.grid.b     $N/regional.grid.b
cp  ${newregionpath}/topo/depth_${U}_${T}.a   $N/regional.depth.a
cp  ${newregionpath}/topo/depth_${U}_${T}.b   $N/regional.depth.b
#
#   idm, jdm, & kdm from target xperiment blkdat.input
#
target_idm=$(blkdat_get $N/regional.grid.b idm)
target_jdm=$(blkdat_get $N/regional.grid.b jdm)
kdmold=$(blkdat_get ${EXPTNEW}/blkdat.input kdm)
kdmnew=${kdmold}
#
# define number of layered archives in subregion folder  before
# applying triming
#
export L="_L${target_kdm}"
#
prog_subreg=${HYCOM_ALL}/subregion/src/isubaregion
prog_nemo=${HYCOM_ALL}/relax/src/nemo_archvz
#
chmod a+rx ${prog_subreg}
chmod a+rx ${prog_nemo}
#
echo
source_archv_i=$2

if [ ! -f ${source_archv_i}.a -o ! -f ${source_archv_i}.b ] ; then
    echo "Source file ${source_archv_i}.[ab] does not exist"
    continue
fi
target_archv=${source_archv_i}
touch ${N}/${target_archv}${L}.a && rm ${N}/${target_archv}${L}.a
touch ${N}/${target_archv}${L}.b && rm ${N}/${target_archv}${L}.b
    
logfile=${N}/isubaregion.log
touch $logfile && rm $logfile
#
# --- 'flnm_reg'  = target sub-region grid       filename
# --- 'flnm_map'  = target sub-region grid map   filename
# --- 'flnm_top'  = target bathymetry filename, or 'NONE'
# --- 'flnm_tin'  = input  bathymetry filename, or 'NONE'
# --- 'flnm_in'   = input  archive    filename
# --- 'flnm_out'  = output archive    filename
# --- 'cline_out' = output title line (replaces preambl(5))
#
echo ${N}/regional.depth.a
echo "Processing ${N}/${target_archv}"

${prog_subreg}   <<EOF
${N}/regional.grid.a
${N}/regional.gmap.a
${N}/regional.depth.a
${D}/regional.depth.a
${D}/${source_archv_i}.a
${N}/${target_archv}${L}.a
${RR} interpolated to ${U}
${target_idm}	  'idm   ' = target longitudinal array size
${target_jdm}	  'jdm   ' = target latitudinal  array size
0	        'iceflg' = ice in output archive flag (0=none,1=energy loan model)
0	        'smooth' = smooth interface depths    (0=F,1=T)
$iscan      'iscan'
EOF
#
echo "Log can be found in $logfile"
echo "    ($h) *** prepare for nesting folder"
if [ ! -f ${N}/${target_archv}${L}.b -o ! -f ${N}/${target_archv}${L}.a ]; then
    echo "missing archive file: " ${N}/${target_archv}${L}
    exit 2
fi
#
# --- change current directory to nesting folder of inner region
#
cd $N
touch   ${target_archv}.a && rm ${target_archv}.a
touch   ${target_archv}.b && rm ${target_archv}.b
logfile=${N}/nemo_archv.log
touch $logfile && rm $logfile



# Retrieve blkdat.input and create subset
touch blkdat.subset
rm    blkdat.subset
#echo "" > blkdat.subset
echo "NEMO Relaxation fields"                              > blkdat.subset
echo "  $SIGVERSN        'sigver ' = Version of eqn of state  "                        >> blkdat.subset
echo "  1        'levtop ' = top level of input clim. to use (optional, default 1)"  >> blkdat.subset
egrep "'iversn'"    blkdat.input >> blkdat.subset
egrep "'iexpt '"    blkdat.input >> blkdat.subset
egrep "'mapflg'"    blkdat.input >> blkdat.subset
egrep "'yrflag'"    blkdat.input >> blkdat.subset
egrep "'idm   '"    blkdat.input >> blkdat.subset
egrep "'jdm   '"    blkdat.input >> blkdat.subset
echo "  0        'jdw    ' = width of zonal average (optional, default 0)"  >> blkdat.subset
echo "  -1       'itest  ' = grid point where detailed diagnostics are desired"  >> blkdat.subset
echo "  -1       'jtest  ' = grid point where detailed diagnostics are desired"  >> blkdat.subset
egrep "'kdm   '"    blkdat.input >> blkdat.subset
egrep "'nhybrd'"    blkdat.input >> blkdat.subset
egrep "'nsigma'"    blkdat.input >> blkdat.subset
egrep "'isotop'"    blkdat.input >> blkdat.subset
egrep "'dp00  '"    blkdat.input >> blkdat.subset
egrep "'dp00x '"    blkdat.input >> blkdat.subset
egrep "'dp00f '"    blkdat.input >> blkdat.subset
egrep "'ds00  '"    blkdat.input >> blkdat.subset
egrep "'ds00x '"    blkdat.input >> blkdat.subset
egrep "'ds00f '"    blkdat.input >> blkdat.subset
egrep "'ds0k  '"    blkdat.input >> blkdat.subset
egrep "'dp0k  '"    blkdat.input >> blkdat.subset
egrep "'dp00i '"    blkdat.input >> blkdat.subset
egrep "'thflag'"    blkdat.input >> blkdat.subset
egrep "'thbase'"    blkdat.input >> blkdat.subset
egrep "'vsigma'"    blkdat.input >> blkdat.subset
egrep "'sigma '"    blkdat.input >> blkdat.subset
egrep "'thkmin'"    blkdat.input >> blkdat.subset
if [ ! -s  blkdat.subset ] ; then
echo "Couldnt get blkdat.input " ; exit 1 ;
fi



rm -rf fort.*
mv blkdat.subset fort.99

#
# --- 'flnm_i' = name of source archive file
# --- 'flnm_o' = name of target  archive file
# --- 'NZ ' = Number of vertical z levels (for MERCATOR sets to 50)
# --- 'flnm_z' = text file containing the vertical grids in z-level
# --- 'flag_t   ' flag to activate temperature for vertical interpolation
# --- 'flag_s   ' flag to activate salinity for vertical interpolation
# --- 'flag_th   ' flag to activate thickness for vertical interpolation
# --- 'flag_u   ' flag to activate u-componen of current for vertical interpolation
# --- 'flag_v   ' flag to activate v-component of current for vertical interpolation
#
touch ${NEST}/${target_archv}.a
touch ${NEST}/${target_archv}.b
rm -rf ${NEST}/${target_archv}.*
# --- 'flnm_i' = name of source archive file
# --- 'flnm_o' = name of target  archive file
# --- 'NZ ' = Number of vertical z levels (for MERCATOR sets to 50)
# --- 'flnm_z' = text file containing the vertical grids in z-level
# --- 'flag_t   ' flag to activate temperature for vertical interpolation
# --- 'flag_s   ' flag to activate salinity for vertical interpolation
# --- 'flag_th   ' flag to activate thickness for vertical interpolation
# --- 'flag_u   ' flag to activate u-componen of current for vertical interpolation
# --- 'flag_v   ' flag to activate v-component of current for vertical interpolation
#
${prog_nemo}  <<EOF
${N}/${target_archv}${L}.a
${NEST}/${target_archv}.a
50
${N}/ZL50.txt
T
T
T
T
T
EOF

cd $D
echo
echo "Succesfully created archive file: $2"
echo
echo
#
