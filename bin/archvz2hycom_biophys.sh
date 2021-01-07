#!/bin/bash
#
# --- This code mapping 2D horizontal mapping between outer and inner regions
# --- for a nesting scenario. The script use two steps:
# --- (1) 2D horizontal interpolation from the outer region to inner region
# --- (2) vertical trimming and discarding of unnecessary fields from (1)
# ---     in order to be used in nesting.
#
# --- Author: Mostafa Bakhoday-Paskyabi, Ocean Modeling group, NERSC, Bergen
# --- Mostafa.Bakhoday@nersc.no
# --- M. Bakhoday-Paskyabi 9 September 2017.
# --- M. Bakhoday-Paskyabi 11 July 2019.
#


iscan=15
usage="
   In order to use this routine, you need to first do mapping between outter and inner regions,
   using ../bin/isuba_gmapi.sh new region. Note that your current directory is in outer region
   experiment folder.

   Usage:
      $(basename $0) [-s iscan]  new_experiment_path archive1 archive2 ....

   Example:
      $(basename $0) /work/$USER/hycom/TP4a0.12/expt_03.1 archv.2013_003_12.a archv.2013_004.a
      $(basename $0) -s $iscan /work/$USER/hycom/TP4a0.12/expt_03.1 archv.2013_003_12.a

   Optional argument 'iscan' has default value of 15. This is the distance that will be scanned on
   this region grid to find  a sea point for the new region grid points.
"
grid_type=native
bio_flag=0

# This will process optional arguments
options=$(getopt -o b:m -- "$@")
[ $? -eq 0 ] || {
    echo "$usage"
    echo "Error: Incorrect options provided"
    exit 1
}
eval set -- "$options"
while true; do
    case "$1" in
    -b)
        shift;
        bio_flag=1
        ;;
    -m)
        grid_type=regular
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done
echo $grid_type
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
T=$T
export X=`echo ${E} | awk '{printf("%04.1f", $1*0.1)}'`
#
export D=${BASEDIR}/expt_${RX}/data
export N=${BASEDIR}/subregion/${E}
export NEST=${newregionpath}/nest/${E}
export EXPTNEW=${newregionpath}/expt_${X}
#
mkdir -p $N
mkdir -p $NEST
mkdir -p $D
#
# Change current directory to the outer region data directory
#
cd $D
echo ${BASEDIR}/topo/depth_${RR}_${RT}.b
echo ${newregionpath}/topo/depth_${U}_${T}.b
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
echo "Inner domain topo: depth_${U}_${T}.b" 
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
if [ "${grid_type}" == "native"  ] ; then
   prog_subreg=${HYCOM_ALL}/subregion/src/isubaregion_modified
   prog_nemo=${HYCOM_ALL}/relax/src/nemo_archvz_biophys
else
#   prog_subreg=${HYCOM_ALL}/subregion/src/isubaregion_modified
   prog_subreg=${HYCOM_ALL}/subregion/src/isubaregion
   prog_nemo=${HYCOM_ALL}/relax/src/nemo_archvz_biophys
fi
#
#
#
export hinterp_method=1
export vinterp_method=1
#
#
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

if [ ${hinterp_method} -eq 1 ]; then

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

##${prog_subreg}   <<EOF
${prog_subreg} >> $logfile  <<EOF
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

fi # hinterp_method

#
# --- change current directory to nesting folder of inner region
#

if [ ${vinterp_method} -eq 1 ]; then


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

logfile=${NEST}/nest_archvz.log
touch $logfile && rm $logfile

echo $logfile
touch ${NEST}/${target_archv}.a
touch ${NEST}/${target_archv}.b
rm -rf ${NEST}/${target_archv}.*


if [[ "${bio_flag}" -eq 0 ]] ; then

${prog_nemo}  >> $logfile  <<EOF
${N}/${target_archv}${L}.a
${NEST}/${target_archv}.a
50
${N}/ZL50.txt
T
T
T
T
T
NONE
NONE
NONE
EOF
else
test1=${target_archv:0:5}
test2=${target_archv:5:16}
target_bioarchv=${test1}_fabm${test2}
echo ${NEST}/${target_bioarchv}.a
echo ${NEST}/${target_archv}.a
touch ${NEST}/${target_bioarchv}.a
touch ${NEST}/${target_bioarchv}.b
rm -rf ${NEST}/${target_bioarchv}.*

${prog_nemo}  >> $logfile  <<EOF
${N}/${target_archv}${L}.a
${NEST}/${target_archv}.a
50
${N}/ZL50.txt
T
T
T
T
T
T
T
T
EOF
fi


fi

cd $D
echo
if [ ! -f ${NEST}/${target_archv}.a -o ! -f ${NEST}/${target_archv}.b ] ; then
    echo "Source file ${NEST}/${target_archv}.[ab] do not exist"
else
    echo "Removing all temporary [ab] files ..."
    touch ${N}/${target_archv}${L}.a
    touch ${N}/${target_archv}${L}.b
    touch ${D}/${source_archv_i}.a
    touch ${D}/${source_archv_i}.b
    rm -rf ${N}/${target_archv}${L}.*
    rm -rf ${D}/${source_archv_i}.*
    echo "Succesfully created archive file: $2"

    cd $newexptpath
   ## python ${BINDIR}/fix_montg1.py ${NEST}/${target_archv}.a ${INPUTDIR}montg_regress.pckl ./ 
    python ${BINDIR}/fix_montg1.py ${NEST}/${target_archv}.a ${INPUTDIR}TP2_montg_regress.pckl ./ 
    cd $BASEDIR
fi
echo
#

