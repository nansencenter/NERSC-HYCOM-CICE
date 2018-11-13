#!/bin/bash
#
## --- form the sst correction, -45*(EXPT-PF).
if [ $# -lt 7 ] ; then
   echo "This script will set up the offlux for HYCOM input sst [ab] files"
   echo "it takes the model long-term mean and also climotologyical SST data (or whatelse)"
   echo "To first calculate the annomally and then offlux. "
   echo
   echo "Example:"
   echo " ../bin/make_offlux.sh 800 760 -0.5 2.75 -45.0 model_sst_path clim_sst"
   echo " ../bin/make_offlux.sh idm jdm minc maxc offlux_corr model_sst_path clim_sst.a"
   echo "Will create offlux_${E}.[ab] and move it on relax/expt folder"
   echo
   echo " A. Wallcraft"
   echo " Mostafa Bakhoday-Paskyabi"
   echo " A. Alfatih "
   exit 1
else
#  sst_model_avg and sst_clim_avg must be a files
   idm=$1
   jdm=$2
   minlim=$3
   maxlim=$4
   off_corr=$5
   sst_model_avg=$6
   sst_clim_avg=$7
fi

ZERO=0.0
Y=00005


# Dont create region if one already exists
# Must be in expt or region dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
elif [ -f REGION.src ] ; then
   export BASEDIR=$(pwd)
else
   echo "Could not find EXPT.src or REGION.src. This script must be run in expt  opr region dir"
   exit 1
fi
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }

# create average fields
files=$(ls $data_files/archm*.a)

prog=${HYCOM_ALL}/bin/hycom_expr
chmod a+x $prog
$prog  $sst_model_avg $sst_clim_avg  $idm $jdm 1  -1 offlux_${E}_std.a | head -1 >> offlux_${E}_std.b


prog_clip=${HYCOM_ALL}/bin/hycom_clip
chmod a+x $prog_clip
$prog_clip  offlux_${E}_std.a $idm $jdm $minlim $maxlim offlux_${E}_std_cl.a  >> offlux_${E}_std_cl.b

# You may use following command insead of hycom_smooth
#prog_bin=${HYCOM_ALL}/bin/hycom_binning_fld
#chmod a+x $prog_bin
#$prog_bin offlux_${E}_std_cl.a 0.5 ../topo/regional.grid.a ../topo/regional.depth.a offlux_${E}_std_cl_bin.a #>! offlux_${E}_std_cl_bin.b

prog_bin=${HYCOM_ALL}/bin/hycom_smooth
chmod a+x $prog_bin
$prog_bin offlux_${E}_std_cl.a $idm $jdm   2  offlux_${E}_std_cl_bin.a # >! offlux_${E}_std_cl_bin.b

touch offlux_${E}.b && rm offlux_${E}.b
echo "Heat Flux Offset (W/m^2)"           >> offlux_${E}.b
echo "$off_corr * (model_sst-climo)"      >> offlux_${E}.b
echo "model mean from year ${Y}"          >> offlux_${E}.b
echo "anomaly between ${minlim} and ${maxlim} degC" >> offlux_${E}.b
echo "i/jdm = $idm $jdm"                  >> offlux_${E}.b


prog_expr=${HYCOM_ALL}/bin/hycom_expr
chmod a+x $prog_expr
touch offlux_${E}.a && rm offlux_${E}.a
$prog_expr offlux_${E}_std_cl_bin.a ONE $idm $jdm $off_corr $ZERO offlux_${E}.a  |  head -1 >> offlux_${E}.b
#
rm offlux_${E}_std*.*
#
prog_avg=${HYCOM_ALL}/bin/hycom_zonal_lat
chmod a+x $prog_avg
$prog_avg offlux_${E}.a 1 0.50 ../topo/regional.grid.a >! offlux_${E}.txt

mv offlux_${E}.a ../relax/${E}/offlux.a
mv offlux_${E}.b ../relax/${E}/offlux.b



