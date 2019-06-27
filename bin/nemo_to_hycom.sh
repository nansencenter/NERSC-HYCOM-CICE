#! /bin/bash
# Name: make_nemo_archives.py
# Purpose: Convert MMERCATOR data to HYCOM archive files in isopycnal coordinates
# Author: Mostafa Bakhoday-Paskyabi (Mostafa.Bakhoday@nersc.no)
# Created: October 2017
# Copyright: (c) NERSC Norway 2017
# Licence:
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# History
# (1) December 2018: Mostafa Bakhoday-Paskyabi
# (2) May 2019: some correction on biophys [ab] files, Mostafa Bakhoday-Paskyabi
#
#
module load basemap/1.0.7-intel-2017a-Python-2.7.13

options=$(getopt -o b:m  -- "$@")
maxinc=50
bio_path=""
grid_type=native
eval set -- "$options"
while true; do
    case "$1" in
    -b)
       shift;
       bio_path=$1
       #grid_type=$1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done

if [ $# -lt 2 ] ; then
    echo "This script will set up the final nesting files from MERCATOR 1/12 degree to be used by HYCOM."
    echo "The code contains the following steps:"
    echo "(1) Create archive [ab] files from the MERCATOR netcdf file."
    echo "    The final archive files includes all 2D fields (filled with zero except baraotropic velocities)."
    echo "    and 3D fields of temperaure, salinity, thickness, and two components of velocities."
    echo "(2) Based on generated archive files in (1) the grid and topography files are generated."
    echo "    Note that the grids are non-native and interpolated into a rectilinear mercator grids horizontally."
    echo "(3) "
    echo "You must input climatology (phc or levitus)"
    echo "when running this script"
    echo ""
    echo "Example:"
    echo "   ../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_01.0/ ../MERCATOR-PHY-24-2011-01-02*.nc"
    echo " ../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_01.2/  /nird/projects/nird/NS9481K/MERCATOR_DATA/PHY/2007/ext-GLORYS12V1_1dAV_20070302_20070303_grid2D_R20070307.nc -m native"
    echo "../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_01.0/ /nird/projects/nird/NS9481K/MERCATOR_DATA/PHY/2013/ext-GLORYS12V1_1dAV_2013110*_grid2D*.nc -b /nird/projects/nird/NS2993K/MERCATOR_DATA/BIO/201"
    echo " NOTE YOU NEED TO RUN THIS SCRIPT WITHIN THE NEMO EXPERIMENT FOLDER"
    exit 1
fi


# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
    export BASEDIR=$(cd .. && pwd)
else
    echo "Could not find EXPT.src. This script must be run in expt dir"
    exit 1
fi
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
#
N="depth_${R}_${T}"
export CDF_NEMO=$1
expt_path=${BASEDIR}/expt_$X
data_path=$expt_path/data
#
export nest_expt=$1
export ncfiles=$2
export mercator_gridfiles=${MERCATOR_NATIVE_MESH}   #${MERCATOR_OLD_PATH}/GRID/GLO_MFC_001_24_MESH.nc
#
iexpt=01        #$T
iversn=22
yrflag=3
idm=4318   #$(egrep "'idm'"  ${data_path}/archv.  | sed "s/.thflag.*$//" | tr -d "[:blank:]")
jdm=1186     #$(egrep "'jdm'"  ${data_path}/expt_${X}/blkdat.input  | sed "s/.iversn.*$//" | tr -d "[:blank:]")
#
# Reading file name using time field of netcdf file. 
#
function model_datetime() {
fname="$1" python - <<END
import cfunits
import datetime
import numpy
import netCDF4
ncid0=netCDF4.Dataset("$1","r")
time=ncid0.variables["time_counter"][0]
unit=ncid0.variables["time_counter"].units
tmp=cfunits.Units(unit)
refy,refm,refd=(1950,1,1)
tmp2=cfunits.Units("hours since %d-%d-%d 00:00:00"%(refy,refm,refd))
tmp3=int(numpy.round(cfunits.Units.conform(time,tmp,tmp2)))
fnametemplate="archv.%Y_%j"
deltat=datetime.datetime(refy,refm,refd,0,0,0)+datetime.timedelta(hours=tmp3)
oname=deltat.strftime(fnametemplate)+"_00"
print oname
END
}
#
#
#
for source_archv in $@ ; do
   # TODO: following 2 lines are for native grid. It needs simply to be modified to be general.
   fn=$(echo ${source_archv:${#source_archv}-32})
   [[ $source_archv != *"ext-GLORYS12V1_1dAV_"* ]] && continue
   echo $source_archv

   filename=$source_archv 
   #
   # (1) Create archive [ab] files from the MERCATOR netcdf file.
   #
   ########################
   if [ ${grid_type} == "native" ] ; then
      if [[ "${bio_path}" == "" ]] ; then
      echo "${BASEDIR}/bin/nemo2archvz_native.py $mercator_gridfiles $source_archv --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}"
      ${BASEDIR}/bin/nemo2archvz_native.py $mercator_gridfiles $source_archv --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}
      ########################
      #
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      ########################
      ${BASEDIR}/bin/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$filename")
      ########################
   
      else
      ${BASEDIR}/bin/nemo2archvz_native.py $mercator_gridfiles $source_archv --bio_path=${bio_path}  --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}
      ########################
      #
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      ########################
      ${BASEDIR}/bin/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$filename") -b 1  
      ########################

      fi
   else
      ${BASEDIR}/bin/nemo2archvz.py $1 $source_archv --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}      
      ########################
      #
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      ########################
      ${BASEDIR}/bin/archvz2hycom.sh $nest_expt $(model_datetime) "$filename" 
      ########################

   fi
done
