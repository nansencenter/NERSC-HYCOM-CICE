#! /bin/bash
# Name: noresm2archvz.py
# Purpose: Convert NORESM data to HYCOM archive files in isopycnal coordinates
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
# (3) July 9 2019, accounting for both bio & phy nesting.
# (4) July 11 2019, further imporvment.
# (5) August 2022, Modify to work with output from NORESM.

options=$(getopt -o b:m -- "$@")
[ $? -eq 0 ] || {
    echo "$usage"
    echo "Error: Incorrect options provided"
    exit 1
}
grid_type=native
bio_path=""
eval set -- "$options"
while true; do
    case "$1" in
    -b)
       shift;  
       bio_path=$1
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

if [ $# -lt 2 ] ; then
    echo "This script will set up the final nesting files from MERCATOR 1/12 degree to be used by HYCOM."
    echo "The code contains the following steps:"
    echo "(1) Create archive [ab] files from the MERCATOR netcdf file."
    echo "    The final archive files includes all 2D fields (filled with zero except baraotropic velocities and sea surface height)."
    echo "    and 3D fields of temperaure, salinity, thickness, and two components of velocities."
    echo ""
    echo "Example:"
    echo " ../bin/Nesting_noresm/noresm_to_hycom.sh ../../TZ2a0.10/expt_06.0/  ../../NORESM_Nesting/thetao_Omon_NorESM2-MM_historical_r1i1p1f1_gr_194001_extrap.nc"
    echo " NOTE YOU NEED TO RUN THIS SCRIPT WITHIN THE ESMa1.00 EXPERIMENT FOLDER"
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
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
#
N="depth_${R}_${T}"
export CDF_NORESM=$1
expt_path=${BASEDIR}/expt_$X
data_path=$expt_path/data
#
export nest_expt=$1
export ncfiles=$2
export noresm_gridfile=${NORESM_NATIVE_MESH}   #areacello_Ofx_NorESM2-MM_historical_r1i1p1f1_gn.nc
#
iexpt=01        #$T
iversn=22
yrflag=3
idm=$(blkdat_get blkdat.input idm)
jdm=$(blkdat_get blkdat.input jdm)
#
# Reading file name using time field of netcdf file. 
#
function model_datetime() {
fname="$1" python - <<END
#import netCDF4
from netCDF4 import Dataset, num2date
ncid=Dataset("$1","r")
time=ncid.variables["time"][0]
tunit=ncid.variables["time"].units                                                                                                  
t_cal=ncid.variables["time"].calendar
ncid.close()
date=num2date(time,units = tunit,calendar = t_cal)
fnametemplate="archv.%Y_%j"
oname=date.strftime(fnametemplate)+"_00"
print(oname)
END
}
#
#
for source_archv in $@ ; do
   fn=$(echo ${source_archv})
   [[ $source_archv != *"Omon_NorESM2-MM"*"extrap.nc" ]] && continue
   filename=$source_archv 
   echo $filename 
   echo $BASEDIR 
   #
   # (1) Create archive [ab] files from the MERCATOR netcdf file.
   #
   ########################
   if [[ "${bio_path}" == "" ]] ; then
      ${BASEDIR}/bin/Nesting_noresm/noresm2archvz.py ${noresm_gridfile} $source_archv \
          --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}
      ########################
      #                                       
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      echo  $(model_datetime "$source_archv")                                              
      ########################                                                            
      ${BASEDIR}/bin/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$source_archv")
      ########################
   else
      ${BASEDIR}/bin/Nesting_noresm/noresm2archvz.py ${noresm_gridfile} $source_archv \
          --bio_path=${bio_path} --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}
      ########################                                                            
      #                                                                                     
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #                                                           
      echo  $(model_datetime "$source_archv")
      ########################                                                                                             
      ${BASEDIR}/bin/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$source_archv") -b 1   
      ########################                                                                                 
   fi
done
