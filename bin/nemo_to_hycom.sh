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
# (3) July 9 2019, accounting for both bio & phy nesting.
# (4) July 11 2019, further imporvment.
OPTSTRING="b:d:g:i:m:n:h"
echo "Set defaults that assume grid_type=native"
grid_type=native
maxinc=50
bio_file=""
while getopts ${OPTSTRING} opt; do
  case ${opt} in
    b)
      echo "Option -b biofile was triggered, Argument: ${OPTARG}"
      echo "Default empty string (no biofile)"
      bio_file=${OPTARG}
      ;;
    d)
      echo "Option -d Destination experiment triggered, Argument: ${OPTARG}"
      echo "Mandatory" 
      nest_expt=${OPTARG}
      ;;
    g)
      echo "Option -g grid_type was triggered, Argument: ${OPTARG}"
      echo "Default native Should be either native or regular"
      grid_type=${OPTARG}
      ;;
    i)
      echo "Option -i maxinc - search distance for water points - is triggered, Argument: ${OPTARG}"
      echo "Default 50"
      maxinc=${OPTARG}
      ;;

    m)
      echo "Option -m mercator_mesh was triggered, Argument: ${OPTARG}"
      echo "Default includes path and filename. See script"
      mercator_mesh=${OPTARG}
      ;;
    n)
      echo "Option -n netcdf pattern of input netcdf file was triggered, Argument: ${OPTARG}"
      echo "Mandatory path+name of input file. This accepts wildcard in filename"
      export ncfile="${OPTARG}"
      echo ${OPTARG}
      ;;
    h) 
     echo "This script will set up the final nesting files from MERCATOR 1/12 degree to be used by HYCOM."
     echo "The code contains the following steps:"
     echo "(1) Create archive [ab] files from the MERCATOR netcdf file."
     echo "    The final archive files includes all 2D fields (filled with zero except barotropic velocities)."
     echo "    and 3D fields of temperaure, salinity, thickness, and two components of velocities."
     echo "(2) Based on generated archive files in (1) the grid and topography files are generated."
     echo "    Note that the grids are non-native and interpolated into a rectilinear mercator grids horizontally."
     echo ""
     echo "Example:"
     echo "   pathtobin/nemo_to_hycom.sh -d ../../TP5a0.06/expt_01.0/ -n "/nird/projects/NS9481K/MERCATOR_DATA/PHY/2011/MERCATOR-PHY-24-2011-01-02*.nc" -g regular"
     echo "   pathtobin/nemo_to_hycom.sh -d ../../TP5a0.06/expt_01.2/ -n "/nird/projects/NS9481K/MERCATOR_DATA/PHY/2007/ext-GLORYS12V1_1dAV_20070302_20070303_grid2D_R20070307.nc" -g native"
     echo "../bin/nemo_to_hycom.sh -d../../TP5a0.06/expt_01.0/ -n /nird/projects/NS9481K/MERCATOR_DATA/PHY/2013/ext-GLORYS12V1_1dAV_2013110*_grid2D*.nc -b /nird/projects/NS9481K/MERCATOR_DATA/BIO/DAILY/2013/global_analysis_forecast_bio_2013110*.nc"
     echo " NOTE YOU NEED TO RUN THIS SCRIPT WITHIN THE NEMO EXPERIMENT FOLDER"
     echo " The following arguments are valid:"
     echo "-b: bio file including path"
     echo "-d: Mandatory, Destination experiment including path e.g. TP5a0.06/expt_xx.x/"
     echo "-g: grid_type. Either native or regular"
     echo "-h: This message (help)"
     echo "-i: maxinc. Default 50. Distance where the algorithm search for water points"
     echo "-m: mercator_mesh file. Default (native):  /nird/projects/NS9481K/MERCATOR_DATA/GRID_COORD/ext-GL12V1_mesh_zgr.nc"
     echo "    if -g is set to regular default will change to /nird/projects/NS9481K/MERCATOR_DATA/REGULAR_GRID_COORD/GLO_MFC_001_24_MESH.nc"
     echo "-n: Mandatory: Path and Pattern of Mercator input netCDF files."
     exit 1
     ;;
    :)
      echo "Option -${OPTARG} requires an argument."
      exit 1
      ;;
    ?)
      echo "Invalid option: -${OPTARG}."
      exit 1
      ;;
  esac
done
# Check sanity of input
#if [ -z "$nest_expt" ] || [ -z "$user" ]; then
if [ -z "$nest_expt" ]; then 
        echo "Missing Mandatory argument -d - Path to local domain experiment"
	echo " See nemo_to_hycom -h"
        exit 1
fi
if [ -z "$ncfile" ]; then
        echo "Missing Mandatory argument -n - Path and filename of Mercator input file"
        echo " See nemo_to_hycom -h"
        exit 1
fi
if [ "$grid_type" != "native" ] && [ "$grid_type" != "regular" ]; then
	echo "grid_type (-g) can only be native or regular"
	echo " See nemo_to_hycom -h"
	exit
fi
# timevar hardcoded for now. It may be better to change this /Till 01/05/2025
# Reading file name using time field of netcdf file. 
#
if [ ${grid_type} == "native" ] ; then
        timevar="time_counter"
	if [ -z "$mercator_mesh" ]; then
		export mercator_mesh="/nird/projects/NS9481K/MERCATOR_DATA/GRID_COORD/ext-GL12V1_mesh_zgr.nc"
		echo "mercator_mesh has not been set."
		echo "set to default (native) /nird/projects/NS9481K/MERCATOR_DATA/GRID_COORD/ext-GL12V1_mesh_zgr.nc"
        fi
elif [ ${grid_type} == "regular" ] ; then
        timevar="time"
        if [ -z "$mercator_mesh" ]; then
                export mercator_mesh="/nird/projects/NS9481K/MERCATOR_DATA/REGULAR_GRID_COORD/GLO_MFC_001_24_MESH.nc"
                echo "mercator_mesh has not been set."
		echo "set to default (regular) /nird/projects/NS9481K/MERCATOR_DATA/REGULAR_GRID_COORD/GLO_MFC_001_24_MESH.nc"
        fi
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
#export CDF_NEMO=$1
expt_path=${BASEDIR}/expt_$X
data_path=$expt_path/data
#
iexpt=01        #$T
iversn=22
yrflag=3
idm=$(blkdat_get blkdat.input idm)
jdm=$(blkdat_get blkdat.input jdm)
#

function model_datetime() {
fname="$1" python - <<END
import cfunits
import datetime
import numpy
import netCDF4
ncid0=netCDF4.Dataset("$1","r")
time=ncid0.variables["${timevar}"][0]
unit=ncid0.variables["${timevar}"].units
tmp=cfunits.Units(unit)
refy,refm,refd=(1950,1,1)
tmp2=cfunits.Units("hours since %d-%d-%d 00:00:00"%(refy,refm,refd))
tmp3=int(numpy.round(cfunits.Units.conform(time,tmp,tmp2)))
fnametemplate="archv.%Y_%j"
deltat=datetime.datetime(refy,refm,refd,0,0,0)+datetime.timedelta(hours=tmp3)
oname=deltat.strftime(fnametemplate)+"_00"
print(oname)
END
}

for filename in $ncfile ; do
   #
   # (1) Create archive [ab] files from the MERCATOR netcdf file.
   #
   ########################
   echo $filename
   if [ ${grid_type} == "native" ] ; then
      if [[ "${bio_file}" == "" ]] ; then
      ${BINDIR}/nemo2archvz_native.py $mercator_mesh $filename --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag} --makegrid
      ########################
      #
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      ########################
      ${BINDIR}/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$filename")
      ########################
      else
      ${BINDIR}/nemo2archvz_native.py $mercator_mesh $filename --bio_file=${bio_file}  --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}
      ########################
      #
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      ########################
      ${BINDIR}/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$filename") -b 1  
      ########################

      fi
   else
      if [[ "${bio_file}" == "" ]] ; then
      ${BINDIR}/nemo2archvz_regular.py $mercator_mesh $filename --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}
      ########################
      #
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      ########################
      ${BINDIR}/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$filename") -g regular
      ########################
      else

      ${BINDIR}/nemo2archvz_regular.py $mercator_mesh $filename --bio_file=${bio_file}  --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}
      ########################
      #
      # (2) Based on generated archive files in (1) the grid and topography files are generated.
      #
      ########################
      ${BINDIR}/archvz2hycom_biophys.sh $nest_expt $(model_datetime "$filename") -b 1 -g regular
      ########################

   fi
  fi
done
