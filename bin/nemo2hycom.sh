#! /bin/bash
# Name: make_nemo_archives.py
# Purpose: Convert MMERCATOR data to HYCOM archive files in isopycnal coordinates
# Author: Mostafa Bakhoday-Paskyabi (Mostafa.Bakhoday@nersc.no)
# Created: November 2017
# Copyright: (c) NERSC Norway 2017
# Licence:
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


options=$(getopt -o:  -- "$@")
maxinc=50
eval set -- "$options"
while true; do
    case "$1" in
    -m)
       shift;
       maxinc=$1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done

if [ $# -lt 3 ] ; then
    echo "This script will set up the final nesting files from MERCATOR 1/12 degree to be used by HYCOM."
    echo "The code contains the following steps:"
    echo "(1) Create archive [ab] files from the MERCATOR netcdf file."
    echo "    The final archive files includes all 2D fields (filled with zero except baraotropic velocities)."
    echo "    and 3D fields of temperaure, salinity, thickness, and two components of velocities."
    echo "(2) Based on generated archive files in (1) the grid and topography files are generated."
    echo "    Note that the grids are non-native and interpolated into a rectilinear mercator grids horizontally."
    echo ""
    echo "Example:"
    echo "   ../bin/make_nemo_archives.sh ../GLO_MFC_001_24_MESH.nc ../../TP5a0.06/expt_01.0/ ../MERCATOR-PHY-24-2011-01-02*.nc"
    exit 1
fi

export nest_expt=$2
export ncfiles=$3


echo $ncfiles
echo $@

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
HYCOM_ALL="/Users/mba098/Documents/Mostafa/Software/models/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL"
#
iexpt=10
iversn=22
yrflag=3
idm=4320   #$(egrep "'idm'"  ${data_path}/archv.  | sed "s/.thflag.*$//" | tr -d "[:blank:]")
jdm=721     #$(egrep "'jdm'"  ${data_path}/expt_${X}/blkdat.input  | sed "s/.iversn.*$//" | tr -d "[:blank:]")
#
#
#
for source_archv in $@ ; do
fn=$(echo ${source_archv:${#source_archv}-32})
#[ "${fn:0:8}" != "MERCATOR"  ] && continue
[[ $source_archv != *"MERCATOR-PHY-24-"* ]] && continue
#[ echo "$source_archv" | grep -q "MERCATOR" ] && continue

year=$(echo ${source_archv} | awk -F '-' '{print $4}')
month=$(echo $source_archv | awk -F '-' '{print $5}')
vars=$(echo $source_archv | awk -F '-' '{print $6}')
day=$(echo $vars | awk -F '.' '{print $1}')
#( $year, $month, $day ) = ( $source_archv =~ MERCATOR-PHY-24-(\d{4})(\d{2})(\d{2})/ );

filename="MERCATOR-PHY-24-"$year"-"$month"-"$day".nc"

#Use the date command and the %j option. %j denotes day of year as a decimal number, includes a leading zero. (001 to 366)
#doy=$(date -d "${year}-${month}-${day} 00:00:00 UTC"+%j)
#start_oday=$(date -d "${year}-${month}-${day} 00:00:00" +%j)
#
# (1) Create archive [ab] files from the MERCATOR netcdf file.
#
chmod a+x ${BASEDIR}/bin/nemo2archvz.py
${BASEDIR}/bin/nemo2archvz.py $1 $source_archv --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag}

#
# (2) Based on generated archive files in (1) the grid and topography files are generated.
#

${BASEDIR}/bin/remap_nemo.sh $nest_expt $(cat archvname.txt)

done
