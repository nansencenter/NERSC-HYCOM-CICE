


#! /bash/bin


#
# --- This code creates regional and bathymetry files required for offline nesting
# --- Author: Mostafa Bakhoday-Paskyabi, Ocean Modeling group, NERSC, Bergen
# --- Mostafa.Bakhoday@nersc.no
# --- 9 Jan 2018.
#


iscan=15
usage="
In order to use this routine, you need to first do create MESH/MASK file for the subdomain,
used to download MERCATOR data. For Agulhas domain it is:
lonmin=-4; lonmax=70.;
latmin=-53;latmax=-3;

You have to run the script from the experiment directory, and the output are copied in
the topo directory.

Usage:
$(basename $0) PATH/MASK_netcdf_file

Example:
$(basename $0) GLO_MFC_001_24_MASK.nc

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


#export CDF_NEMO="../GLO_MFC_001_24_MASK.nc"

export CDF_NEMO=$1
echo $CDF_NEMO

N="depth_${R}_$T"


#idm=4320    #$(egrep "'idm'"  ${data_path}/archv.  | sed "s/.thflag.*$//" | tr -d "[:blank:]")
#jdm=721     #$(egrep "'jdm'"  ${data_path}/expt_${X}/blkdat.input  | sed "s/.iversn.*$//" | tr -d "[:blank:]")



idm=889     #$(egrep "'idm'"  ${data_path}/archv.  | sed "s/.thflag.*$//" | tr -d "[:blank:]")
jdm=601     #$(egrep "'jdm'"  ${data_path}/expt_${X}/blkdat.input  | sed "s/.iversn.*$//" | tr -d "[:blank:]")



rm -rf regional.grid.[ab] fort.*
echo "${idm}    'idm   ' = longitudinal array size" >> regional.grid.b
echo "${jdm}    'jdm   ' = latitudinal  array size" >> regional.grid.b

#
prog=${HYCOM_ALL}/nemo/src/nemo2grid
#
${prog}
mv fort.12 ../topo/${N}.b
mv fort.012a ../topo/${N}.a
mv fort.011a regional.grid.a
rm fort.11


export FOR051=regional.grid.b
export FOR051A=regional.grid.a


prog=${HYCOM_ALL}/topo/src/grid_lonlat_2d

ln -s regional.grid.a fort.51a
ln -s regional.grid.b fort.51


${prog} < regional.grid.b

mv fort.61 ../topo/regional.grid.b
mv fort.061a ../topo/regional.grid.a

rm -rf fort*

prog=${HYCOM_ALL}/topo/src/topo_landfill

cd ../topo/
ln -s ${N}.b fort.51
ln -s ${N}.a fort.51A

${prog}

rm -rf ../topo/fort*
