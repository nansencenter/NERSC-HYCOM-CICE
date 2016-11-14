#!/bin/bash
# Wrapper file to call ensstat_field repeatedly and produce a 3D
# description of "field" correlations. NB: only outputs the 3D

usage="
This script is a wrapper around ensstat_point to create 3D version
of field correlations between 3D variable field1 and variable field2
at layer index k2, the latter at point xpoint,ypoint.
It calls ensstat_point several times to produce the 3D field. See 
ensstat_point for more info on that routine. 

The routine requires the nco tools to be installed, which does some 
renaming and the final concatenation iof the individual files produced
by ensstat_point

usage :
   $(basename $0) field1 field2 k2 xpoint ypoint hycomfiles
Example:
   $(basename $0) saln saln 1 400 400 hycomfiles
Will create a 3D netcdf file with correlation of salinity in layer
1 (at point 400, 400) with salinity in all layers (at all locations)

NB: Only the 3D variables are present in the final netcdf file
"

if [ $# -lt 6 ] ; then
   echo -e "$usage"
   echo
   echo "Error - too few input arguments"
   exit 1
fi


if ! which ncecat &> /dev/null; then
   echo "Can not find ncecat - make sure nco tools are in your path"
   exit 1
fi
if ! which ncrename &> /dev/null; then
   echo "Can not find ncrename - make sure nco tools are in your path"
   exit 1
fi

# Assumes script is in same location as ensstat prog
ENSSTAT=$(dirname $0)/ensstat_point

field1=$1
field2=$2
kind2=$3
ckind2=$(printf "%2.2d" $kind2)
xind=$4
yind=$5
shift 5

if [ ! -f blkdat.input ]  ; then
   echo "blkdat.input is missing - I quit !"
   exit 1
fi
kdm=$(cat blkdat.input  | grep kdm | sed "s/^[ ]*//" | sed "s/[ \t].*//")

flist=""
#for k in $(seq 1 2) ; do
for k in $(seq 1 $kdm) ; do

   # Rune ensstat_field
   ck=$(printf "%2.2d" $k)
   $ENSSTAT $field1 $k $field2 $kind2 $xind $yind  $@

   # move to temporary file
   mv ensstat_point.nc ensstat_point$ck.nc 

   # Rename variable names sot that ncecat can concatenate them
   # this means changing variable name for first input field by
   # removing layer index of the 3D variable
   #NB: if you change variable names in ensstat_field you will have to change this
   ncrename -v ave_a${field1}$ck,ave_a${field1} ensstat_point$ck.nc  
   ncrename -v var_a${field1}$ck,var_a${field1} ensstat_point$ck.nc 
   ncrename -v cov_a${field1}${ck}_b${field2}$ckind2,cov_a${field1}_b${field2}$ckind2 \
      ensstat_point$ck.nc 
   ncrename -v corr_a${field1}${ck}_b${field2}$ckind2,corr_a${field1}_b${field2}$ckind2 \
   ensstat_point$ck.nc 
   flist="$flist ensstat_point$ck.nc"
done
ncecat -O -u layer -v corr_a${field1}_b${field2}$ckind2,cov_a${field1}_b${field2}$ckind2,\
ave_a${field1},var_a${field1} $flist ensstat_point3D.nc

rm $flist 
echo 
echo "3D fields in ensstat_point3D.nc"

