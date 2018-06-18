#!/bin/bash

year=$1
[ $# -eq 0 ] && { echo "No year specified" ; exit 1 ; }

FTPHOST="ftp.cdc.noaa.gov"
FTPPATH="/Datasets/ncep.reanalysis/surface_gauss/"
FTPPATH2="/Datasets/ncep.reanalysis/other_gauss/"


# Translate from field to filenames
filetair="air.2m.gauss.${year}.nc"
fileshum="shum.2m.gauss.${year}.nc"
fileicec="icec.sfc.gauss.${year}.nc"
fileprec="prate.sfc.gauss.${year}.nc"
fileuwnd="uwnd.10m.gauss.${year}.nc"
filevwnd="vwnd.10m.gauss.${year}.nc"
fileuflx="uflx.sfc.gauss.${year}.nc"
filevflx="vflx.sfc.gauss.${year}.nc"
filepres="pres.sfc.gauss.${year}.nc"
filessrd="dswrf.sfc.gauss.${year}.nc"
filelmsk="land.sfc.gauss.nc"
filehgt='hgt.sfc.gauss.nc'

fileccov="tcdc.eatm.gauss.${year}.nc"


# These are in surface_gauss directory
for i in $filetair $fileshum $fileicec $fileprec $fileuwnd $filessrd \
   $filevwnd $filepres $fileuflx $filevflx $filelmsk $filehgt ; do
   echo "ncftpget ftp://$FTPHOST/$FTPPATH/$i"
   ncftpget ftp://$FTPHOST/$FTPPATH/$i
done

# These are in other_gauss directory
echo "ncftpget ftp://$FTPHOST/$FTPPATH2/$fileccov"
ncftpget ftp://$FTPHOST/$FTPPATH2/$fileccov
