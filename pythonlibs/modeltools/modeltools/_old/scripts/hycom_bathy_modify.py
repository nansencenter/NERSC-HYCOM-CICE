#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
from mpl_toolkits.basemap import Basemap
import netCDF4
import logging
import sys
import re

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False


def main(intopo) :

   bathy_threshold=0. # TODO

   # Read plon,plat
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetri
   bfile=abfile.ABFileBathy(intopo,"r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   in_depth_m=bfile.read_field("depth")
   bfile.close()

   # Modify basin 
   in_depth=numpy.ma.filled(in_depth_m,bathy_threshold)
   depth=numpy.copy(in_depth)

   # Open stdin
   logger.info("Reading input data")
   for line in sys.stdin.readlines() :
      logger.info(line.strip())
      ifirst,ilast,jfirst,jlast,d = line.split()
      ifirst=int(ifirst)-1   # Fortran indexing -> C
      ilast =int(ilast)      # Fortran indexing -> C
      jfirst=int(jfirst)-1   # Fortran indexing -> C
      jlast =int(jlast)      # Fortran indexing -> C
      d =float(d)
      logger.info("ifirst=%5d, ilast=%5d, jfirst=%5d, jlast=%5d : depth=%6.2f"%(ifirst,ilast,jfirst,jlast,d))
      depth[jfirst:jlast,ifirst:ilast]=d

   # Mask data where depth below thresholddata = sys.stdin.readlines()

   depth_m=numpy.ma.masked_where(depth<=bathy_threshold,depth)

   # Create netcdf file with all  stages for analysis
   logger.info("Writing bathymetry to file bathy_modify.nc")
   ncid = netCDF4.Dataset("bathy_modify.nc","w")
   ncid.createDimension("idm",depth.shape[1])
   ncid.createDimension("jdm",depth.shape[0])
   ncid.createVariable("lon","f8",("jdm","idm"))
   ncid.createVariable("lat","f8",("jdm","idm"))
   ncid.createVariable("old","f8",("jdm","idm"))
   ncid.createVariable("old_masked","f8",("jdm","idm"))
   ncid.createVariable("new","f8",("jdm","idm"))
   ncid.createVariable("new_masked","f8",("jdm","idm"))
   ncid.createVariable("modified","i4",("jdm","idm"))
   ncid.variables["lon"][:]=plon
   ncid.variables["lat"][:]=plat
   ncid.variables["old"][:]=in_depth
   ncid.variables["old_masked"][:]=in_depth_m
   ncid.variables["new"][:]=depth
   ncid.variables["new_masked"][:]=depth_m
   modmask=numpy.abs(in_depth-depth)>.1
   ncid.variables["modified"][:]=modmask.astype("i4")
   ncid.close()


   # Print to HYCOM and CICE bathymetry files
   abfile.write_bathymetry("MODIFIED",0,depth,bathy_threshold)




if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Ensure consistenct of HYCOM bathy files')
   parser.add_argument('intopo', type=str)
   args = parser.parse_args()

   
   m = re.match("(.*)(\.[ab]{1})$",args.intopo) 
   if m :
      intopo=m.group(1)
   else :
      intopo=args.intopo

   main(intopo)
