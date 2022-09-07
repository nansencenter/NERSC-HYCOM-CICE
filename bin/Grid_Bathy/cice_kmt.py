#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
import abfile.abfile as abf
import modeltools.cice.io
import numpy
import netCDF4
import logging
import sys

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
   # Read plon,plat
   gfile=abf.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetry
   bfile=abf.ABFileBathy(intopo,"r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   in_depth_m=bfile.read_field("depth")
   bfile.close()

   # Print to CICE mask files
   kmt=numpy.where(~in_depth_m.mask,1.,0.)
   modeltools.cice.io.write_netcdf_kmt(kmt,"cice_kmt.nc")


if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Ensure consistenct of HYCOM bathy files')
   parser.add_argument('intopo', type=str)
   args = parser.parse_args()
   main(args.intopo)




