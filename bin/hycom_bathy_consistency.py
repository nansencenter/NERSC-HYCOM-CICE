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


def main(infile,blo,bla,remove_isolated_basins=True,remove_one_neighbour_cells=True,remove_islets=True) :

   bathy_threshold=0. # TODO

   # Read plon,plat
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetri
   bfile=abfile.ABFileBathy(infile,"r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   in_depth_m=bfile.read_field("depth")
   print "in_depth_m type, min, max:",type(in_depth_m),in_depth_m.min(),in_depth_m.max()
   bfile.close()



   # Modify basin 
   in_depth=numpy.ma.filled(in_depth_m,bathy_threshold)
   depth=numpy.copy(in_depth)
   print "depth min max",depth.min(),depth.max()
   it=1
   while it==1 or numpy.count_nonzero(numpy.abs(depth-depth_old)) > 0 :
      depth_old = numpy.copy(depth)
      logger.info("Basin modifications ... pass %d"%(it))
      if remove_isolated_basins     : depth=modeltools.tools.remove_isolated_basins(plon,plat,depth,blo,bla,threshold=bathy_threshold)
      if remove_islets              : depth=modeltools.tools.remove_islets(depth,threshold=bathy_threshold)
      if remove_one_neighbour_cells : depth=modeltools.tools.remove_one_neighbour_cells(depth,threshold=bathy_threshold)
      logger.info("Modified %d points "%numpy.count_nonzero(depth-depth_old) )
      it+=1
   w5=numpy.copy(depth)


   w5[:,0]=bathy_threshold
   w5[:,-1]=bathy_threshold
   w5[0,:]=bathy_threshold
   w5[-1,:]=bathy_threshold
   print "w5 type min max",type(w5),w5.min(),w5.max()

   # Mask data where depth below threshold
   w5_m=numpy.ma.masked_where(w5<=bathy_threshold,w5)

   # Create netcdf file with all  stages for analysis
   logger.info("Writing bathymetry to file bathy_consistency.nc")
   ncid = netCDF4.Dataset("bathy_consistency.nc","w")
   ncid.createDimension("idm",w5.shape[1])
   ncid.createDimension("jdm",w5.shape[0])
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
   ncid.variables["new"][:]=w5
   ncid.variables["new_masked"][:]=w5_m
   modmask=numpy.abs(in_depth-depth)>.1
   ncid.variables["modified"][:]=modmask.astype("i4")
   ncid.close()
   
   logger.info("Writing bathymetry plot to file newbathy.png")
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolormesh(w5_m)
   figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=w5_m.min(), vmax=w5_m.max()))
   I,J=numpy.where(numpy.abs(modmask)>.1)
   ax.scatter(J,I,20,"r")
   figure.canvas.print_figure("newbathy.png")


   # Print to HYCOM and CICE bathymetry files
   abfile.write_bathymetry("CONSISTENT",0,w5,bathy_threshold)


if __name__ == "__main__" :
   class PointParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values[0].split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       tmp1= getattr(args, self.dest)
       tmp1.append(tmp)
       setattr(args, self.dest, tmp1)

   parser = argparse.ArgumentParser(description='Ensure consistenct of HYCOM bathy files')
   parser.add_argument('--no_remove_isolated_basins'    , action="store_true", default=False)
   parser.add_argument('--no_remove_islets'             , action="store_true", default=False)
   parser.add_argument('--no_remove_one_neighbour_cells', action="store_true", default=False)
   parser.add_argument('--basin_point', nargs="*", action=PointParseAction,default=[])
   parser.add_argument('infile', type=str)
   args = parser.parse_args()
   
   m = re.match("(.*)(\.[ab]{1})$",args.infile) 
   if m :
      infile=m.group(1)
   else :
      infile=args.infile

   if args.basin_point :
      blo = [elem[0] for elem in args.basin_point]
      bla = [elem[1] for elem in args.basin_point]
   else :
      blo=[]
      bla=[]
   main(infile,blo,bla,
         remove_isolated_basins    =not args.no_remove_isolated_basins,
         remove_islets             =not args.no_remove_islets        ,
         remove_one_neighbour_cells=not args.no_remove_one_neighbour_cells,
         )
