#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
#import modeltools.hycom.io
import modeltools.tools
import abfile
import numpy
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

"""
Example: ../bin/hycom_bathy.py --cutoff 0.5 regional.grid.a
"""


def main(infile,blo,bla,shapiro_passes,resolution=None,cutoff=5.) :

   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   scpx=gfile.read_field("scpx")
   scpy=gfile.read_field("scpy")
   width=numpy.median(scpx)/1000.0
   logger.info("Grid median resolution:%8.2f km "%(width))
   # give the bathy directory here
   bathyDir="/cluster/projects/nn2993k/ModelInput/bathymetry/GEBCO_2014/"
   # GEBCO only - TODO: move logic to gebco set
   if resolution is None  :
      if width  > 20 :
         dfile=bathyDir+"GEBCO_2014_2D_median20km.nc"
      elif width > 8 :
         dfile=bathyDir+"GEBCO_2014_2D_median8km.nc"
      elif width > 4 :
         dfile=bathyDir+"GEBCO_2014_2D_median4km.nc"
      else :
         dfile=bathyDir+"GEBCO_2014_2D.nc"
      logger.info ("Source resolution not set - choosing datafile %s"%dfile)
   else :
      dfile=bathyDir+"GEBCO_2014_2D_median%dkm.nc" % resolution
      logger.info ("Source resolution set to %d - trying to use datafile %s"%dfile)
   gebco = modeltools.forcing.bathy.GEBCO2014(filename=dfile)
   w2=gebco.regrid(plon,plat,width=width)
   #print w2.min(),w2.max()

   # Run shapiro filter on interpolated data to remove 2 DeltaX noise
   w3=numpy.copy(w2)
   for i in range(shapiro_passes):
      logger.info("Shapiro filter ... pass %d"%(i+1))
      w3=modeltools.tools.shapiro_filter(w3,threshold=cutoff)
   #print w3.min(),w3.max()

   # Modify basin 
   w4=numpy.copy(-w3)
   it=1
   while it==1 or numpy.count_nonzero(w4-w4old) > 0 :
      w4old = numpy.copy(w4)
      logger.info("Basin modifications ... pass %d"%(it))
      w4=modeltools.tools.remove_isolated_basins(plon,plat,w4,blo,bla,threshold=cutoff)
      w4=modeltools.tools.remove_islets(w4,threshold=cutoff)
      w4=modeltools.tools.remove_one_neighbour_cells(w4,threshold=cutoff)
      logger.info("Modified %d points "%numpy.count_nonzero(w4-w4old) )
      it+=1
   w5=numpy.copy(w4)
   #print w5.min(),w5.max()

   # Mask data where depth below threshold
   w5=numpy.ma.masked_where(w5<=cutoff,w5)

   # Create netcdf file with all  stages for analysis
   logger.info("Writing bathymetry to file hycom_bathymetry.nc")
   ncid = netCDF4.Dataset("hycom_bathymetry.nc","w")
   ncid.createDimension("idm",w5.shape[1])
   ncid.createDimension("jdm",w5.shape[0])
   ncid.createVariable("lon","f8",("jdm","idm"))
   ncid.createVariable("lat","f8",("jdm","idm"))
   ncid.createVariable("h_1","f8",("jdm","idm"))
   ncid.createVariable("h_2","f8",("jdm","idm"))
   ncid.createVariable("h_3","f8",("jdm","idm"))
   ncid.variables["lon"][:]=plon
   ncid.variables["lat"][:]=plat
   ncid.variables["h_1"][:]=w2
   ncid.variables["h_2"][:]=w3
   ncid.variables["h_3"][:]=w5
   ncid.close()

   # Print to HYCOM and CICE bathymetry files
   abfile.write_bathymetry("bathy",1,w5,cutoff)
   
   # Show some grid statistics 
   sx = (w2[1:,:-1]-w2[:-1,:-1])/scpy[1:,:-1]
   sy = (w2[:-1,1:]-w2[:-1,:-1])/scpx[:-1,1:]
   grad = sx + sy
   slopefac=numpy.sqrt(sx**2+sy**2)
   logger.info("Maximum slope factor after interpolation =%.4f"%slopefac.max())

   sx = (w5[1:,:-1]-w5[:-1,:-1])/scpy[1:,:-1]
   sy = (w5[:-1,1:]-w5[:-1,:-1])/scpx[:-1,1:]
   grad = sx + sy
   slopefac=numpy.sqrt(sx**2+sy**2)
   logger.info("Maximum slope factor after smoothing    =%.4f"%slopefac.max())

   
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolor(slopefac)
   figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=w5.min(), vmax=w5.max()))
   #ax.contour(w5)#,[-10.,-100.,-500.,-1000.])
   #ax.set_title("Slope fac in color, depth contours in black")
   logger.info("Slope factor in slopefac.png")
   figure.canvas.print_figure("slopefac.png")



if __name__ == "__main__" :

   class PointParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values[0].split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       tmp1= getattr(args, self.dest)
       tmp1.append(tmp)
       setattr(args, self.dest, tmp1)

   parser = argparse.ArgumentParser(description='Prepare HYCOM bathy files')
   parser.add_argument('--basin_point'   , nargs="*", action=PointParseAction,default=[])
   parser.add_argument('--shapiro_passes', type=int, default=1, help="Number of shapiro passes to apply")
   parser.add_argument('--resolution'    , type=int, default=None)
   parser.add_argument('--cutoff'        , type=float, default=5.0)
   parser.add_argument('regional_grid_file', type=str)
   args = parser.parse_args()

   if args.basin_point :
      blo = [elem[0] for elem in args.basin_point]
      bla = [elem[1] for elem in args.basin_point]
   else :
      blo=[]
      bla=[]


   m = re.match("(.*)(\.[ab]{1})$",args.regional_grid_file) 
   if m :
      regfile=m.group(1)
   else :
      regfile=args.regional_grid_file

   main(regfile,blo,bla,args.shapiro_passes,args.resolution,args.cutoff)

