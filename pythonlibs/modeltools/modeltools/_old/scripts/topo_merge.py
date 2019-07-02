#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
#import modeltools.hycom.io
import abfile
import modeltools.cice.io
import numpy
from mpl_toolkits.basemap import Basemap
import netCDF4
import logging

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


def main(infile_coarse,infile_fine,ncells_linear=20,ncells_exact=3) :

   bathy_threshold=0. # TODO

   # Read plon,plat
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetry - fine version
   bfile=abfile.ABFileBathy(infile_fine,"r",idm=gfile.idm,jdm=gfile.jdm)
   fine_depth_m=bfile.read_field("depth")
   fine_depth_m=numpy.ma.masked_where(fine_depth_m<=bathy_threshold,fine_depth_m)
   fine_depth=numpy.ma.filled(fine_depth_m,bathy_threshold)
   bfile.close()

   # Read input bathymetry - coarse version
   bfile=abfile.ABFileBathy(infile_coarse,"r",idm=gfile.idm,jdm=gfile.jdm)
   coarse_depth_m=bfile.read_field("depth")
   coarse_depth_m=numpy.ma.masked_where(coarse_depth_m<=bathy_threshold,coarse_depth_m)
   coarse_depth=numpy.ma.filled(coarse_depth_m,bathy_threshold)
   bfile.close()

   # create relaxation mask (rmu)
   tmp=numpy.linspace(0.,1.,ncells_linear)
   tmp=numpy.concatenate((numpy.zeros((ncells_exact,)),tmp)) # ie: hree first cells will match outer bathymetry
   ncells=ncells_linear+ncells_exact
   rmu=numpy.ones(coarse_depth.shape)
   rmu[:,0:ncells] = numpy.minimum(tmp,rmu[:,0:ncells])
   rmu[0:ncells,:] = numpy.minimum(tmp,rmu[0:ncells,:].transpose()).transpose()
   rmu[:,-ncells:] = numpy.minimum(tmp[::-1],rmu[:,-ncells:])
   rmu[-ncells:,:] = numpy.minimum(tmp[::-1],rmu[-ncells:,:].transpose()).transpose()

   ## Only allow points where both models are defined in the boundaruy
   rmumask=fine_depth_m.mask
   rmumask[:,0:ncells] = numpy.logical_or(rmumask[:,0:ncells],coarse_depth_m.mask[:,0:ncells])
   rmumask[0:ncells,:] = numpy.logical_or(rmumask[0:ncells,:],coarse_depth_m.mask[0:ncells,:])
   rmumask[:,-ncells:] = numpy.logical_or(rmumask[:,-ncells:],coarse_depth_m.mask[:,-ncells:])
   rmumask[-ncells:,:] = numpy.logical_or(rmumask[-ncells:,:],coarse_depth_m.mask[-ncells:,:])

   
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolormesh(rmu)
   figure.colorbar(P)#,norm=matplotlib.colors.LogNorm(vmin=mask.min(), vmax=mask.max()))
   figure.canvas.print_figure("tst.png")


   # Modify bathy in mask region
   newbathy = (1.-rmu) * coarse_depth + rmu * fine_depth
   newbathy[rmumask] = bathy_threshold
   newbathy[:,0]=bathy_threshold
   newbathy[:,-1]=bathy_threshold
   newbathy[0,:]=bathy_threshold
   newbathy[-1,:]=bathy_threshold
   #print newbathy.min(),newbathy.max()

   # Mask data where depth below threshold
   newbathy_m=numpy.ma.masked_where(newbathy<=bathy_threshold,newbathy)


   # Create netcdf file with all  stages for analysis
   logger.info("Writing bathymetry to diagnostic file bathy_merged.nc")
   ncid = netCDF4.Dataset("bathy_merged.nc","w")
   ncid.createDimension("idm",newbathy.shape[1])
   ncid.createDimension("jdm",newbathy.shape[0])
   ncid.createVariable("lon","f8",("jdm","idm"))
   ncid.createVariable("lat","f8",("jdm","idm"))
   ncid.createVariable("coarse","f8",("jdm","idm"))
   ncid.createVariable("coarse_masked","f8",("jdm","idm"))
   ncid.createVariable("fine","f8",("jdm","idm"))
   ncid.createVariable("fine_masked","f8",("jdm","idm"))
   ncid.createVariable("final","f8",("jdm","idm"))
   ncid.createVariable("final_masked","f8",("jdm","idm"))
   ncid.createVariable("rmu","f8",("jdm","idm"))
   ncid.createVariable("modified","f8",("jdm","idm"))
   ncid.variables["lon"][:]=plon
   ncid.variables["lat"][:]=plat
   ncid.variables["coarse"][:]=coarse_depth
   ncid.variables["coarse_masked"][:]=coarse_depth_m
   ncid.variables["fine"][:]=fine_depth
   ncid.variables["fine_masked"][:]=fine_depth_m
   ncid.variables["final"][:]=newbathy
   ncid.variables["final_masked"][:]=newbathy_m
   modmask=newbathy-fine_depth
   ncid.variables["modified"][:] = modmask
   ncid.variables["rmu"][:] = rmu
   ncid.close()
   
   logger.info("Writing bathymetry plot to file newbathy.png")
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolormesh(newbathy)
   figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=newbathy.min(), vmax=newbathy.max()))
   I,J=numpy.where(numpy.abs(modmask)>.1)
   ax.scatter(J,I,20,"r")
   figure.canvas.print_figure("newbathy.png")

   # Print to HYCOM
   abfile.write_bathymetry("MERGED",0,newbathy,bathy_threshold)





if __name__ == "__main__" :


   parser = argparse.ArgumentParser(description='Ensure consistenct of HYCOM bathy files')
   parser.add_argument('--ncells_linear', type=int,default=20)
   parser.add_argument('--ncells_exact', type=int,default=3)
   parser.add_argument('infile_coarse', type=str)
   parser.add_argument('infile_fine', type=str)
   args = parser.parse_args()

   main(args.infile_coarse,args.infile_fine,ncells_linear=args.ncells_linear,ncells_exact=args.ncells_exact)
