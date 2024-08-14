#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
import modeltools.tools
import abfile.abfile as abf
import modeltools.cice.io
import numpy as np
import netCDF4
import logging
import re
from scipy.interpolate import griddata
import cartopy.crs as ccrs

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
Example (from the topo-folder of the inner grid - TOPAZ2, merging with NEMO-grid):
../expt_03.1/bin/Grid_Bathy/hycom_bathy_merge_grids.py ../../NMOa0.08/topo/depth_NMOa0.08_01.a ../../NMOa0.08/topo/regional.grid.a depth_TP2a0.10_01.a
To see all options:
../expt_03.1/bin/Grid_Bathy/hycom_bathy_merge_grids.py --help
"""


def main(infile_coarse,gridfile_coarse,infile_fine,
      ncells_linear=20,
      ncells_exact=3,
      check_consistency=False,
      bathy_threshold=0.) :

   #bathy_threshold=0. # TODO
   logger.info("Bathy threshold is %12.4f"%bathy_threshold)

   # Read plon,plat
   gfile=abf.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plon[np.where(plon<=-180)]=plon[np.where(plon<=-180)]+360.0
   plat=gfile.read_field("plat")
   gfile.close()

   # Read input bathymetry - fine version
   m=re.match( "^(.*)(\.[ab])", infile_fine)
   if m : infile_fine=m.group(1)
   bfile=abf.ABFileBathy(infile_fine,"r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   fine_depth_m=bfile.read_field("depth")
   fine_depth_m=np.ma.masked_where(fine_depth_m<=bathy_threshold,fine_depth_m)
   fine_depth=np.ma.filled(fine_depth_m,bathy_threshold)
   bfile.close()

   # Read plon and plat from coarse grid
   gfile=abf.ABFileGrid(gridfile_coarse,"r")
   c_plon=np.ndarray.flatten(gfile.read_field("plon"))
   c_plon[np.where(c_plon>=180)]=c_plon[np.where(c_plon>=180)]-360.0
   c_plat=np.ndarray.flatten(gfile.read_field("plat"))
   gfile.close()

   # Read input bathymetry - coarse version
   m=re.match( "^(.*)(\.[ab])", infile_coarse)
   if m : infile_coarse=m.group(1)
   bfile=abf.ABFileBathy(infile_coarse,"r",idm=gfile.idm,jdm=gfile.jdm)
   coarse_depth_m=np.ndarray.flatten(bfile.read_field("depth"))
   coarse_depth_m[coarse_depth_m>1.0e30]=0.0
   bfile.close()

   proj=ccrs.NorthPolarStereo(central_longitude=-40.0)
   pxy = proj.transform_points(ccrs.PlateCarree(),plon, plat)
   px=pxy[:,:,0]
   py=pxy[:,:,1]

   c_pxy = proj.transform_points(ccrs.PlateCarree(),c_plon, c_plat)
   c_px=c_pxy[:,0]
   c_py=c_pxy[:,1]

   # interpolate coarse bathymetry to the finer grid
   coarse_depth_int=np.ma.zeros(fine_depth.shape)
   coarse_depth_int=griddata((c_px,c_py),coarse_depth_m, (px,py), method='linear')
   del coarse_depth_m
   coarse_depth_m=np.ma.array(coarse_depth_int)
   coarse_depth_m.mask=False
   coarse_depth_m=np.ma.masked_where(coarse_depth_m<=bathy_threshold,coarse_depth_m)
   coarse_depth=np.ma.filled(coarse_depth_m,bathy_threshold)

   # create relaxation mask (rmu)
   tmp=np.linspace(0.,1.,ncells_linear)
   tmp=np.concatenate((np.zeros((ncells_exact,)),tmp)) # ie: hree first cells will match outer bathymetry
   ncells=ncells_linear+ncells_exact
   rmu=np.ones(fine_depth.shape)
   rmu[:,0:ncells] = np.minimum(tmp,rmu[:,0:ncells])
   rmu[0:ncells,:] = np.minimum(tmp,rmu[0:ncells,:].transpose()).transpose()
   rmu[:,-ncells:] = np.minimum(tmp[::-1],rmu[:,-ncells:])
   rmu[-ncells:,:] = np.minimum(tmp[::-1],rmu[-ncells:,:].transpose()).transpose()

   ## Only allow points where both models are defined in the boundaruy
   rmumask=fine_depth_m.mask
   rmumask[:,0:ncells] = np.logical_or(rmumask[:,0:ncells],coarse_depth_m.mask[:,0:ncells])
   rmumask[0:ncells,:] = np.logical_or(rmumask[0:ncells,:],coarse_depth_m.mask[0:ncells,:])
   rmumask[:,-ncells:] = np.logical_or(rmumask[:,-ncells:],coarse_depth_m.mask[:,-ncells:])
   rmumask[-ncells:,:] = np.logical_or(rmumask[-ncells:,:],coarse_depth_m.mask[-ncells:,:])

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
#   print(newbathy.min(),newbathy.max())



   # Make call to consistency routine
   if check_consistency :
      logger.info("Passing merged bathymetry to consistency check ")
      import hycom_bathy_consistency # Normally in same dir as this python routine, so ok
      newbathy=hycom_bathy_consistency.main("",[],[],
            remove_isolated_basins=True,
            remove_one_neighbour_cells=True,
            remove_islets=True,
            remove_inconsistent_nesting_zone=True,
            inbathy=np.ma.masked_where(newbathy<=bathy_threshold,newbathy),
            write_to_file=False)


   # Mask data where depth below threshold
   newbathy_m=np.ma.masked_where(newbathy<=bathy_threshold,newbathy)

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
   I,J=np.where(np.abs(modmask)>.1)
   ax.scatter(J,I,20,"r")
   figure.canvas.print_figure("newbathy.png")



   # Print to HYCOM
   abf.write_bathymetry("MERGED",0,newbathy,bathy_threshold)




if __name__ == "__main__" :


   parser = argparse.ArgumentParser(description='Ensure consistenct of HYCOM bathy files')
   parser.add_argument('--check-consistency', action="store_true",default=False,
         help="Will pass bathymetry to consistency check before finishing")
   parser.add_argument('--ncells-linear', type=int,default=20,
         help="Number of cells in transition zone when going from coarse to fine grid")
   parser.add_argument('--ncells-exact', type=int,default=3,  
         help="Number of grid cells in transition zone having same value  as coarse grid")
   parser.add_argument('--bathy-threshold', type=float,default=0.,
         help="depths shallower than this are marked as land. input >0 !")
   parser.add_argument('infile_coarse', type=str,
         help="bathymetry values  use near edge of domain (normally from a coarse model ...)")
   parser.add_argument('gridfile_coarse', type=str,
         help="grid information near edge of domain (normally from a coarse model ...)")
   parser.add_argument('infile_fine', type=str,
         help="bathymetry values  use in the interior  of domain (normally from a higher resolution model ...)")
   args = parser.parse_args()

   main(args.infile_coarse,args.gridfile_coarse,args.infile_fine,
         ncells_linear=args.ncells_linear,
         ncells_exact=args.ncells_exact,
         check_consistency=args.check_consistency,
         bathy_threshold = args.bathy_threshold,
         )
