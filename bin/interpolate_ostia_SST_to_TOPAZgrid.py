from netCDF4 import Dataset
import time as ttimm
import numpy as np
import matplotlib
import abfile
#import mod_reading as mr
from pprint import pprint
import argparse
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import logging
#import cmocean 
import pyproj as pyproj
import matplotlib.patheffects as PathEffects
import os.path
#import scipy.interpolate
from scipy.interpolate import griddata as grd
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
Script to intrepolate OSTAIA global SST into TOPAZ grid 
#usage:
python ./interpolate_ostia_SST_to_TOPAZgrid.py --filename ./global_SST_OSTIA_20161007.nc

"""

def interp2points(fobj,varname,target_lonlats,mapping=None,latlon=True,**kwargs):
    """
    interp2points(fobj,varname,target_lonlats,time_index=0,mapping=None,**kwargs):
    * fobj is a file object eg from netcdf
    * varname is a string (name of variable in file object)
    * latlon (bool)
      if latlon:
         target_lonlats = [target_lon,target_lat], target_lon/lat are numpy arrays
      else:
         target_lonlats = [target_x,target_y], target_x/y are numpy arrays
         need to provide mapping in this case to map lon,lat to x,y
    * time_index (integer) - for multi-record files
    * mapping is a pyproj.Proj object to project form lat/lon to x,y
        (stereographic projection)
    """

    #lon,lat  = fobj.get_lonlat()
    glon = ncfile1.variables['lon'][:] 
    glat = ncfile1.variables['lat'][:] 
    Z0    = ncfile1.variables[varname][:]
    Z=np.transpose(np.squeeze(Z0))
    #mask1  = np.logical_not(np.isfinite(Z))
    #mask1=np.isnan(Z)
    #Z = np.ma.array(Z, mask=mask1)
    print 'Z.shape=', Z.shape
    print 'Z0.shape=', Z0.shape
    print 'lon.shape=', glon.shape
    print 'lon.ndim=', glon.ndim
    print 'lat.shape=',glat.shape
    lon,lat  = np.meshgrid(glon,glat,indexing='ij')
    #lon,lat  = np.meshgrid(glon,glat,indexing='ij')
    #print 'plon=', plon.shape
    #print 'plat=', plat.shape
    # do interpolation in stereographic projection
    if mapping is None:
        if not latlon:
            raise ValueError('need to provide mapping if latlon=False')
        souths    = lat[lat<0]
        lat_0     = 90.
        lon_0     = -45.
        if len(souths)==len(lat):
            # only use south pole projection if all points in southern hemisphere
            lat_0     = -90.
        srs = ('+proj=stere'
                  +' +lon_0='+str(lon_0)
                  +' +lat_0='+str(lat_0)
                  +' +lat_ts='+str(lat_0)
                  +' +ellps=WGS84')

        mapping  = pyproj.Proj(srs)

    #source
    X,Y = mapping(lon,lat)
    #Z    = fobj.get_var(varname,time_index=time_index).values #numpy masked array
    #Z     = ncfile1.variables[varname][:]
    #target
    if latlon:
         lons,lats = target_lonlats
         x,y         = mapping(lons,lats)
    else:
         x,y = target_lonlats

    print 'X.shape=', X.shape
    print 'Y.shape=', Y.shape
    print 'Z.shape=', Z.shape
    print 'x.shape=', x.shape
    print 'y.shape=', y.shape
    # do interpolation
    fvals = reproj_mod2obs(X,Y,Z,x,y,**kwargs) #numpy masked array

    return fvals

#-------------------------------------------------------------------------
# Function that reprojects model into observational grid
def reproj_mod2obs(X1,Y1,Z1,X2,Y2,method='nearest',mask=None):
    # input coords from X1,Y1; Z1 is array to interp; X2,Y2 are output matrices

    # reduce size of source grid
    bbox = X2.min(), X2.max(), Y2.min(), Y2.max()
    X1_reduced, Y1_reduced, Z1_reduced = reduce_grid(X1, Y1, Z1, bbox)
    X1d = X1_reduced.flatten()
    Y1d = Y1_reduced.flatten()
    Z1d = Z1_reduced.flatten()

    # Interpolation
    # - can be done with other methods ('nearest','linear','cubic'<--doesn't work for our data)
    C  = np.array([X1d, Y1d]).T
    Z2 = grd(C, Z1d, (X2, Y2), method=method)

    # check output and add mask
    mask2 = np.isnan(Z2)
    if mask is not None:
        # apply union of mask and model nans
       mask2 = np.logical_or(mask1, mask)

    Z2 = np.ma.array(Z2, mask=mask2)
    return(Z2)

#-------------------------------------------------------------------------
def reduce_grid(X, Y, Z, bbox):
    """
    Reduce the grid so don't have to interp to as many points for example
    Parameters:
    -----------
    X : numpy.ndarray
        array with source x coordinates
    Y : numpy.ndarray
        array with source y coordinates
    Z : numpy.ma.core.MaskedArray
        array with data to be interpolated
    bbox : list(float)
        xmin, xmax, ymin, ymax
    Returns:
    --------
    X_reduced : numpy.ndarray
    Y_reduced : numpy.ndarray
    Z_reduced : numpy.ndarray
        masked values filled with numpy.nan
    """
    print('reduce grid bbox: ', bbox)
    xmin, xmax, ymin, ymax = bbox
    id = np.where( (X>=xmin) & (X<=xmax) &
                 (Y>=ymin)  & (Y<=ymax) )
    imin = id[0].min()
    jmin = id[1].min()
    imax = id[0].max() +1
    jmax = id[1].max() +1
    #Z_ = Z.filled(np.nan).astype(float)
    if hasattr(Z, 'mask'):
	Z_ = Z.filled(np.nan).astype(float)
    else:
	Z_ = Z
    return(
            X[imin: imax, jmin: jmax],
            Y[imin: imax, jmin: jmax],
            Z_[imin: imax, jmin: jmax],
            )
#-------------------------------------------------------------------------

if __name__ == "__main__" :
   class ClimParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       setattr(args, self.dest, tmp)
   class WindowParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [int(elem) for elem in tmp[0:4]]
       setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
   #parser.add_argument('--clim',       action=ClimParseAction,default=None)
   #parser.add_argument('--cmap',       type=str,default="jet",help="matplotlib colormap to use")
   #parser.add_argument('--window',     action=WindowParseAction, help='firsti,firstj,lasti,lastj', default=None)
   #parser.add_argument('--idm',     type=int, help='Grid dimension 1st index []')
   #parser.add_argument('--jdm',     type=int, help='Grid dimension 2nd index []')
   #parser.add_argument('--fieldname',  type=str)
   #parser.add_argument('--fieldlevel', type=int)
   parser.add_argument('--filename', help="",nargs='+')
   #parser.add_argument('records',  nargs="+",    type=int)
   
   args = parser.parse_args()
   print args.filename
   #Example
   #python ./interpolate_sstncof2TP5.py --filename ../ncof_sst/ncof_sst_20*.nc

   ### fh = Dataset('mld_dr003_l3.nc', mode='r')
   ### fh_lons = fh.variables['lon'][:]
   ### fh_lats = fh.variables['lat'][:]
   ### 
   ### fh_time = fh.variables['time'][:]
   ### fh_mld_clim = fh.variables['mld_dr003_rmoutliers_smth_okrg'][:]
   ### fh.close()
   
   gfile = abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()
   dpfile = abfile.ABFileBathy("regional.depth","r",idm=gfile.idm,jdm=gfile.jdm)
   depth=dpfile.read_field("depth")
   dpfile.close()
   target_lon=plon
   target_lat=plat
   target_lonlats = [target_lon,target_lat]

   # compute for TPZ files
   counter=0
   file_count=0
   if args.filename:
      for ncfile0 in args.filename :
            #nci = mr.nc_getinfo(ncfile0)
            ncfile1 = Dataset(ncfile0,'r',format="NETCDF4")
            logger.info("ostia sst data: Now processing  %s"%ncfile0)
            # interpolate ncof_sst in tp5 grid
###             ncfile1 = Dataset(ncfile0,'r',format="NETCDF4")
###             glon = ncfile1.variables['lon'][:] 
###             glat = ncfile1.variables['lat'][:] 
###             gl_mask1 = ncfile1.variables['mask'][:]
###             ncfile1.close()
###             gl_mask = gl_mask1[0,:,:]
###             glonn,glatt=np.meshgrid(glon,glat)
###             print 'glonn=', glonn.shape
###             print 'glatt=',glatt.shape
###             print 'mask=',gl_mask.shape
###             print 'plon=', plon.shape
###             print 'plat=', plat.shape
           #start_time = ttimm.time()
            #print 'start time=', start_time
###         #    newMask=scipy.interpolate.griddata( (glonn.flatten(),glatt.flatten()),gl_mask.flatten(),(plon,plat),'nearest')
###         #    newMask=np.ma.masked_where(np.ma.getmask(depth),newMask)
            #print 'elapsed time=',start_time - ttimm.time()
            #newMask=scipy.interpolate.griddata( (glon.flatten(),glat.flatten()),gl_mask.flatten(),(plon,plat),'nearest')
            start_time = ttimm.time()
            print 'start time=', start_time
            newMask = interp2points(ncfile1,'mask',target_lonlats,mapping=None)
            sst_anlys = interp2points(ncfile1,'analysed_sst',target_lonlats,mapping=None)
            print 'elapsed time=',start_time - ttimm.time()
            # Create scalar field for vectors
            fld=sst_anlys - 273.1
            fld=np.ma.masked_where(np.ma.getmask(depth),fld)
            fld=np.ma.masked_invalid(fld)
            #fld=np.ma.masked_where(fld<-1.8,fld)
            print 'mn,mx  data=',fld.min(),fld.max()
            print 'fldin_nansum=', np.nansum(fld)
            #dt_cl[counter]=numpy.nanmean(fld)
            counter=counter+1
            file_count=file_count+1
            # write to nc file
            suff=os.path.basename(ncfile0)
            ncfilename='TOPAZ_'+suff #.replace('.nc','')
            rootgrp = Dataset(ncfilename, "w", format="NETCDF4")
            logger.info("output to ncfile in  %s"%ncfilename)
            #dimension
            lat = rootgrp.createDimension("lat", gfile.jdm)
            lon = rootgrp.createDimension("lon", gfile.idm)
            time = rootgrp.createDimension("time", None)
            #variable
            times = rootgrp.createVariable("time","f8",("time",))
            ostia_sst = rootgrp.createVariable("analysed_sst","f4",("time","lat","lon",))
            sst_mask = rootgrp.createVariable("ostia_mask","f4",("time","lat","lon",))
            print 'ostia_sst.shape=',ostia_sst.shape
            print 'sst_mask.shape=',sst_mask.shape
            times[:] = 1
            ostia_sst[0,:,:]=fld[:,:]
            sst_mask[0,:,:]=np.array(newMask[:,:])

            rootgrp.close()
            ### End file_intloop
            print 'Computing the avearge of file_counter= ', file_count, 'counter=',counter

