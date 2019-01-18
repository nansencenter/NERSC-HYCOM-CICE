#!/usr/bin/env python
import modeltools.nemo
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
import abfile
import numpy
import netCDF4
import logging
import re
import pyproj
import os.path
import time
import cfunits
import nemo_mesh_to_hycom
import sys
import shutil
import glob
from matplotlib import pyplot as plt


# Hycom-ified NMO files. Approach:
# 1) Create hycom archv files and topo/region files using this routine
# 2) Interpolate NEMO topo file to target region
# 3) Merge target region topography and new region. Set up experiment to use new topo
# 4) Interpolate archive file to new region / experiment
# 5) Remap archive files vertically
# 6) Add montgomery potential in 1st layer to file
# 
# Mostafa Bakhoday-Paskyabi, December 2018

# Set up logger
_loglevel=logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False



idm=4320
jdm=1188
sj_sel=slice(0,jdm)
si_sel=slice(0,idm)
latmin=38
latmax=90
lonmin=-180
lonmax=180
mask_method=2
spval=2.0**100.
timeavg_method   = 0             # (1) time average of two consecutive netcdf files; and other values no temporal averaging.


def maplev(a):
    # gapfilling method
    jm,im=a.shape
    J,I=numpy.where(~numpy.isnan(a))
    with numpy.errstate(invalid='ignore'):
        av=numpy.nansum(a[J,I])/(len(I)*len(J))
    J,I=numpy.where(numpy.isnan(a))
    a[J,I]=av
    b=a
    lpp=1  # it is better to be set to 100, but for practial reasnon, we keep it very small for now
    i=range(1,im-1)
    j=range(1,jm-1)
    ip1=range(2,im)
    jp1=range(2,jm)
    im1=range(0,im-2)
    jm1=range(0,jm-2)
    
    cc=numpy.zeros(a.shape)
    for k in range(lpp):
        cc[1:-2,1:-2]=b[1:-2,1:-2]+.5/4*( b[1:-2,2:-1]+b[0:-3,1:-2]+b[1:-2,0:-3]+b[2:-1,1:-2]-4.*b[1:-2,1:-2] )
        cc[:,0]=cc[:,1]
        cc[:,im-1]=cc[:,im-2]
        cc[jm-1,:]=cc[jm-2,:]
        b[J,I]=cc[J,I]

    a[J,I]=cc[J,I]
    J,I=numpy.where(numpy.isnan(a))
    a[J,I]=0.
    return a


def u_to_hycom_u(field2d)  :
   return numpy.roll(field2d,1,axis=1)


# nemo V at i,j -> HYCOM U at i,j+1. 
def v_to_hycom_v(field2d,extrapolate="none")  :
   myfield=numpy.copy(field2d)
   # For now the bottom row is just replicated
   myfield[1:,:] = myfield[:-1,:]
   
   return myfield

# NB: field.shape=(jdm,idm)
def f_to_hycom_q(field,extrapolate="none")  :
   myfield=u_to_hycom_u(field)
   myfield=v_to_hycom_v(myfield)
   return myfield


def periodic_i_shift_right(field,istep) :
   # shift field left by istep steps
   field2  = numpy.roll(field,istep,axis=1)
   return field2


def arctic_patch_shift_up(field,jstep) :
   # shift field down
   if jstep <> 1 :
      raise NameError,"Arctic_patch_shift only with jstep=1 for now"
   field2 = numpy.copy(field)
   field2[1:,:] = field2[0:-1,:] # Shift up
   # NB:  row 0 is same as row 1 (extrapolated
   return field2

def arctic_patch_shift_down(field,jstep) :
   # shift field down
   if jstep <> 1 :
      raise NameError,"Arctic_patch_shift only with jstep=1 for now"
   field2 = numpy.copy(field)
   field2[0:-1,:] = field2[1:,:] # Shift down
   tmp=field2[-1,:]              # Top row as top ...
   field2[-1,:] = tmp[::-1]      # .. but reversed direction
   return field2


def depth_u_points(depth,mbathy,gdepw) :
   depthip1  = periodic_i_shift_right(depth ,-1)    # nemo values at cell i+1
   mbathyip1 = periodic_i_shift_right(mbathy,-1)    # nemo values at cell i+1
   depthu  =numpy.minimum(depth , depthip1)
   mbathy_u=numpy.minimum(mbathy,mbathyip1)
   e3u_ps  =numpy.zeros(depthu.shape)
   J,I= numpy.where(mbathy_u>-1)
   e3u_ps[J,I] = depthu[J,I] - gdepw[mbathy_u[J,I]]
   return mbathy_u,e3u_ps,depthu


      
def depth_v_points(depth,mbathy,gdepw) :
   depthjp1  = arctic_patch_shift_down(depth ,1)    # nemo values at cell j+1
   mbathyjp1 = arctic_patch_shift_down(mbathy,1)    # nemo values at cell j+1 
   depthv  =numpy.minimum(depth ,depthjp1)
   mbathy_v=numpy.minimum(mbathy,mbathyjp1)
   e3v_ps  =numpy.zeros(depthv.shape)
   I= numpy.where(mbathy_v>-1)
   e3v_ps[I] = depthv[I] - gdepw[mbathy_v[I]]
   return mbathy_v,e3v_ps,depthv

def sliced(field2d) :
   return field2d[sj_sel,si_sel]



def read_mesh(filemesh):
   meshfile_hgr=filemesh[:-6]+"hgr.nc"
   meshfile_zgr=filemesh[:-6]+"zgr.nc"
   maskfile=filemesh[:-11]+"mask.nc"

   ncidh=netCDF4.Dataset(meshfile_hgr,"r")
   ncidz=netCDF4.Dataset(meshfile_zgr,"r")
   ncidm=netCDF4.Dataset(maskfile,"r")

   depth = sliced(ncidz.variables["hdepw" ][0,:,:]) 
   gdept = ncidz.variables["gdept_0"][0,:]       # Depth of t points
   gdepw  = ncidz.variables["gdepw_0"][0,:]       # Depth of w points
   e3t_ps = sliced(ncidz.variables["e3t_ps"] [0,:,:])     # Partial steps of t cell
   e3w_ps = sliced(ncidz.variables["e3w_ps"] [0,:,:])     # Partial steps of w cell
   mbathy = sliced(ncidz.variables["mbathy"] [0,:,:])     # bathy index
   hdepw  = sliced(ncidz.variables["hdepw"]  [0,:,:])     # Total depth of w points
   plon = sliced(ncidz.variables["nav_lon"][:,:])     # bathy index
   plat  = sliced(ncidz.variables["nav_lat"][:,:])     # Total depth of w points

   return gdept,gdepw,e3t_ps,e3w_ps,mbathy,hdepw,depth,plon,plat



def make_grid(filemesh):

   meshfile_hgr=filemesh[:-6]+"hgr.nc"
   meshfile_zgr=filemesh[:-6]+"zgr.nc"
   maskfile=filemesh[:-11]+"mask.nc"

   ncidh=netCDF4.Dataset(meshfile_hgr,"r")
   ncidz=netCDF4.Dataset(meshfile_zgr,"r")
   ncidm=netCDF4.Dataset(maskfile,"r")

   # Now acquire the data. P-cell data
   plon = sliced(ncidh.variables["glamt"][0,:,:])
   plat = sliced(ncidh.variables["gphit"][0,:,:])
   scpx = sliced(ncidh.variables["e1t"]  [0,:,:])
   scpy = sliced(ncidh.variables["e2t"]  [0,:,:])

   # U-cell data. 
   ulon = sliced(u_to_hycom_u( ncidh.variables["glamu"][0,:,:]))
   ulat = sliced(u_to_hycom_u( ncidh.variables["gphiu"][0,:,:]))
   scux = sliced(u_to_hycom_u( ncidh.variables["e1u"  ][0,:,:]))
   scuy = sliced(u_to_hycom_u( ncidh.variables["e2u"  ][0,:,:]))

   # V-cell data.
   # TODO: Proper extrapolation of data on grid edges (mainly bottom row)
   vlon = sliced(v_to_hycom_v( ncidh.variables["glamv"][0,:,:]))
   vlat = sliced(v_to_hycom_v( ncidh.variables["gphiu"][0,:,:]))
   scvx = sliced(v_to_hycom_v( ncidh.variables["e1v"]  [0,:,:]))
   scvy = sliced(v_to_hycom_v( ncidh.variables["e2v"]  [0,:,:]))

   # Q-cell data
   # TODO: Proper extrapolation of data on grid edges (mainly bottom row)
   qlon = sliced(f_to_hycom_q( ncidh.variables["glamf"][0,:,:] ))
   qlat = sliced(f_to_hycom_q( ncidh.variables["gphif"][0,:,:] ))
   scqx = sliced(f_to_hycom_q( ncidh.variables["e1f"  ][0,:,:] ))
   scqy = sliced(f_to_hycom_q( ncidh.variables["e2f"  ][0,:,:] ))

   # Angle used for rotation
   ulon_rgt = numpy.copy(ulon)
   ulat_rgt = numpy.copy(ulat)
   ulon_lft = numpy.roll(ulon,1,axis=1)
   ulat_lft = numpy.roll(ulat,1,axis=1)
   pang = modeltools.tools.p_azimuth(ulon_lft,ulat_lft,ulon_rgt,ulat_rgt)

   # Aspect ratio
   asp = numpy.where(scpy==0.,99.0,scpx/scpy)

   # Coriolis
   corio = numpy.sin(numpy.radians(qlat)) * 4. * numpy.pi / 86164.0 # Sidereal day

   # Put inside datadict for abfile writing
   ddict={}
   ddict["plon"]=plon
   ddict["plat"]=plat
   ddict["ulon"]=ulon
   ddict["ulat"]=ulat
   ddict["vlon"]=vlon
   ddict["vlat"]=vlat
   ddict["qlon"]=qlon
   ddict["qlat"]=qlat
   #
   ddict["scpx"]=scpx
   ddict["scpy"]=scpy
   ddict["scux"]=scux
   ddict["scuy"]=scuy
   ddict["scvx"]=scvx
   ddict["scvy"]=scvy
   ddict["scqx"]=scqx
   ddict["scqy"]=scqy
   #
   ddict["cori"]=corio
   ddict["pang"] =pang
   ddict["pasp"] =asp
   abfile.write_regional_grid(ddict)

   # Bathymetry
   hdepw  = sliced(ncidz.variables["hdepw"] [0,:,:])
   mbathy = sliced(ncidz.variables["mbathy"][0,:,:])

   ncidt=netCDF4.Dataset(maskfile,"r")
   tmask = sliced(ncidt.variables["tmaskutil"][0,:,:])
   tmp2 = numpy.where( mbathy>0.  ,hdepw,0.)
   abfile.write_bathymetry("bathy",1,tmp2,0.)
 
   shutil.move("depth_bathy_01.a",'../topo/depth_NMOa0.08_01.a')
   shutil.move("depth_bathy_01.b",'../topo/depth_NMOa0.08_01.b')
   shutil.move("regional.grid.a","../topo/regional.grid.a")
   shutil.move("regional.grid.b","../topo/regional.grid.b")

def search_biofile(bio_path,dt):

   logger.info("BIO")
   # filename format MERCATOR-BIO-14-2013-01-05-00
   lst=glob.glob(bio_path+"/MERCATOR-BIO-14-%s*.nc"%str(dt[:-2]))
   df = numpy.zeros(len(lst))*numpy.nan 
   for num,myfile in enumerate(lst):
      tmp = datetime.datetime.strptime(myfile[-16:-6],'%Y-%M-%d')-datetime.datetime.strptime(dt,'%Y-%M-%d')
      df[num]=tmp.days
   val, idx = min((val, idx) for (idx, val) in enumerate(numpy.abs(df)))
   return idx,lst[idx]

def check_inputs(x, y, Z, points, mode, bounds_error):
    """Check inputs for interpolate2d function
    """

    msg = 'Only mode "linear" and "constant" are implemented. I got %s' % mode
    if mode not in ['linear', 'constant']:
        raise RuntimeError(msg)

    try:
        x = numpy.array(x)
    except Exception, e:
        msg = ('Input vector x could not be converted to numpy array: '
               '%s' % str(e))
        raise Exception(msg)

    try:
        y = numpy.array(y)
    except Exception, e:
        msg = ('Input vector y could not be converted to numpy array: '
               '%s' % str(e))
        raise Exception(msg)

    msg = ('Input vector x must be monotoneously increasing. I got '
           'min(x) == %.15f, but x[0] == %.15f' % (min(x), x[0]))
    if not min(x) == x[0]:
        raise RuntimeError(msg)

    msg = ('Input vector y must be monotoneously increasing. '
           'I got min(y) == %.15f, but y[0] == %.15f' % (min(y), y[0]))
    if not min(y) == y[0]:
        raise RuntimeError(msg)

    msg = ('Input vector x must be monotoneously increasing. I got '
           'max(x) == %.15f, but x[-1] == %.15f' % (max(x), x[-1]))
    if not max(x) == x[-1]:
        raise RuntimeError(msg)

    msg = ('Input vector y must be monotoneously increasing. I got '
           'max(y) == %.15f, but y[-1] == %.15f' % (max(y), y[-1]))
    if not max(y) == y[-1]:
        raise RuntimeError(msg)

    try:
        Z = numpy.array(Z)
        m, n = Z.shape
    except Exception, e:
        msg = 'Z must be a 2D numpy array: %s' % str(e)
        raise Exception(msg)

    Nx = len(x)
    Ny = len(y)
    msg = ('Input array Z must have dimensions %i x %i corresponding to the '
           'lengths of the input coordinates x and y. However, '
           'Z has dimensions %i x %i.' % (Nx, Ny, m, n))
    if not (Nx == m and Ny == n):
        raise RuntimeError(msg)

    # Get interpolation points
    points = numpy.array(points)
    xi = points[:, 0]
    eta = points[:, 1]

    if bounds_error:
        msg = ('Interpolation point %f was less than the smallest value in '
               'domain %f and bounds_error was requested.' % (xi[0], x[0]))
        if xi[0] < x[0]:
            raise Exception(msg)

        msg = ('Interpolation point %f was greater than the largest value in '
               'domain %f and bounds_error was requested.' % (xi[-1], x[-1]))
        if xi[-1] > x[-1]:
            raise Exception(msg)

        msg = ('Interpolation point %f was less than the smallest value in '
               'domain %f and bounds_error was requested.' % (eta[0], y[0]))
        if eta[0] < y[0]:
            raise Exception(msg)

        msg = ('Interpolation point %f was greater than the largest value in '
               'domain %f and bounds_error was requested.' % (eta[-1], y[-1]))
        if eta[-1] > y[-1]:
            raise Exception(msg)

    return x, y, Z, xi, eta


def interpolate2d(x, y, Z, points, mode='linear', bounds_error=False):
    """Fundamental 2D interpolation routine
    Input
        x: 1D array of x-coordinates of the mesh on which to interpolate
        y: 1D array of y-coordinates of the mesh on which to interpolate
        Z: 2D array of values for each x, y pair
        points: Nx2 array of coordinates where interpolated values are sought
        mode: Determines the interpolation order. Options are
              'constant' - piecewise constant nearest neighbour interpolation
              'linear' - bilinear interpolation using the four
                         nearest neighbours (default)
        bounds_error: Boolean flag. If True (default) an exception will
                      be raised when interpolated values are requested
                      outside the domain of the input data. If False, nan
                      is returned for those values
    Output
        1D array with same length as points with interpolated values
    Notes
        Input coordinates x and y are assumed to be monotonically increasing,
        but need not be equidistantly spaced.
        Z is assumed to have dimension M x N, where M = len(x) and N = len(y).
        In other words it is assumed that the x values follow the first
        (vertical) axis downwards and y values the second (horizontal) axis
        from left to right.
        If this routine is to be used for interpolation of raster grids where
        data is typically organised with longitudes (x) going from left to
        right and latitudes (y) from left to right then user
        interpolate_raster in this module
    """

    # Input checks
    x, y, Z, xi, eta = check_inputs(x, y, Z, points, mode, bounds_error)

    # Identify elements that are outside interpolation domain or NaN
    outside = (xi < x[0]) + (eta < y[0]) + (xi > x[-1]) + (eta > y[-1])
    outside += numpy.isnan(xi) + numpy.isnan(eta)

    inside = ~outside
    xi = xi[inside]
    eta = eta[inside]

    # Find upper neighbours for each interpolation point
    idx = numpy.searchsorted(x, xi, side='left')
    idy = numpy.searchsorted(y, eta, side='left')

    # Internal check (index == 0 is OK)
    msg = ('Interpolation point outside domain. This should never happen. '
           'Please email Ole.Moller.Nielsen@gmail.com')
    if len(idx) > 0:
        if not max(idx) < len(x):
            raise RuntimeError(msg)
    if len(idy) > 0:
        if not max(idy) < len(y):
            raise RuntimeError(msg)

    # Get the four neighbours for each interpolation point
    x0 = x[idx - 1]
    x1 = x[idx]
    y0 = y[idy - 1]
    y1 = y[idy]

    z00 = Z[idx - 1, idy - 1]
    z01 = Z[idx - 1, idy]
    z10 = Z[idx, idy - 1]
    z11 = Z[idx, idy]

    # Coefficients for weighting between lower and upper bounds
    oldset = numpy.seterr(invalid='ignore')  # Suppress warnings
    alpha = (xi - x0) / (x1 - x0)
    beta = (eta - y0) / (y1 - y0)
    numpy.seterr(**oldset)  # Restore

    if mode == 'linear':
        # Bilinear interpolation formula
        dx = z10 - z00
        dy = z01 - z00
        z = z00 + alpha * dx + beta * dy + alpha * beta * (z11 - dx - dy - z00)
    else:
        # Piecewise constant (as verified in input_check)

        # Set up masks for the quadrants
        left = alpha < 0.5
        right = -left
        lower = beta < 0.5
        upper = -lower

        lower_left = lower * left
        lower_right = lower * right
        upper_left = upper * left

        # Initialise result array with all elements set to upper right
        z = z11

        # Then set the other quadrants
        z[lower_left] = z00[lower_left]
        z[lower_right] = z10[lower_right]
        z[upper_left] = z01[upper_left]

    # Self test
    if len(z) > 0:
        mz = numpy.nanmax(z)
        mZ = numpy.nanmax(Z)
        msg = ('Internal check failed. Max interpolated value %.15f '
               'exceeds max grid value %.15f ' % (mz, mZ))
        if not(numpy.isnan(mz) or numpy.isnan(mZ)):
            if not mz <= mZ:
                raise RuntimeError(msg)

    # Populate result with interpolated values for points inside domain
    # and NaN for values outside
    r = numpy.zeros(len(points))
    r[inside] = z
    r[outside] = numpy.nan

    return r



def main(filemesh,grid2dfiles,first_j=0,mean_file=False,iexpt=10,iversn=22,yrflag=3,makegrid=None,bio_path=None) :

   if mean_file :
      fnametemplate="archm.%Y_%j_%H"
   else :
      fnametemplate="archv.%Y_%j_%H"
   itest=1
   jtest=200
   gdept,gdepw,e3t_ps,e3w_ps,mbathy,hdepw,depth,plon,plat=read_mesh(filemesh)
   if makegrid is not None: 
      logger.info("Making NEMO grid & bathy [ab] files ...")
      make_grid(filemesh)
   print plon.min(),plon.max()
   mbathy = mbathy -1                       # python indexing starts from 0
   nlev   = gdept.size

   mbathy_u,e3u_ps,depthu=depth_u_points(depth,mbathy,gdepw)
   mbathy_v,e3v_ps,depthv=depth_v_points(depth,mbathy,gdepw)
   #
   mbathy_u=sliced(u_to_hycom_u(mbathy_u))
   e3u_ps  =sliced(u_to_hycom_u(e3u_ps  ))
   depthu  =sliced(u_to_hycom_u(depthu  ))
   #
   mbathy_v=sliced(v_to_hycom_v(mbathy_v))
   e3v_ps  =sliced(v_to_hycom_v(e3v_ps  ))
   depthv  =sliced(v_to_hycom_v(depthv  ))

   # Thickness of t layers (NB: 1 less than gdepw dimension)
   dt = gdepw[1:] - gdepw[:-1]

   # Loop over input files. All must be in same directory
   for file2d in grid2dfiles : 

      # See if actually a grid2D file
      dirname=os.path.dirname(file2d)
      m=re.match("(.*_)(grid2D)(_.*\.nc)",os.path.basename(file2d))
      if not m :
         msg="File %s is not a grid2D file, aborting"%file2d
         logger.error(msg)
         raise ValueError,msg

      # Construct remaining files
      filet  =os.path.join(dirname,m.group(1) + "gridT" + m.group(3))
      files  =os.path.join(dirname,m.group(1) + "gridS" + m.group(3))
      fileu  =os.path.join(dirname,m.group(1) + "gridU" + m.group(3))
      filev  =os.path.join(dirname,m.group(1) + "gridV" + m.group(3))
      filew  =os.path.join(dirname,m.group(1) + "gridW" + m.group(3))
      fileice=os.path.join(dirname,m.group(1) + "icemod" + m.group(3))
      logger.info("grid2D file: %s"%file2d)

      # P-points
      logger.info("gridS  file: %s"%files)
      logger.info("gridT  file: %s"%filet)
      ncids=netCDF4.Dataset(files,"r")
      ncidt=netCDF4.Dataset(filet,"r")

      # time from gridT file. 
      time = ncidt.variables["time_counter"][0]
      tunit = ncidt.variables["time_counter"].units
      tmp=cfunits.Units(tunit)
      refy, refm, refd=(1958,1,1)
      tmp2=cfunits.Units("hours since %d-%d-%d 00:00:00"%(refy,refm,refd))            # Units from CF convention
      tmp3=cfunits.Units.conform(time,tmp,tmp2)                                       # Transform to new new unit 
      tmp3=int(numpy.round(tmp3))
      mydt = datetime.datetime(refy,refm,refd,0,0,0) + datetime.timedelta(hours=tmp3) # Then calculate dt. Phew!
      if bio_path:
         jdm,idm=numpy.shape(plon)
         points = numpy.transpose(((plat.flatten(),plon.flatten())))
         delta = mydt.strftime( '%Y-%m-%d')
         # filename format MERCATOR-BIO-14-2013-01-05-00
         idx,biofname=search_biofile(bio_path,delta)
         if idx >7: 
            msg="No available BIO file within a week difference with PHY"
	    logger.error(msg)
            raise ValueError,msg
         logger.info("BIO file %s reading & interpolating to 1/12 deg grid cells ..."%biofname)
         ncidb=netCDF4.Dataset(biofname,"r")
         blon=ncidb.variables["longitude"][:];
         blat=ncidb.variables["latitude"][:]
         no3=ncidb.variables["NO3"][0,:,:,:];
         no3[numpy.abs(no3)>1e+10]=numpy.nan
         po4=ncidb.variables["PO4"][0,:,:,:]
         si=ncidb.variables["Si"][0,:,:,:]
         po4[numpy.abs(po4)>1e+10]=numpy.nan
         si[numpy.abs(si)>1e+10]=numpy.nan
         # TODO: Ineed to improve this part
         nz,ny,nx=no3.shape
         dummy=numpy.zeros((nz,ny,nx+1))
         dummy[:,:,:nx]=no3;dummy[:,:,-1]=no3[:,:,-1]
         no3=dummy
         dummy=numpy.zeros((nz,ny,nx+1))
         dummy[:,:,:nx]=po4;dummy[:,:,-1]=po4[:,:,-1]
         po4=dummy
         dummy=numpy.zeros((nz,ny,nx+1))
         dummy[:,:,:nx]=si;dummy[:,:,-1]=si[:,:,-1]
         si=dummy
         dummy=numpy.zeros((nx+1))
         dummy[:nx]=blon
         blon=dummy
         blon[-1]=-blon[0]
      # Read and calculculate U in hycom U-points. 
      logger.info("gridU, gridV, gridT & gridS  file")
      ncidu=netCDF4.Dataset(fileu,"r")
      u=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
      ncidv=netCDF4.Dataset(filev,"r")
      v=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
      udummy=ncidu.variables["vozocrtx"][:,:,:,:] 
      vdummy=ncidv.variables["vomecrty"][:,:,:,:]
      tdummy=ncidt.variables["votemper"][:,:,:,:]
      tdummy_fill=ncidt.variables["votemper"]._FillValue
      sdummy=ncids.variables["vosaline"][:,:,:,:]
      sdummy_fill=ncids.variables["vosaline"]._FillValue

      for k in range(nlev) : 
         u[k,:,:] = sliced(u_to_hycom_u(udummy[0,k,:,:] ))   # Costly, make more efficient if needed
         v[k,:,:] = sliced(v_to_hycom_v(vdummy[0,k,:,:] ))   # Costly, make more efficient if needed

      u = numpy.where(numpy.abs(u)<1e10,u,0.)
      v = numpy.where(numpy.abs(v)<1e10,v,0.)
      logger.info("Calculate barotropic velocities ...")

      #Calculate barotropic and baroclinic u
      usum=numpy.zeros(u.shape[-2:])
      dsumu=numpy.zeros(u.shape[-2:])
      vsum=numpy.zeros(v.shape[-2:])
      dsumv=numpy.zeros(v.shape[-2:])

      for k in range(u.shape[0]-1) : # Dont include lowest layer
         J,I = numpy.where(mbathy_u>k) 
         usum[J,I] = usum[J,I] + u[k,J,I]*dt[k]
         dsumu[J,I] = dsumu[J,I] + dt[k]
         J,I = numpy.where(mbathy_v>k)
         vsum[J,I] = vsum[J,I] + v[k,J,I]*dt[k]
         dsumv[J,I] = dsumv[J,I] + dt[k]
      J,I = numpy.where(mbathy>=0)
      usum[J,I] = usum[J,I] + u[mbathy_u[J,I],J,I]*e3u_ps[J,I]
      dsumu[J,I] = dsumu[J,I] + e3u_ps[J,I]
      dsumu=numpy.where(abs(dsumu)<1e-2,0.05,dsumu)
      ubaro=numpy.where(dsumu>0.1,usum/dsumu,0.)
      J,I = numpy.where(mbathy_v>=0)
      vsum[J,I] = vsum[J,I] + v[mbathy_v[J,I],J,I]*e3v_ps[J,I]
      dsumv[J,I] = dsumv[J,I] + e3v_ps[J,I]
      dsumv=numpy.where(abs(dsumv)<1e-2,0.05,dsumv)
      vbaro=numpy.where(dsumv>.1,vsum/dsumv,0.)
      fnametemplate="archv.%Y_%j"
      deltat=datetime.datetime(refy,refm,refd,0,0,0)+datetime.timedelta(hours=tmp3)
      oname=deltat.strftime(fnametemplate)+"_00"

      # model day
      refy, refm, refd=(1900,12,31)
      model_day= deltat-datetime.datetime(refy,refm,refd,0,0,0)
      model_day=model_day.days
      logger.info("Model day in HYCOM is %s"%str(model_day))

      # Masks (land:True)
      if mask_method == 1 :
         ip = mbathy   == -1
         iu = mbathy_u == -1
         iv = mbathy_v == -1
      else :
         ip = depth   == 0
         iu = depthu  == 0
         iv = depthv  == 0

      flnm = open('archvname.txt', 'w')
      flnm.write(oname)
      flnm.close()

      # 2D data
      ncid2d=netCDF4.Dataset(file2d,"r")
      ssh          = sliced(ncid2d.variables["sossheig"][0,:,:])
      ssh = numpy.where(ssh==ncid2d.variables["sossheig"]._FillValue,0.,ssh)
      ssh = numpy.where(ssh>1e10,0.,ssh*9.81) # NB: HYCOM srfhgt is in geopotential ...
      montg1=numpy.zeros(ssh.shape)

      # Write to abfile
      outfile = abfile.ABFileArchv("./data/"+oname,"w",iexpt=iexpt,iversn=iversn,yrflag=yrflag,)

      logger.info("Writing 2D variables")
      outfile.write_field(montg1,                ip,"montg1"  ,0,model_day,1,0)
      outfile.write_field(ssh,                   ip,"srfhgt"  ,0,model_day,0,0)
      outfile.write_field(numpy.zeros(ssh.shape),ip,"surflx"  ,0,model_day,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),ip,"salflx"  ,0,model_day,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),ip,"bl_dpth" ,0,model_day,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),ip,"mix_dpth",0,model_day,0,0) # Not used
      outfile.write_field(ubaro                 ,iu,"u_btrop" ,0,model_day,0,0) # u: nemo in cell i is hycom in cell i+1
      outfile.write_field(vbaro                 ,iv,"v_btrop" ,0,model_day,0,0) # v: nemo in cell j is hycom in cell j+1
      for k in numpy.arange(u.shape[0]) :
         if bio_path:
	    no3k=interpolate2d(blat, blon, no3[k,:,:], points)
            no3k = maplev(numpy.reshape(no3k,[jdm,idm]))
            po4k=interpolate2d(blat, blon, po4[k,:,:], points)
            po4k = maplev(numpy.reshape(po4k,[jdm,idm]))
            si_k=interpolate2d(blat, blon, si[k,:,:], points)
            si_k = maplev(numpy.reshape(si_k,[jdm,idm]))
#            fig = plt.figure(figsize=(16,8));plt.subplot(121);plt.pcolormesh(no3k,vmin=0,vmax=20);plt.colorbar();plt.show();
#            plt.subplot(122);plt.pcolormesh(no3[0,:,:],vmin=0,vmax=20);plt.colorbar();plt.show();fig.savefig('intrplolated_no3.png')
#            exit(0)
         if bio_path :
            if k%10==0 : logger.info("Writing 3D variables including BIO, level %d of %d"%(k+1,u.shape[0]))
         else:
            if k%10==0 : logger.info("Writing 3D variables, level %d of %d"%(k+1,u.shape[0]))
         ul = numpy.squeeze(u[k,:,:]) - ubaro # Baroclinic velocity
         vl = numpy.squeeze(v[k,:,:]) - vbaro # Baroclinic velocity

         # Layer thickness
         dtl=numpy.zeros(ul.shape)
         if k < u.shape[0]-1 :
            J,I = numpy.where(mbathy>k)
            dtl[J,I] = dt[k]
            J,I = numpy.where(mbathy==k)
            dtl[J,I] = e3t_ps[J,I]
         else:
            J,I = numpy.where(mbathy==k)
            dtl[J,I] = e3t_ps[J,I]

         tmpfill=sdummy_fill#ncids.variables["vosaline"]._FillValue
         sl = sliced(sdummy[0,k,:,:])
         tmpfill=tdummy_fill#ncidt.variables["votemper"]._FillValue
         tl = sliced(tdummy[0,k,:,:])
         sl = numpy.where(numpy.abs(sl)<1e2,sl,numpy.nan)
         sl = numpy.minimum(numpy.maximum(maplev(sl),25),80.)
         tl = numpy.where(numpy.abs(tl)<=5e2,tl,numpy.nan)
         tl = numpy.minimum(numpy.maximum(maplev(tl),-5.),50.)

         # Fill empty layers with values from above
         if k > 0 :
            K= numpy.where(dtl < 1e-4)

            tl[K] = tl_above[K]


         onem=9806.
         outfile.write_field(ul      ,iu,"u-vel.",0,model_day,k+1,0) # u: nemo in cell i is hycom in cell i+1
         outfile.write_field(vl      ,iv,"v-vel.",0,model_day,k+1,0) # v: nemo in cell j is hycom in cell j+1
         outfile.write_field(dtl*onem,ip,"thknss",0,model_day,k+1,0)
         outfile.write_field(tl      ,ip,"temp"  ,0,model_day,k+1,0)
         outfile.write_field(sl      ,ip,"salin" ,0,model_day,k+1,0)
         if bio_path :
            outfile.write_field(no3k      ,ip,"Eco_no3"  ,0,model_day,k+1,0)
            outfile.write_field(po4k      ,ip,"Eco_pho" ,0,model_day,k+1,0)
            outfile.write_field(si_k      ,ip,"Echo_sil" ,0,model_day,k+1,0)

         tl_above=numpy.copy(tl)
         sl_above=numpy.copy(sl)
         

      # TODO: Process ice data
      ncid2d.close()
      outfile.close()
      ncidt.close()
      ncids.close()
      ncidu.close()
      ncidv.close()
      ncidb.close()
   nemo_mesh = []





if __name__ == "__main__" :

   parser = argparse.ArgumentParser(
         description='This tool will convert NEMO netcdf files to hycom archive files. It will also create grid and topo files for hycom.'
         )
   parser.add_argument('--first-j',   type=int,default=0,help="first j-index to process. Defaults to 0")
   parser.add_argument('--mean',   action="store_true",default=False,help="if mean flag is set, a mean archive will be created")
   parser.add_argument('meshfile',   type=str,help="NEMO mesh file in netcdf format")
   parser.add_argument('grid2dfile', type=str, nargs="+",help="NEMO 2D data file in netcdf format")
   parser.add_argument('--iexpt',    type=int,default=10,  help="    ")
   parser.add_argument('--makegrid',    type=int,  help="    ")
   parser.add_argument('--iversn',   type=int,default=22,  help="    ")
   parser.add_argument('--yrflag',   type=int,default=3,   help="    ")
   parser.add_argument('--bio_path',   type=str,   help="    ")

   args = parser.parse_args()
   main(args.meshfile,args.grid2dfile,first_j=args.first_j,mean_file=args.mean,iexpt=args.iexpt,makegrid=args.makegrid,iversn=args.iversn,yrflag=args.yrflag,bio_path=args.bio_path)
