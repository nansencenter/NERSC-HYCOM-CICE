#!/usr/bin/env python 

# Name: nemo2archvz.py
# Purpose: Convert MMERCATOR data to HYCOM archive files in z coordinates
# Author: Mostafa Bakhoday-Paskyabi (Mostafa.Bakhoday@nersc.no)
# Created: 8 December 2017
# Copyright: (c) NERSC Norway 2017
# Licence:
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
# http://www.gnu.org/licenses/gpl-3.0.html
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# USAGE:
#  WITH BIO
#
# ../bin/nemo2archvz_regular.py ../../MERCATOR_REGULAR/GLO_MFC_001_24_MESH.nc ../../MERCATOR_REGULAR/MERCATOR-PHY-24-2018-01-01-12.nc  --bio_path ./
# ../bin/nemo2archvz_regular.py PATH-ON-FRAM/GLO_MFC_001_24_MESH.nc PHY-NCFILE-WITH-NC-EXTENSION  --bio_path BIO-FILE-ONLY-PATH
#
# WITHOUT BIO
# ../bin/nemo2archvz_regular.py ../../MERCATOR_REGULAR/GLO_MFC_001_24_MESH.nc ../../MERCATOR_REGULAR/MERCATOR-PHY-24-2018-01-01-12.nc
# ../bin/nemo2archvz_regular.py PATH-ON-FRAM/GLO_MFC_001_24_MESH.nc PHY-NCFILE-WITH-NC-EXTENSION


# Refrences:
# G. Madec and the NEMO team, 2012, NEMO ocean engine.
# M. Bakhoday-Paskyabi et al., 2017 (under preparation), Effects of nesting and open boundary conditions: a comparative study between TOPAZ4 and TOPAZ5 systems
# M. Bakhoday-Paskyabi, 11 July 2019, adding bio fields
#
from   matplotlib import pyplot as plt
import abfile
import numpy
import numpy.ma as ma
#from   mpl_toolkits.basemap import Basemap, cm
import sys
import os.path
import re
import argparse
import netCDF4
import cfunits
import time
import datetime
import modeltools.forcing.bathy
import modeltools.nemo
import logging
import glob


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


# input variables
grav  = 9.81                     # acceleration due to gravity (m.s-2)
onem  = 9806.
spval = 2.0**100.0
file_pre="global_analysis_forecast_phy_"

fillgap_method   = 1             # (1) using maplev for temperature & salinity; (2) partial use of maplev
thickness_method = 2             # (1) based on depth; and (2) based on partial steps
mbathy_method    = 1             # (1) using mbathy; (2) using mask
timeavg_method   = 0             # (1) time average of two consecutive netcdf files; and other values no temporal averaging.

# Note that using fill-gapping method will induce small discprepancies with original quantity (in the order of 1e-7, 1e-6, 1e-8 for temperature, salinity, and velocity, respectively). The discrepency is substantailly reduced if the value of "lpp" increase to more than, for example, 30. However, the method is more reliable than replacing the bad-flagged points with some constant values.
#


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



def calc_uvbaro(uo,vo,e3t,iu,iv) :
    #
    # Calculate barotropic velocity. uo and vo
    # are 3D velocities in P-grid positions.
    #
    nlev=uo.shape[0]
    utemp=numpy.zeros(uo.shape)
    vtemp=numpy.zeros(uo.shape)
    ubaro=numpy.zeros(uo.shape[-2:])
    vbaro=numpy.zeros(uo.shape[-2:])
    
    for k in range(nlev) :
        
        utmp = uo[k,:,:]
        vtmp = vo[k,:,:]
        utmp = numpy.where(abs(utmp)>10,0.,utmp)
        vtmp = numpy.where(abs(vtmp)>10,0.,vtmp)
        utemp[k,:,:]=e3t[k,:,:]*utmp[:,:]
        vtemp[k,:,:]=e3t[k,:,:]*vtmp[:,:]
    #
    dsum=numpy.nansum(e3t,axis=0)
    # Avoid divide by zero warning
    dsum=numpy.where(abs(dsum)<1e-2,0.05,dsum)
    ubaro=numpy.nansum(utemp,axis=0)/dsum
    vbaro=numpy.nansum(vtemp,axis=0)/dsum
    # convert to U- and V-grid
    ubaro = p2u_2d(numpy.where(abs(ubaro)>10,0.,ubaro))
    vbaro = p2v_2d(numpy.where(abs(vbaro)>10,0.,vbaro))
    # Land mask on U- and V-grid points 
    ubaro[iu]=spval
    vbaro[iv]=spval
    return ubaro,vbaro


def p2u_2d(var_p) :
    # Locationa of grid points in MERCATOR products
    # are defined on T-point (in Arakawa C-grid arcitecture)
    # here we assume T-grid of the NEMO outputs on regular grid
    # corresponds to P-grid of HYCOM, and convert velocity
    # field to the HYCOM u-grid (first row needs to be extrapolated).
    #
    var_v=numpy.zeros(var_p.shape)
    var_v[1:,:]=0.5*(var_p[1:,:]+var_p[:-1,:])
    # extrapolation
    var_v[0,:] = 2.0*var_v[1,:] - var_v[2,:]
    return var_v

def p2v_2d(var_p) :
    # Locationa of grid points in MERCATOR products
    # are defined on T-point (in Arakawa C-grid arcitecture)
    # here we assume T-grid of the NEMO outputs on regular grid
    # corresponds to P-grid of HYCOM, and convert velocity
    # field to the HYCOM v-grid (first column needs to be extrapolated).
    #
    var_u=numpy.zeros(var_p.shape)
    var_u[:,1:]=0.5*(var_p[:,:-1]+var_p[:,1:])
    # extrapolation
    var_u[:,0] = 2.0*var_u[:,1] - var_u[:,2]
    return var_u


def read_grid(filemesh) :
    
    ncid0=netCDF4.Dataset(filemesh[:-7]+"COORD.nc","r")
    numpy.seterr(invalid='ignore')
    e3t=ncid0.variables["e3t"][:,:,:]
    ncid0.close()
    
    ncid0=netCDF4.Dataset(filemesh,"r")
    plon=ncid0.variables["lon"][:,:]
    plat=ncid0.variables["lat"][:,:]
    with numpy.errstate(invalid='ignore'):
        hdept=ncid0.variables["Bathymetry"][:,:]
        gdept=ncid0.variables["depth"][:]
        mask=ncid0.variables["mask"][:,:]
        mbathy=numpy.int8(ncid0.variables["mbathy"][:,:])
        mbathy_u=numpy.int8(ncid0.variables["mbathy_u"][:,:])
        mbathy_v=numpy.int8(ncid0.variables["mbathy_v"][:,:])

    mbathy   = mbathy  -1
    mbathy_u = mbathy_u-1
    mbathy_v = mbathy_v-1
    ncid0.close()
    
    return hdept,gdept,mbathy,mbathy_u,mbathy_v,mask,e3t,plon,plat


def main(meshfile,file,iexpt=10,iversn=22,yrflag=3,bio_path=None) :
    
    #
    # Trim input netcdf file name being appropriate for reading
    #
    meshfile=str(meshfile)[2:-2]
    logger.info("Reading mesh information from %s."%(meshfile))
    #
    # Read mesh file containing grid and coordinate information.
    # Note that for now, we are using T-grid in vertical which may need
    # to be improved by utilizing W-point along the vertical axis.
    #
    hdept,gdept,mbathy,mbathy_u,mbathy_v,mask,e3t,plon,plat=read_grid(meshfile)
    logger.warning("Reading grid information from regional.grid.[ab] (not completed)")
    #
    # Convert from P-point (i.e. NEMO grid) to U and V HYCOM grids
    #
    mask_u=p2u_2d(mask)
    mask_v=p2v_2d(mask)
    #
    # Read regional.grid.[ab]
    # Grid angle is not used for this product because all quantities are
    # on regular rectangular grid points.
    #
    angle=numpy.zeros(plon.shape)
    #
    # Number vertical layers in T-point.
    #
    nlev=gdept.size
    #
    # layer thickness in the absence of layer partial steps.
    #
    dt = gdept[1:] - gdept[:-1]
    #
    # Prepare/read input data file (in netcdf format). Reference time is 1950-01-01
    #
    logger.info("Reading data files.")
    file=str(file).strip()[2:-2]
    dirname=os.path.dirname(file)
    logger.debug("file name is {}".format(file))
    logger.debug("dirname is {}".format(dirname))
    logger.debug("basename is {}".format(os.path.basename(file)))
    m=re.match("(MERCATOR-PHY-24-)(.*\.nc)",os.path.basename(file))
    logger.debug("file prefix is {}".format(file_pre))
###    m=re.match(file_pre,os.path.basename(file))
    if not m:
        msg="File %s is not a grid2D file, aborting"%file
        logger.error(msg)
        raise ValueError,msg
    
    #fileinput0=os.path.join(dirname+"/"+"MERCATOR-PHY-24-"+m.group(2))
    file_date=file[-16:-6]
    fileinput0=file
    print file_date,file
    next_day=datetime.datetime.strptime(file_date, '%Y-%m-%d')+datetime.timedelta(days=1)
    fileinput1=datetime.datetime.strftime(next_day,'%Y%m%d')
    fileinput1=os.path.join(dirname+"/"+file_pre+fileinput1+'.nc')
    
    logger.info("Reading from %s"%(fileinput0))
    ncid0=netCDF4.Dataset(fileinput0,"r")
    if timeavg_method==1 and os.path.isfile(fileinput1) :
        
        logger.info("timeavg_method=1, Reading from %s"%(fileinput1))
        ncid1=netCDF4.Dataset(fileinput1,"r")
        #
        # Calculate temporal averaged temperature, salinity, and velocity
        #
        uo =   0.5*(ncid0.variables["uo"][0,:,:,:]+    ncid1.variables["uo"][0,:,:,:])
        vo =   0.5*(ncid0.variables["vo"][0,:,:,:]+    ncid1.variables["vo"][0,:,:,:])
        salt = 0.5*(ncid0.variables["so"][0,:,:,:]+    ncid1.variables["so"][0,:,:,:])
        temp = 0.5*(ncid0.variables["thetao"][0,:,:,:]+ncid1.variables["thetao"][0,:,:,:])
        ssh = numpy.squeeze(0.5*(ncid0.variables["zos"][0,:,:]+ncid1.variables["zos"][0,:,:]))
    
    else:
        #
	# Set variables based on current file when timeavg_method ~=1 or the next netcdf file is not available
        logger.debug("time average method set to {}".format(timeavg_method))
        uo =   ncid0.variables["uo"][0,:,:,:]
        vo =   ncid0.variables["vo"][0,:,:,:]
        salt = ncid0.variables["so"][0,:,:,:]
        temp = ncid0.variables["thetao"][0,:,:,:]
        ssh = numpy.squeeze(ncid0.variables["zos"][0,:,:])
    #
    # I will account these values afterward. Because in the current version, I am accounting for missing values using a gap-filling methodology.
    #	
    logger.debug("getting _FillValue")
    uofill=ncid0.variables["uo"]._FillValue
    vofill=ncid0.variables["vo"]._FillValue
    slfill=ncid0.variables["so"]._FillValue
    tlfill=ncid0.variables["thetao"]._FillValue
    shfill=ncid0.variables["zos"]._FillValue

    # Set time
    logger.info("Set time.")
    time=ncid0.variables["time"][0]
    unit=ncid0.variables["time"].units
    tmp=cfunits.Units(unit)
    refy,refm,refd=(1950,1,1)
    tmp2=cfunits.Units("hours since %d-%d-%d 00:00:00"%(refy,refm,refd))
    tmp3=int(numpy.round(cfunits.Units.conform(time,tmp,tmp2)))
    mydt = datetime.datetime(refy,refm,refd,0,0,0) + datetime.timedelta(hours=tmp3) # Then calculate dt. Phew!

    if timeavg_method==1 and os.path.isfile(fileinput1)  :
        fnametemplate="archv.%Y_%j_%H"
        deltat=datetime.datetime(refy,refm,refd,0,0,0) + \
              datetime.timedelta(hours=tmp3) + \
              datetime.timedelta(hours=12)
        oname=deltat.strftime(fnametemplate)
    else:
        #
        # I am assuming that daily mean can be set at 00 instead of 12
        # for cases that there is no information of next day.
        #
        fnametemplate="archv.%Y_%j"
        deltat=datetime.datetime(refy,refm,refd,0,0,0) + \
              datetime.timedelta(hours=tmp3)
        oname=deltat.strftime(fnametemplate) + '_00'

    # model day
    refy, refm, refd=(1900,12,31)
    model_day= deltat-datetime.datetime(refy,refm,refd,0,0,0)
    model_day=model_day.days
    logger.info("Model day in HYCOM is %s"%str(model_day))
    if bio_path:
       jdm,idm=numpy.shape(plon)
       points = numpy.transpose(((plat.flatten(),plon.flatten())))
       delta = mydt.strftime( '%Y-%m-%d')
       # filename format MERCATOR-BIO-14-2013-01-05-00
       print bio_path,delta
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
       po4 = po4 * 106.0 * 12.01
       si = si * 6.625 * 12.01
       no3 = no3 * 6.625 * 12.01
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




    logger.debug("filename prefix is {}".format(oname))
    logger.debug("Hours = {}".format(tmp3))
    logger.info("Read, trim, rotate NEMO velocities.")
    u=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
    v=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
    utmp=numpy.zeros((mbathy.shape))
    vtmp=numpy.zeros((mbathy.shape))
    #
    # Metrices to detect carrefully bottom at p-, u-, and v-grid points.While I have used 3D, mask data,following methods are good enough for now.
    #
    if mbathy_method  ==  1 :
        ip = mbathy   == -1
        iu = mbathy_u == -1
        iv = mbathy_v == -1
    else:
        ip = mask   == 0
        iu = mask_u == 0
        iv = mask_v == 0
    #
    # Read 3D velocity field to calculate barotropic velocity
    #
    # Estimate barotropic velocities using partial steps along the vertical axis. Note that for the early version of this code, 
    # I used dt = gdept[1:] - gdept[:-1] on NEMO t-grid. Furthermore, you may re-calculate this part on vertical grid cells for future. 
    #
    logger.info("Calculate barotropic velocities.")
    ubaro,vbaro=calc_uvbaro(uo,vo,e3t,iu,iv)
    #
    # Save 2D fields (here only ubaro & vbaro)
    #
    zeros=numpy.zeros(mbathy.shape)
    #flnm = open(oname+'.txt', 'w')
    #flnm.write(oname)
    #flnm.close()
    ssh = numpy.where(numpy.abs(ssh)>1000,0.,ssh*9.81) # NB: HYCOM srfhgt is in geopotential ...
    #
    outfile = abfile.ABFileArchv("./data/"+oname,"w",iexpt=iexpt,iversn=iversn,yrflag=yrflag,)
    outfile.write_field(zeros,                   ip,"montg1"  ,0,model_day,1,0)
    outfile.write_field(ssh,                     ip,"srfhgt"  ,0,model_day,0,0)
    outfile.write_field(zeros,                   ip,"surflx"  ,0,model_day,0,0) # Not used
    outfile.write_field(zeros,                   ip,"salflx"  ,0,model_day,0,0) # Not used
    outfile.write_field(zeros,                   ip,"bl_dpth" ,0,model_day,0,0) # Not used
    outfile.write_field(zeros,                   ip,"mix_dpth",0,model_day,0,0) # Not used
    outfile.write_field(ubaro,                   iu,"u_btrop" ,0,model_day,0,0)
    outfile.write_field(vbaro,                   iv,"v_btrop" ,0,model_day,0,0)
    #
    if bio_path:
       logger.info("Calculate baroclinic velocities, temperature, and salinity data as well as BIO field.")
    else:
       logger.info("Calculate baroclinic velocities, temperature, and salinity data.")
    for k in numpy.arange(u.shape[0]) :
        if bio_path:
           no3k=interpolate2d(blat, blon, no3[k,:,:], points)
           no3k = maplev(numpy.reshape(no3k,[jdm,idm]))
           po4k=interpolate2d(blat, blon, po4[k,:,:], points)
           po4k = maplev(numpy.reshape(po4k,[jdm,idm]))
           si_k=interpolate2d(blat, blon, si[k,:,:], points)
           si_k = maplev(numpy.reshape(si_k,[jdm,idm]))
        #

        #
        uo[k,:,:]=numpy.where(numpy.abs(uo[k,:,:])<10,uo[k,:,:],0)
        vo[k,:,:]=numpy.where(numpy.abs(vo[k,:,:])<10,vo[k,:,:],0)

        # Baroclinic velocity (in HYCOM U- and V-grid)
        ul = p2u_2d(numpy.squeeze(uo[k,:,:])) - ubaro
        vl = p2v_2d(numpy.squeeze(vo[k,:,:])) - vbaro
        ul[iu]=spval
        vl[iv]=spval
        
        # Layer thickness
        
        dtl=numpy.zeros(mbathy.shape)
        # Use dt for the water column except the nearest cell to bottom 
        if thickness_method==1:
            if k < u.shape[0]-1 :
                J,I = numpy.where(mbathy>k)
                e3=(e3t[k,:,:])
                dtl[J,I]=dt[k]
                J,I = numpy.where(mbathy==k)
                dtl[J,I]=e3[J,I]
            else:
                e3=(e3t[k,:,:])
                J,I = numpy.where(mbathy==k)
                dtl[J,I]=e3[J,I]
	# Use partial cells for the whole water column.
        else :
            J,I = numpy.where(mbathy>=k)
            dtl[J,I]=e3t[k,J,I]

        # Salinity
        sl = salt[k,:,:]

        # Temperature
        tl = temp[k,:,:]
        # Need to be carefully treated in order to minimize artifacts to the resulting [ab] files.
        if fillgap_method==1:
            J,I= numpy.where(mbathy<k)
            sl = maplev(numpy.where(numpy.abs(sl)<1e2,sl,numpy.nan))
            sl[J,I]=spval
            J,I= numpy.where(mbathy<k)
            tl = maplev(numpy.where(numpy.abs(tl)<1e2,tl,numpy.nan))
            tl[J,I]=spval
        else:
            sl = numpy.where(numpy.abs(sl)<1e2,sl,numpy.nan)
            sl = numpy.minimum(numpy.maximum(maplev(sl),25),80.)
            tl = numpy.where(numpy.abs(tl)<=5e2,tl,numpy.nan)
            tl = numpy.minimum(numpy.maximum(maplev(tl),-5.),50.)

        # Thickness
        dtl = maplev(dtl)
        if k > 0 :
            with numpy.errstate(invalid='ignore'):
                K= numpy.where(dtl < 1e-4)
            sl[K] = sl_above[K]
            tl[K] = tl_above[K]
            if k%10==0:
                logger.info("Save to the output %s at level %s"%(deltat.strftime(fnametemplate),str(k)))
        #
        sl[ip]=spval
        tl[ip]=spval

        # Save 3D fields
        outfile.write_field(ul      ,iu,"u-vel.",0,model_day,k+1,0)
        outfile.write_field(vl      ,iv,"v-vel.",0,model_day,k+1,0)
        outfile.write_field(dtl*onem,ip,"thknss",0,model_day,k+1,0)
        outfile.write_field(tl      ,ip,"temp" , 0,model_day,k+1,0)
        outfile.write_field(sl      ,ip,"salin" ,0,model_day,k+1,0)
        if bio_path :
           outfile.write_field(no3k      ,ip,"ECO_no3" ,0,model_day,k+1,0)
           outfile.write_field(po4k      ,ip,"ECO_pho" ,0,model_day,k+1,0)
           outfile.write_field(si_k      ,ip,"ECO_sil" ,0,model_day,k+1,0)
                
        tl_above=numpy.copy(tl)
        sl_above=numpy.copy(sl)
    
    outfile.close()
    ncid0.close()
    if timeavg_method==1 and os.path.isfile(fileinput1)  :
        ncid1.close()
    if bio_path :
       ncidb.close()


if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('meshfile',   type=str, nargs="+",  help="    ")
    parser.add_argument('file',       type=str, nargs="+",  help="    ")
    parser.add_argument('--iexpt',    type=int,default=10,  help="    ")
    parser.add_argument('--iversn',   type=int,default=22,  help="    ")
    parser.add_argument('--yrflag',   type=int,default=3,   help="    ")
    parser.add_argument('--bio_path', type=str,             help="    ")
    args = parser.parse_args()
    main(args.meshfile,args.file,iexpt=args.iexpt,iversn=args.iversn,yrflag=args.yrflag,bio_path=args.bio_path)

