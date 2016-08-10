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
import numpy
import netCDF4
import logging
import re
import pyproj

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

def p_azimuth(plon1,plat1,plon2,plat2) :
   tmp = fwd_azimuth( plon1,plat1,plon2,plat2)
   # Compass angle -> angle rel lat line
   tmp = numpy.radians(90-tmp)
   # Put in range [-pi,pi)
   tmp = numpy.where(tmp >= numpy.pi,tmp-2*numpy.pi,tmp) # Safest option
   tmp = numpy.where(tmp < -numpy.pi,tmp+2*numpy.pi,tmp) # Safest option
   #mult = max(-tmp.min()/numpy.pi,1)
   #mult = 2*(mult/2) + 1                                # use any odd number as long as tmp+mult*numpy.pi > 0...
   #tmp = numpy.fmod(tmp+9*numpy.pi,2*numpy.pi)-numpy.pi # use any odd number as long as tmp+mult*numpy.pi > 0...
   return tmp

def fwd_azimuth(lon1,lat1,lon2,lat2) : 
   geod=pyproj.Geod(ellps="sphere")
   az1,az2,dist = geod.inv(lon1,lat1,lon2,lat2)
   return az1


def haversine(lon1,lat1,lon2,lat2) :
   deg2rad = numpy.pi / 180.
   dlon=lon2-lon1
   dlat=lat2-lat1
   dlon=dlon*deg2rad
   dlat=dlat*deg2rad
   a=numpy.sin(dlat/2.)**2 + numpy.cos(lat1*deg2rad) * numpy.cos(lat2*deg2rad) * numpy.sin(dlon/2.)**2
   #c=2.*numpy.arctan2(numpy.sqrt(a),numpy.sqrt(1.-a))
   c=2.*numpy.arcsin(numpy.sqrt(a))
   d=6371000.*c
   return d
   

def main(meshfile) :

   ncid =netCDF4.Dataset(meshfile,"r")
   plon =ncid.variables["glamt"][0,:,:]
   plat =ncid.variables["gphit"][0,:,:]
   depth=ncid.variables["hdept"][0,:,:]
   jdm,idm=plon.shape


   #1) Check if periodic in i. Done by checking if there is a "open" grid cell along left and right boundaries
   #periodic_i=False
   #for j in range(jdm) :
   #   periodic_i=periodic_i or depth[j,0] > 0 or depth[j,-1] > 0

   #2) Check if periodic in i. Done by checking if there is a overlap in mesh. Currently there is a two point  overlap 
   periodic_i=False
   for j in range(jdm) :
      periodic_i=periodic_i or (plon[j,0] == plon[j,-2]  and plat[j,0] == plat[j,-2])

   # Check if arctic patch. Currently this is done by checking if there is a "open" grid cell along top boundary
   arctic_patch=False
   for i in range(idm) :
      arctic_patch=arctic_patch or depth[-1,i] > 0.
   logger.info("Grid periodic in i  :%s"%str(periodic_i))
   logger.info("Arctic Patch        :%s"%str(arctic_patch))

   # Consistency check for arctic patch. Only diag for now. if Max distance is too big, something is wrong ....
   maxdist=-1e8
   mindist=1e8
   if arctic_patch :
      for i in range(idm) :
         i2 = idm - i -1
         #print i,i2,plon[-1,i],plon[-1,i2], plat[-1,i],plat[-1,i2],haversine( plon[-1,i], plat[-1,i],plon[-1,i2],plat[-1,i2])
         d=haversine( plon[-1,i], plat[-1,i],plon[-1,i2],plat[-1,i2])
         maxdist=max(maxdist,d)
         mindist=min(mindist,d)
   logger.info("Arctic Patch maxdist: %s meters across 'cap'"%str(maxdist))
   logger.info("Arctic Patch mindist: %s meters across 'cap'"%str(mindist))

   #  NEMO stagger with same i,j index:
   #  -------------
   #       V----F
   #       |    |
   #       |    |
   #       P----U
   #
   #  HYCOM stagger with same i,j index:
   #  -------------
   #       
   #  U----P
   #  |    |
   #  |    |
   #  Q----V   
   #
   # Procedure:
   #  - Shift u points 1 cell to left 
   #  - Shift v points 1 cell down
   #  - Shift q points 1 cell to left and 1 cell down


   # Must have periodic for now
   if not periodic_i :
      msg="Only periodic in i supported for now"
      logger.error(msg)
      raise ValueError,msg

   # remove wrapping for periodic grid.
   # nemo T at i,j -> HYCOM T at i,j.
   idm=idm-2 
   plon =plon[:,1:-1]       
   plat =plat[:,1:-1]
   scpx = ncid.variables["e1t"][0,:,1:-1]
   scpy = ncid.variables["e2t"][0,:,1:-1]

   # nemo U at i,j -> HYCOM U at i+1,j. 
   ulon = ncid.variables["glamu"][0,:,:-2]
   ulat = ncid.variables["gphiu"][0,:,:-2]
   scux = ncid.variables["e1u"]  [0,:,:-2]
   scuy = ncid.variables["e2u"]  [0,:,:-2]

   # nemo V at i,j -> HYCOM V at i,j+1. V extrapolation below
   vlon = ncid.variables["glamv"][0,:,1:-1]
   vlat = ncid.variables["gphiu"][0,:,1:-1]
   scvx = ncid.variables["e1v"][0,:,1:-1]
   scvy = ncid.variables["e2v"][0,:,1:-1]
   dlon=numpy.mod(vlon[0,:] - numpy.roll(vlon[0,:],1,axis=0) + 360.,360.)
   dlat=vlat[1,:]-vlat[0,:]

   # v extrapolation: TODO: Make more generic approach if necessary
   if  not (dlon.min() == dlon.max() and dlat.min() == dlat.max()) :
      msg="Lower row extrapolation assumes constant lat and regular lon spacing in 1st row."
      logger.error(msg)
      raise ValueError,msg
   else :
      dlon=dlon[0]
      dlat=dlat[0]
      vlon[1:,:] = vlon[:-1,:]
      #
      vlat[1:,:] = vlat[:-1,:]
      vlat[ 0,:] = vlat[  1,:] - dlat
      #
      scvx[1:,:] = scvx[:-1,:]
      scvx[0,:]  = scvx[0,:] * numpy.cos(numpy.rad2deg(vlat[1,0])) / numpy.cos(numpy.rad2deg(vlat[0,0]))  # Correct for latitude 
      #
      scvy[1:,:] = scvy[:-1,:]
      scvy[0,:]  = scvy[1,:]


   # q shifted one grid cell down and one cell to the left in hycom (rel to nemo). Extrapolation below
   qlon = ncid.variables["glamf"][0,0:,2:]
   qlat = ncid.variables["gphif"][0,0:,2:]
   scqx = ncid.variables["e1f"][0,0:,2:]
   scqy = ncid.variables["e2f"][0,0:,2:]
   dlon=numpy.mod(qlon[0,:] - numpy.roll(qlon[0,:],1,axis=0) + 360.,360.)
   dlat=qlat[1,:]-qlat[0,:]

   # q extrapolation. TODO: Make more generic approach if necessary
   if  not (dlon.min() == dlon.max() and dlat.min() == dlat.max()) :
      msg="Lower row extrapolation assumes constant lat and regular lon spacing in 1st row."
      logger.error(msg)
      raise ValueError,msg
   else :
      dlon=dlon[0]
      dlat=dlat[0]
      qlon[1:,:] = qlon[:-1,:]
      #
      qlat[1:,:] = qlat[:-1,:]
      qlat[ 0,:] = qlat[  1,:] - dlat
      #
      scqx[1:,:] = scqx[:-1,:]
      scqx[0,:]  = scqx[0,:] * numpy.cos(numpy.deg2rad(qlat[1,0])) / numpy.cos(numpy.rad2deg(qlat[0,0]))  # Correct for latitude 
      #
      scqy[1:,:] = scqy[:-1,:]
      scqy[0,:]  = scqy[1,:]





   # Angle used for rotation
   #vlon_top = ncid.variables["glamv"][0,:,2:]
   #vlat_top = ncid.variables["gphiu"][0,:,2:]
   #vlon_bot = numpy.copy(vlon)
   #vlat_bot = numpy.copy(vlat)
   #pang = p_azimuth(
   ulon_rgt = numpy.copy(ulon)
   ulat_rgt = numpy.copy(ulat)
   ulon_lft = numpy.roll(ulon,1,axis=1) 
   ulat_lft = numpy.roll(ulat,1,axis=1) 
   pang = p_azimuth(ulon_lft,ulat_lft,ulon_rgt,ulat_rgt)

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


   # Total depth calc
   hdept  = ncid.variables["hdept"] [0,:,1:-1]
   hdepw  = ncid.variables["hdepw"] [0,:,1:-1]
   mbathy = ncid.variables["mbathy"][0,:,1:-1]
   e3t_ps = ncid.variables["e3t_ps"][0,:,1:-1]
   e3w_ps = ncid.variables["e3w_ps"][0,:,1:-1]
   nav_lev= ncid.variables["nav_lev"]  [:]
   gdept_0= ncid.variables["gdept_0"][0,:]
   gdepw_0= ncid.variables["gdepw_0"][0,:]
   e3t_0  = ncid.variables["e3t_0"]  [0,:]
   itest=0
   jtest=300


   # KAL - my conclusions  (Double check with Gilles)
   #   hdepw = Total water depth in t-cell
   #   gdepw_0[mbathy-1] Gives Water depth of full vertical cells 1 up to and including cell mbathy-1
   #   e3t_ps            Is the fraction of cell "mbathy" that is used in the calculations'
   #   hdepw = gdepw_0[mbathy-1] + e3t_ps
   tmp1 = numpy.where(mbathy>0,gdepw_0[mbathy-1] +  e3t_ps,0)
   tmp2 = numpy.where(mbathy>0,hdepw,0)
   #print tmp1[jtest,itest],tmp2[jtest,itest]
   #print numpy.abs(tmp1-tmp2).max()
   #print numpy.count_nonzero(numpy.where(numpy.abs(tmp2-tmp1)>1e-1))

   # Diag plot of depths
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   ax.hold(True)
   ax.plot(nav_lev,"r-",label="nav_lev")
   ax.plot(gdept_0,"b.",label="gdept_0")
   ax.plot(gdepw_0,"g.",label="gdepw_0")
   ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*hdepw[jtest,itest],"c",lw=2,label="hdepw")
   ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*hdept[jtest,itest],"m",lw=2,label="hdept")
   #ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*gdept_0[mbathy[jtest,itest]],"k",lw=2,label="gdept_0[mbathy[jtest,itest]]")
   #ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*gdepw_0[mbathy[jtest,itest]],"y",lw=2,label="gdepw_0[mbathy[jtest,itest]]")
   #ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*gdepw_0[mbathy[jtest,itest]-1],".5",lw=2,label="gdepw_0[mbathy[jtest,itest]-1]")
   ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*gdepw_0[mbathy[jtest,itest]-1]+e3t_ps[jtest,itest],"m.",lw=2,label="gdepw_0[mbathy[jtest,itest]-1]+e3t_ps")
   #ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*gdepw_0[mbathy[jtest,itest]-1]+e3w_ps[jtest,itest],"m.",lw=2,label="gdepw_0[mbathy[jtest,itest]-1]+e3w_ps")
   ax.legend()
   ax.set_ylim(hdepw[jtest,itest]-600,hdepw[jtest,itest]+600)
   figure.canvas.print_figure("tst.png")

   # Diag plot of depths
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolormesh(tmp2-tmp1,vmin=-1e-2,vmax=1e-2)
   figure.colorbar(P)
   figure.canvas.print_figure("tst2.png")

   abfile.write_bathymetry("bathy",1,tmp2,0.)









if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('meshfile', type=str,help="NEMO grid mesh file")
   args = parser.parse_args()

   main(args.meshfile)


