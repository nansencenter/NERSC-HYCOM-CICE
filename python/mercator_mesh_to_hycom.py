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

   #2) Check if periodic in i. Done by checking if there is a wrap-around in mesh
   periodic_i=False
   for j in range(jdm) :
      periodic_i=periodic_i or (plon[j,0] == plon[j,-2]  and plat[j,0] == plat[j,-2])

   # Check if arctic patch. Currently this is done by checking if there is a "open" grid cell along top boundary
   arctic_patch=False
   for i in range(idm) :
      arctic_patch=arctic_patch or depth[-1,i] > 0.
   logger.info("Grid periodic in i  :%s"%str(periodic_i))
   logger.info("Arctic Patch        :%s"%str(arctic_patch))

   # Consistency check for arctic patch
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

   # Procedure:
   #  1) Neglect bottom row
   #  2) Shift u points 1 cell to left 
   #  3) Shift v points 1 cell to left 
   #  4) Shift q points 1 cell to left and 1 cell down
   #
   # TODO: Avoid neglecting bottom row by extrapolating qlon and vlon


   # 
   if not periodic_i :
      msg="Only periodic in i supported for now"
      logger.error(msg)
      raise ValueError,msg
   else :

      # remove wrapping for periodic grid
      idm=idm-2
      plon =plon[:,2:]
      plat =plat[:,2:]
      scpx = ncid.variables["e1t"][0,:,2:]
      scpy = ncid.variables["e2t"][0,:,2:]

      # u shifted one grid cell to the left in hycom (rel to nemo)
      ulon = ncid.variables["glamu"][0,:,1:-1]
      ulat = ncid.variables["gphiu"][0,:,1:-1]
      scux = ncid.variables["e1u"][0,:,1:-1]
      scuy = ncid.variables["e2u"][0,:,1:-1]

      # v shifted one grid cell down in hycom (rel to nemo). Extrapolated below
      vlon = ncid.variables["glamv"][0,:,2:]
      vlat = ncid.variables["gphiu"][0,:,2:]
      scvx = ncid.variables["e1v"][0,:,2:]
      scvy = ncid.variables["e2v"][0,:,2:]
      dlon=numpy.mod(vlon[0,:] - numpy.roll(vlon[0,:],1,axis=0) + 360.,360.)
      dlat=vlat[1,:]-vlat[0,:]
      # TODO: Make more generic approach if necessary
      if  not (dlon.min() == dlon.max() and dlat.min() == dlat.max()) :
         msg="Lower row extrapolation assumes constant lat and regular lon spacing."
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
      # TODO: extrapolate
      qlon = ncid.variables["glamf"][0,0:,1:-1]
      qlat = ncid.variables["gphif"][0,0:,1:-1]
      scqx = ncid.variables["e1f"][0,0:,1:-1]
      scqy = ncid.variables["e2f"][0,0:,1:-1]
      dlon=numpy.mod(qlon[0,:] - numpy.roll(qlon[0,:],1,axis=0) + 360.,360.)
      dlat=qlat[1,:]-qlat[0,:]
      # TODO: Make more generic approach if necessary
      if  not (dlon.min() == dlon.max() and dlat.min() == dlat.max()) :
         msg="Lower row extrapolation assumes constant lat and regular lon spacing."
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
      ulon_lft = numpy.copy(ulon)
      ulat_lft = numpy.copy(ulat)
      ulon_rgt = numpy.roll(ulon,1,axis=1) 
      ulat_rgt = numpy.roll(ulat,1,axis=1) 
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

      ddict["scpx"]=scpx
      ddict["scpy"]=scpy
      ddict["scux"]=scux
      ddict["scuy"]=scuy
      ddict["scvx"]=scvx
      ddict["scvy"]=scvy
      ddict["scqx"]=scqx
      ddict["scqy"]=scqy

      ddict["cori"]=corio
      ddict["pang"] =pang
      ddict["pasp"] =asp

      abfile.write_regional_grid(ddict)


      # Total depth calc
      hdept  = ncid.variables["hdept"] [0,1:,1:]
      hdepw  = ncid.variables["hdepw"] [0,1:,1:]
      mbathy = ncid.variables["mbathy"][0,1:,1:]
      e3t_ps = ncid.variables["e3t_ps"][0,1:,1:]
      e3w_ps = ncid.variables["e3w_ps"][0,1:,1:]
      nav_lev= ncid.variables["nav_lev"][:]
      gdept_0= ncid.variables["gdept_0"][0,:]
      gdepw_0= ncid.variables["gdepw_0"][0,:]
      e3t_0  = ncid.variables["e3t_0"][0,:]
      itest=0
      jtest=300


      # KAL - my conclusions  (Double check with Gilles)
      #   hdepw = Total water depth in t-cell
      #   gdepw_0[mbathy-1] Gives Water depth of full vertical cells 1 up tu and including cell mbathy-1
      #   e3t_ps            Is the fraction of cell "mbathy" that is used in the calculations'

      tmp1 = numpy.where(mbathy>0,gdepw_0[mbathy-1] +  e3t_ps,0)
      tmp2 = numpy.where(mbathy>0,hdepw,0)

      print tmp1[jtest,itest],tmp2[jtest,itest]
      print numpy.abs(tmp1-tmp2).max()
      #abfile.write_bathymetry("bathy",1,w5,cutoff)

      print numpy.count_nonzero(numpy.where(numpy.abs(tmp2-tmp1)>1e-1))


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



      print tmp2-tmp1
      figure = matplotlib.pyplot.figure(figsize=(8,8))
      ax=figure.add_subplot(111)
      P=ax.pcolormesh(tmp2-tmp1,vmin=-1e-2,vmax=1e-2)
      figure.colorbar(P)
      figure.canvas.print_figure("tst2.png")

      abfile.write_bathymetry("bathy",1,tmp2,0.)








#   gfile=abfile.ABFileGrid("regional.grid","r")
#   plon=gfile.read_field("plon")
#   plat=gfile.read_field("plat")
#   scpx=gfile.read_field("scpx")
#   scpy=gfile.read_field("scpy")
#   width=numpy.median(scpx)
#   logger.info("Grid median resolution:%8.2f km "%(width/1000.))
#
#   # GEBCO only - TODO: move logic to gebco set
#   if resolution is None  :
#      if width  > 20 :
#         dfile="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D_median20km.nc"
#      elif width > 8 :
#         dfile="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D_median8km.nc"
#      elif width > 4 :
#         dfile="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D_median4km.nc"
#      else :
#         dfile="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D.nc"
#      logger.info ("Source resolution not set - choosing datafile %s"%dfile)
#   else :
#      dfile="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D_median%dkm.nc" % resolution
#      logger.info ("Source resolution set to %d - trying to use datafile %s"%dfile)
#   gebco = modeltools.forcing.bathy.GEBCO2014(filename=dfile)
#   w2=gebco.regrid(plon,plat,width=width)
#   #print w2.min(),w2.max()
#
#   # Run shapiro filter on interpolated data to remove 2 DeltaX noise
#   w3=numpy.copy(w2)
#   for i in range(shapiro_passes):
#      logger.info("Shapiro filter ... pass %d"%(i+1))
#      w3=modeltools.tools.shapiro_filter(w3,threshold=cutoff)
#   #print w3.min(),w3.max()
#
#   # Modify basin 
#   w4=numpy.copy(-w3)
#   it=1
#   while it==1 or numpy.count_nonzero(w4-w4old) > 0 :
#      w4old = numpy.copy(w4)
#      logger.info("Basin modifications ... pass %d"%(it))
#      w4=modeltools.tools.remove_isolated_basins(plon,plat,w4,blo,bla,threshold=cutoff)
#      w4=modeltools.tools.remove_islets(w4,threshold=cutoff)
#      w4=modeltools.tools.remove_one_neighbour_cells(w4,threshold=cutoff)
#      logger.info("Modified %d points "%numpy.count_nonzero(w4-w4old) )
#      it+=1
#   w5=numpy.copy(w4)
#   #print w5.min(),w5.max()
#
#   # Mask data where depth below threshold
#   w5=numpy.ma.masked_where(w5<=cutoff,w5)
#
#   # Create netcdf file with all  stages for analysis
#   logger.info("Writing bathymetry to file hycom_bathymetry.nc")
#   ncid = netCDF4.Dataset("hycom_bathymetry.nc","w")
#   ncid.createDimension("idm",w5.shape[1])
#   ncid.createDimension("jdm",w5.shape[0])
#   ncid.createVariable("lon","f8",("jdm","idm"))
#   ncid.createVariable("lat","f8",("jdm","idm"))
#   ncid.createVariable("h_1","f8",("jdm","idm"))
#   ncid.createVariable("h_2","f8",("jdm","idm"))
#   ncid.createVariable("h_3","f8",("jdm","idm"))
#   ncid.variables["lon"][:]=plon
#   ncid.variables["lat"][:]=plat
#   ncid.variables["h_1"][:]=w2
#   ncid.variables["h_2"][:]=w3
#   ncid.variables["h_3"][:]=w5
#   ncid.close()
#
#   # Print to HYCOM and CICE bathymetry files
#   abfile.write_bathymetry("bathy",1,w5,cutoff)
#   
#   # Show some grid statistics 
#   sx = (w2[1:,:-1]-w2[:-1,:-1])/scpy[1:,:-1]
#   sy = (w2[:-1,1:]-w2[:-1,:-1])/scpx[:-1,1:]
#   grad = sx + sy
#   slopefac=numpy.sqrt(sx**2+sy**2)
#   logger.info("Maximum slope factor after interpolation =%.4f"%slopefac.max())
#
#   sx = (w5[1:,:-1]-w5[:-1,:-1])/scpy[1:,:-1]
#   sy = (w5[:-1,1:]-w5[:-1,:-1])/scpx[:-1,1:]
#   grad = sx + sy
#   slopefac=numpy.sqrt(sx**2+sy**2)
#   logger.info("Maximum slope factor after smoothing    =%.4f"%slopefac.max())
#
#   
#   figure = matplotlib.pyplot.figure(figsize=(8,8))
#   ax=figure.add_subplot(111)
#   P=ax.pcolor(slopefac)
#   figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=w5.min(), vmax=w5.max()))
#   #ax.contour(w5)#,[-10.,-100.,-500.,-1000.])
#   #ax.set_title("Slope fac in color, depth contours in black")
#   logger.info("Slope factor in slopefac.png")
#   figure.canvas.print_figure("slopefac.png")
#


if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('meshfile', type=str,help="NEMO grid mesh file")
   args = parser.parse_args()

   main(args.meshfile)
   

