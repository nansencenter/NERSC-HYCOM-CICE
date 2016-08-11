#!/usr/bin/env python
import argparse
import modeltools.tools
import modeltools.nemo
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
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



def main(meshfile,first_j=0) :

   nemo_mesh=modeltools.nemo.NemoMesh(meshfile,first_j=first_j)

#   ncid =netCDF4.Dataset(meshfile,"r")
#   plon =ncid.variables["glamt"][0,:,:]
#   plat =ncid.variables["gphit"][0,:,:]
#   ulon =ncid.variables["glamu"][0,:,:]
#   ulat =ncid.variables["gphiu"][0,:,:]
#   vlon =ncid.variables["glamv"][0,:,:]
#   vlat =ncid.variables["gphiv"][0,:,:]
#   jdm,idm=plon.shape
#
#   # These tests are on full grid 
#   periodic_i = modeltools.nemo.periodic_i(plon,plat)
#   arctic_patch=modeltools.nemo.arctic_patch(plon,plat)
#   logger.info("Grid periodic in i  :%s"%str(periodic_i))
#   logger.info("Arctic Patch        :%s"%str(arctic_patch))
#
#   # Not difficult to extend, but for now...
#   if not periodic_i and arctic_patch:
#      msg="Only periodic in i  and arctic patchsupported for now"
#      logger.error(msg)
#      raise ValueError,msg

   #  HYCOM stagger with same i,j index:
   #  -------------
   #       
   #  x----x----x
   #  |    |    |
   #  |    |    |
   #  U----P----x
   #  |    |    |
   #  |    |    |
   #  Q----V----x
   #
   # Procedure:
   #  - Shift u points 1 cell to left 
   #  - Shift v points 1 cell down
   #  - Shift q points 1 cell to left and 1 cell down

   # Get slices for the unique domain of the nemo mesh
   #si = nemo_mesh.slicei
   #sj = nemo_mesh.slicej

   # Now acquire the data. P-cell data
   plon = nemo_mesh.sliced(nemo_mesh["glamt"][0,:,:])
   plat = nemo_mesh.sliced(nemo_mesh["gphit"][0,:,:])
   scpx = nemo_mesh.sliced(nemo_mesh["e1t"]  [0,:,:])
   scpy = nemo_mesh.sliced(nemo_mesh["e2t"]  [0,:,:])

   # U-cell data. 
   ulon = nemo_mesh.sliced(nemo_mesh.u_to_hycom_u( nemo_mesh["glamu"][0,:,:]))
   ulat = nemo_mesh.sliced(nemo_mesh.u_to_hycom_u( nemo_mesh["gphiu"][0,:,:]))
   scux = nemo_mesh.sliced(nemo_mesh.u_to_hycom_u( nemo_mesh["e1u"  ][0,:,:]))
   scuy = nemo_mesh.sliced(nemo_mesh.u_to_hycom_u( nemo_mesh["e2u"  ][0,:,:]))

   # V-cell data.
   # TODO: Proper extrapolation of data on grid edges (mainly bottom row)
   vlon = nemo_mesh.sliced(nemo_mesh.v_to_hycom_v( nemo_mesh["glamv"][0,:,:]))
   vlat = nemo_mesh.sliced(nemo_mesh.v_to_hycom_v( nemo_mesh["gphiu"][0,:,:]))
   scvx = nemo_mesh.sliced(nemo_mesh.v_to_hycom_v( nemo_mesh["e1v"]  [0,:,:]))
   scvy = nemo_mesh.sliced(nemo_mesh.v_to_hycom_v( nemo_mesh["e2v"]  [0,:,:]))

   # Q-cell data
   # TODO: Proper extrapolation of data on grid edges (mainly bottom row)
   qlon = nemo_mesh.sliced(nemo_mesh.f_to_hycom_q( nemo_mesh["glamf"][0,:,:] ))
   qlat = nemo_mesh.sliced(nemo_mesh.f_to_hycom_q( nemo_mesh["gphif"][0,:,:] ))
   scqx = nemo_mesh.sliced(nemo_mesh.f_to_hycom_q( nemo_mesh["e1f"  ][0,:,:] ))
   scqy = nemo_mesh.sliced(nemo_mesh.f_to_hycom_q( nemo_mesh["e2f"  ][0,:,:] ))

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
   hdepw  = nemo_mesh.sliced(nemo_mesh["hdepw"] [0,:,:])
   mbathy = nemo_mesh.sliced(nemo_mesh["mbathy"][0,:,:])
   tmp2 = numpy.where(mbathy>0,hdepw,0)
   abfile.write_bathymetry("bathy",1,tmp2,0.)

   # Total depth calc
   hdept  = nemo_mesh.sliced(nemo_mesh["hdept"] [0,:,:])
   hdepw  = nemo_mesh.sliced(nemo_mesh["hdepw"] [0,:,:])
   e3t_ps = nemo_mesh.sliced(nemo_mesh["e3t_ps"][0,:,:])
   e3w_ps = nemo_mesh.sliced(nemo_mesh["e3w_ps"][0,:,:])
   nav_lev= nemo_mesh["nav_lev"]  [:]
   gdept_0= nemo_mesh["gdept_0"][0,:]
   gdepw_0= nemo_mesh["gdepw_0"][0,:]
   e3t_0  = nemo_mesh["e3t_0"]  [0,:]
   tmp1 = numpy.where(mbathy>0,gdepw_0[mbathy-1] +  e3t_ps,0)
   itest=0
   jtest=300


   # KAL - my conclusions  (Double check with Gilles)
   #   hdepw = Total water depth in t-cell
   #   gdepw_0[mbathy-1] Gives Water depth of full vertical cells 1 up to and including cell mbathy-1
   #   e3t_ps            Is the fraction of cell "mbathy" that is used in the calculations'
   #   hdepw = gdepw_0[mbathy-1] + e3t_ps

   # Diag plot of depths
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   ax.hold(True)
   ax.plot(nav_lev,"r-",label="nav_lev")
   ax.plot(gdept_0,"b.",label="gdept_0")
   ax.plot(gdepw_0,"g.",label="gdepw_0")
   ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*hdepw[jtest,itest],"c",lw=2,label="hdepw")
   ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*hdept[jtest,itest],"m",lw=2,label="hdept")
   ax.plot(numpy.arange(nav_lev.size),numpy.ones((nav_lev.size,))*gdepw_0[mbathy[jtest,itest]-1]+e3t_ps[jtest,itest],"m.",lw=2,label="gdepw_0[mbathy[jtest,itest]-1]+e3t_ps")
   ax.legend()
   ax.set_ylim(hdepw[jtest,itest]-600,hdepw[jtest,itest]+600)
   figure.canvas.print_figure("tst.png")

   # Diag plot of depths
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=ax.pcolormesh(tmp2-tmp1,vmin=-1e-2,vmax=1e-2)
   figure.colorbar(P)
   figure.canvas.print_figure("tst2.png")










if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--first_j',   type=int,default=0)
   parser.add_argument('meshfile', type=str,help="NEMO grid mesh file")
   args = parser.parse_args()

   main(args.meshfile,first_j=args.first_j)


