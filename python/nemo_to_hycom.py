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
import os.path
import time
import cfunits

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

#def p_azimuth(plon1,plat1,plon2,plat2) :
#   tmp = fwd_azimuth( plon1,plat1,plon2,plat2)
#   # Compass angle -> angle rel lat line
#   tmp = numpy.radians(90-tmp)
#   # Put in range [-pi,pi)
#   tmp = numpy.where(tmp >= numpy.pi,tmp-2*numpy.pi,tmp) # Safest option
#   tmp = numpy.where(tmp < -numpy.pi,tmp+2*numpy.pi,tmp) # Safest option
#   #mult = max(-tmp.min()/numpy.pi,1)
#   #mult = 2*(mult/2) + 1                                # use any odd number as long as tmp+mult*numpy.pi > 0...
#   #tmp = numpy.fmod(tmp+9*numpy.pi,2*numpy.pi)-numpy.pi # use any odd number as long as tmp+mult*numpy.pi > 0...
#   return tmp
#
#def fwd_azimuth(lon1,lat1,lon2,lat2) : 
#   geod=pyproj.Geod(ellps="sphere")
#   az1,az2,dist = geod.inv(lon1,lat1,lon2,lat2)
#   return az1
#
#
#def haversine(lon1,lat1,lon2,lat2) :
#   deg2rad = numpy.pi / 180.
#   dlon=lon2-lon1
#   dlat=lat2-lat1
#   dlon=dlon*deg2rad
#   dlat=dlat*deg2rad
#   a=numpy.sin(dlat/2.)**2 + numpy.cos(lat1*deg2rad) * numpy.cos(lat2*deg2rad) * numpy.sin(dlon/2.)**2
#   #c=2.*numpy.arctan2(numpy.sqrt(a),numpy.sqrt(1.-a))
#   c=2.*numpy.arcsin(numpy.sqrt(a))
#   d=6371000.*c
#   return d


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


def periodic_i_shift_right(field,istep) :
   # shift field left by istep steps
   field2  = numpy.roll(field,istep,axis=1)
   return field2

   

def main(filemesh,grid2dfiles) :

   fnametemplate="archv.%Y_%j_%H"
   itest=1
   jtest=200

   ncidmesh=netCDF4.Dataset(filemesh,"r")
   gphiu  = ncidmesh["gphiu"][0,:,1:-1]     # Depth of t points
   glamu  = ncidmesh["glamu"][0,:,1:-1]     # Depth of t points
   gdept  = ncidmesh["gdept_0"][0,:]        # Depth of t points
   gdepw  = ncidmesh["gdepw_0"][0,:]        # Depth of w points
   e3t_ps = ncidmesh["e3t_ps"] [0,:,1:-1]   # Partial steps of t cell
   e3w_ps = ncidmesh["e3w_ps"] [0,:,1:-1]   # Partial steps of w cell
   mbathy = ncidmesh["mbathy"] [0,:,1:-1]   # bathy index
   hdepw  = ncidmesh["hdepw"]  [0,:,1:-1]   # Total depth of w points
   mbathy = mbathy -1                       # python indexing starts from 0

   depth     = hdepw
   depthjp1  = arctic_patch_shift_down(depth ,1)    # nemo values at cell j+1
   mbathyjp1 = arctic_patch_shift_down(mbathy,1)    # nemo values at cell j+1 
   depthip1  = periodic_i_shift_right(depth ,-1)    # nemo values at cell i+1
   mbathyip1 = periodic_i_shift_right(mbathy,-1)    # nemo values at cell i+1
   
   # Estimate partial steps and layer in u-points
   depthu  =numpy.minimum(depth , depthip1)
   mbathy_u=numpy.minimum(mbathy,mbathyip1)
   e3u_ps  =numpy.zeros(depthu.shape)
   J,I= numpy.where(mbathy_u>-1)
   e3u_ps[J,I] = depthu[J,I] - gdepw[mbathy_u[J,I]]
   
   # Estimate partial steps and layer in v-points
   depthv  =numpy.minimum(depth ,depthjp1)
   mbathy_v=numpy.minimum(mbathy,mbathyjp1)
   e3v_ps  =numpy.zeros(depthv.shape)
   I= numpy.where(mbathy_v>-1)
   e3v_ps[I] = depthv[I] - gdepw[mbathy_v[I]]

   # Thickness of t layers (NB: 1 less than gdepw dimension)
   dt = gdepw[1:] - gdepw[:-1]
   ncidmesh.close()


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
      logger.info("gridT  file: %s"%files)
      ncids=netCDF4.Dataset(files,"r")
      s = numpy.copy(ncids.variables["vosaline"][0,:,:,1:-1])   # Costly, make more efficient if needed
      s = numpy.where(s<1e30,s,0.)
      s = numpy.where(s==ncids.variables["vosaline"]._FillValue,0.,s)
      ncidt=netCDF4.Dataset(filet,"r")
      t = numpy.copy(ncidt.variables["votemper"][0,:,:,1:-1])   # Costly, make more efficient if needed
      t = numpy.where(t==ncidt.variables["votemper"]._FillValue,0.,t)
      t = numpy.where(t<1e30,t,0.)

      # time from gridT file. 
      time = ncidt.variables["time_counter"][0]
      tunit = ncidt.variables["time_counter"].units
      tmp=cfunits.Units(tunit)
      refy, refm, refd=(1958,1,1)
      tmp2=cfunits.Units("seconds since %d-%d-%d 00:00:00"%(refy,refm,refd))            # Units from CF convention
      tmp3=cfunits.Units.conform(time,tmp,tmp2)                                         # Transform to new new unit 
      mydt = datetime.datetime(refy,refm,refd,0,0,0) + datetime.timedelta(seconds=tmp3) # Then calculate dt. Phew!




      # Read and calculculate U in hycom U-points. Calculate barotropic U
      logger.info("gridU  file: %s"%fileu)
      ncidu=netCDF4.Dataset(fileu,"r")
      u = numpy.copy(ncidu.variables["vozocrtx"][0,:,:,1:-1])   # Costly, make more efficient if needed
      # TODO: Mid-layer depths seem to be undefined - figure out why ...
      u = numpy.where(u<1e10,u,0.)
      u = numpy.where(u>-1e10,u,0.) # Hmm..
      usum=numpy.zeros(mbathy.shape)
      dsum=numpy.zeros(mbathy.shape)
      for k in range(u.shape[0]-1) : # Dont include lowest layer
         logger.debug("k=%3d, u=%10.3g, mbathy_u[jtest,itest]=%3d,gdepw[k]=%8.2f, depthu[jtest,itest]=%8.2f"%(
            k,u[k,jtest,itest],mbathy_u[jtest,itest],gdepw[k],depthu[jtest,itest]))
         J,I = numpy.where(mbathy_u>k) 
         usum[J,I] = usum[J,I] + u[k,J,I]*dt[k]
         dsum[J,I] = dsum[J,I] + dt[k]
      J,I = numpy.where(mbathy>=0)
      usum[J,I] = usum[J,I] + u[mbathy_u[J,I],J,I]*e3u_ps[J,I]
      dsum[J,I] = dsum[J,I] + e3u_ps[J,I]
      ubaro=numpy.where(dsum>0.1,usum/dsum,0.)

      # Read and calculculate V in hycom V-points. Calculate barotropic V
      logger.info("gridV  file: %s"%filev)
      ncidv=netCDF4.Dataset(filev,"r")
      v = numpy.copy(ncidv.variables["vomecrty"][0,:,:,1:-1])   # Costly, make more efficient if needed
      # TODO: Mid-layer depths seem to be undefined - figure out why ...
      v = numpy.where(v<1e30,v,0.)
      v = numpy.where(v>-1e10,v,0.) # Hmm..
      vsum=numpy.zeros(mbathy.shape)
      dsum=numpy.zeros(mbathy.shape)
      for k in range(v.shape[0]-1) : # Dont include lowest layer
         logger.debug("k=%3d, v=%10.3g, mbathy_v[jtest,itest]=%3d,gdepw[k]=%8.2f, depthv[jtest,itest]=%8.2f"%(
            k,v[k,jtest,itest],mbathy_v[jtest,itest],gdepw[k],depthv[jtest,itest]))
         J,I = numpy.where(mbathy_v>k) 
         vsum[J,I] = vsum[J,I] + v[k,J,I]*dt[k]
         dsum[J,I] = dsum[J,I] + dt[k]
      J,I = numpy.where(mbathy_u>=0)
      vsum[J,I] = vsum[J,I] + v[mbathy_u[J,I],J,I]*e3v_ps[J,I]
      dsum[J,I] = dsum[J,I] + e3v_ps[J,I]
      vbaro=numpy.where(dsum>.1,vsum/dsum,0.)

      # Masks (land:True)
      #print mbathy.min(),mbathy.max()
      ip = mbathy   == -1
      iu = mbathy_u == -1
      iv = mbathy_v == -1
      iu = periodic_i_shift_right(iu,1)   # u: nemo in cell i is hycom in cell i+1
      iv = arctic_patch_shift_up(iu,1)    # v: nemo in cell j is hycom in cell j+1
      #raise NameError,"test"

      # 2D data
      ncid2d=netCDF4.Dataset(file2d,"r")
      ssh          = ncid2d.variables["sossheig"][0,:,1:-1]
      ssh = numpy.where(ssh==ncid2d.variables["sossheig"]._FillValue,0.,ssh)
      ssh = numpy.where(ssh>1e30,0.,ssh) # Hmmmmm
      bar_height   = ncid2d.variables["sobarhei"][0,:,1:-1] 
      dyn_height   = ncid2d.variables["sodynhei"][0,:,1:-1] 
      montg1       = ssh * 9.81  #* 1e-3  # Approx
      logger.warning("TODO:montg pot calculation must be checked...")

      # Write to abfile
      outfile = abfile.ABFileArchv(mydt.strftime(fnametemplate),"w",iexpt=10,iversn=22,yrflag=3,)
      logger.info("Writing 2D variables")
      outfile.write_field(montg1,                          ip,"montg1"  ,0,0,1,0)
      outfile.write_field(ssh,                             ip,"srfhgt"  ,0,0,0,0)
      outfile.write_field(numpy.zeros(ssh.shape),          ip,"surflx"  ,0,0,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),          ip,"salflx"  ,0,0,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),          ip,"bl_dpth" ,0,0,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),          ip,"mix_dpth",0,0,0,0) # Not used
      outfile.write_field(periodic_i_shift_right(ubaro,1) ,iu,"u_btrop" ,0,0,0,0) # u: nemo in cell i is hycom in cell i+1
      outfile.write_field(arctic_patch_shift_down(vbaro,1),iv,"v_btrop" ,0,0,0,0) # v: nemo in cell j is hycom in cell j+1
      #outfile.close() ; raise NameError,"test"
      for k in numpy.arange(u.shape[0]) :
         if k%10==0 : logger.info("Writing 3D variables, level %d of %d"%(k+1,u.shape[0]))
         ul = numpy.squeeze(u[k,:,:]) - ubaro # Baroclinic velocity
         vl = numpy.squeeze(v[k,:,:]) - vbaro # Baroclinic velocity
         sl = numpy.squeeze(s[k,:,:])
         tl = numpy.squeeze(t[k,:,:])

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

         onem=9806.
         outfile.write_field(periodic_i_shift_right(ul,1) ,iu,"u-vel.",0,0,k+1,0) # u: nemo in cell i is hycom in cell i+1
         outfile.write_field(arctic_patch_shift_down(vl,1),iv,"v-vel.",0,0,k+1,0) # v: nemo in cell j is hycom in cell j+1
         outfile.write_field(dtl*onem                     ,ip,"thknss",0,0,k+1,0)
         outfile.write_field(sl                           ,ip,"salin" ,0,0,k+1,0)
         outfile.write_field(tl                           ,ip,"temp"  ,0,0,k+1,0)


      # TODO: Process ice data
      ncid2d.close()
      ncids.close()
      ncidt.close()
      ncidu.close()
      ncidv.close()
      outfile.close()

      logger.info("Finished writing %s.[ab] "% mydt.strftime(fnametemplate))





if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('meshfile',   type=str)
   parser.add_argument('grid2dfile', type=str, nargs="+")
   args = parser.parse_args()
   main(args.meshfile,args.grid2dfile)


