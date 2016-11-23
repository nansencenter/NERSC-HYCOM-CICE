#!/usr/bin/env python
import modeltools.nemo
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
import nemo_mesh_to_hycom
import sys

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


#def arctic_patch_shift_up(field,jstep) :
#   # shift field down
#   if jstep <> 1 :
#      raise NameError,"Arctic_patch_shift only with jstep=1 for now"
#   field2 = numpy.copy(field)
#   field2[1:,:] = field2[0:-1,:] # Shift up
#   # NB:  row 0 is same as row 1 (extrapolated
#   return field2
#
#def arctic_patch_shift_down(field,jstep) :
#   # shift field down
#   if jstep <> 1 :
#      raise NameError,"Arctic_patch_shift only with jstep=1 for now"
#   field2 = numpy.copy(field)
#   field2[0:-1,:] = field2[1:,:] # Shift down
#   tmp=field2[-1,:]              # Top row as top ...
#   field2[-1,:] = tmp[::-1]      # .. but reversed direction
#   return field2
#
#def periodic_i_shift_right(field,istep) :
#   # shift field left by istep steps
#   field2  = numpy.roll(field,istep,axis=1)
#   return field2

   

def main(filemesh,grid2dfiles,first_j=0,mean_file=False) :

   if mean_file :
      fnametemplate="archm.%Y_%j_%H"
   else :
      fnametemplate="archv.%Y_%j_%H"
   itest=1
   jtest=200
   logger.info("Mean file:%s"%str(mean_file))
   logger.info("Output file template:%s"%str(fnametemplate))

   # Write regional files
   nemo_mesh_to_hycom.main(filemesh,first_j=first_j)


   nemo_mesh=modeltools.nemo.NemoMesh(filemesh,first_j=first_j)

   #ncidmesh=netCDF4.Dataset(filemesh,"r")
   gdept  = nemo_mesh["gdept_0"][0,:]       # Depth of t points
   gdepw  = nemo_mesh["gdepw_0"][0,:]       # Depth of w points
   e3t_ps = nemo_mesh.sliced(nemo_mesh["e3t_ps"] [0,:,:])     # Partial steps of t cell
   e3w_ps = nemo_mesh.sliced(nemo_mesh["e3w_ps"] [0,:,:])     # Partial steps of w cell
   mbathy = nemo_mesh.sliced(nemo_mesh["mbathy"] [0,:,:])     # bathy index
   hdepw  = nemo_mesh.sliced(nemo_mesh["hdepw"]  [0,:,:])     # Total depth of w points
   mbathy = mbathy -1                       # python indexing starts from 0
   nlev   = gdept.size

   mbathy_u,e3u_ps,depthu=nemo_mesh.depth_u_points()
   mbathy_v,e3v_ps,depthv=nemo_mesh.depth_v_points()
   #
   mbathy_u=nemo_mesh.sliced(nemo_mesh.u_to_hycom_u(mbathy_u))
   e3u_ps  =nemo_mesh.sliced(nemo_mesh.u_to_hycom_u(e3u_ps  ))
   depthu  =nemo_mesh.sliced(nemo_mesh.u_to_hycom_u(depthu  ))
   #
   mbathy_v=nemo_mesh.sliced(nemo_mesh.v_to_hycom_v(mbathy_v))
   e3v_ps  =nemo_mesh.sliced(nemo_mesh.v_to_hycom_v(e3v_ps  ))
   depthv  =nemo_mesh.sliced(nemo_mesh.v_to_hycom_v(depthv  ))

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
      s=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
      for k in range(nlev) : # Dont include lowest layer
         s[k,:,:] = nemo_mesh.sliced(ncids.variables["vosaline"][0,k,:,:])
      s = numpy.where(s<1e30,s,0.)
      s = numpy.where(s==ncids.variables["vosaline"]._FillValue,0.,s)
      ncidt=netCDF4.Dataset(filet,"r")
      t=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
      for k in range(nlev) : # Dont include lowest layer
         t[k,:,:] = nemo_mesh.sliced(ncidt.variables["votemper"][0,k,:,:])
      t = numpy.where(t==ncidt.variables["votemper"]._FillValue,0.,t)
      t = numpy.where(t<1e30,t,0.)

      # time from gridT file. 
      time = ncidt.variables["time_counter"][0]
      tunit = ncidt.variables["time_counter"].units
      #print tunit,time
      tmp=cfunits.Units(tunit)
      refy, refm, refd=(1958,1,1)
      tmp2=cfunits.Units("hours since %d-%d-%d 00:00:00"%(refy,refm,refd))            # Units from CF convention
      tmp3=cfunits.Units.conform(time,tmp,tmp2)                                       # Transform to new new unit 
      #print tmp3,type(tmp3)
      tmp3=int(numpy.round(tmp3))
      #print tmp3
      #print datetime.timedelta(hours=tmp3) # Then calculate dt. Phew!
      mydt = datetime.datetime(refy,refm,refd,0,0,0) + datetime.timedelta(hours=tmp3) # Then calculate dt. Phew!
      logger.info("Valid time from gridT file:%s"%str(mydt))


      # Read and calculculate U in hycom U-points. 
      logger.info("gridU  file: %s"%fileu)
      ncidu=netCDF4.Dataset(fileu,"r")
      u=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
      for k in range(nlev) : 
         u[k,:,:] = nemo_mesh.sliced(nemo_mesh.u_to_hycom_u(ncidu.variables["vozocrtx"][0,k,:,:] ))   # Costly, make more efficient if needed
      u = numpy.where(numpy.abs(u)<1e10,u,0.)

      #Calculate barotropic and baroclinic u
      usum=numpy.zeros(u.shape[-2:])
      dsum=numpy.zeros(u.shape[-2:])
      for k in range(u.shape[0]-1) : # Dont include lowest layer
         # TODO: Mid-layer depths seem to be undefined - figure out why ...
         logger.debug("k=%3d, u=%10.3g, mbathy_u[jtest,itest]=%3d,gdepw[k]=%8.2f, depthu[jtest,itest]=%8.2f"%(
            k,u[k,jtest,itest],mbathy_u[jtest,itest],gdepw[k],depthu[jtest,itest]))
         J,I = numpy.where(mbathy_u>k) 
         usum[J,I] = usum[J,I] + u[k,J,I]*dt[k]
         dsum[J,I] = dsum[J,I] + dt[k]
      J,I = numpy.where(mbathy>=0)
      usum[J,I] = usum[J,I] + u[mbathy_u[J,I],J,I]*e3u_ps[J,I]
      dsum[J,I] = dsum[J,I] + e3u_ps[J,I]
      ubaro=numpy.where(dsum>0.1,usum/dsum,0.)

      # Read and calculculate V in hycom V-points. 
      logger.info("gridV  file: %s"%filev)
      ncidv=netCDF4.Dataset(filev,"r")
      v=numpy.zeros((nlev,mbathy.shape[0],mbathy.shape[1]))
      for k in range(nlev) : 
         v[k,:,:] = nemo_mesh.sliced(nemo_mesh.v_to_hycom_v(ncidv.variables["vomecrty"][0,k,:,:] ))   # Costly, make more efficient if needed
      v = numpy.where(numpy.abs(v)<1e10,v,0.)

      #Calculate barotropic and baroclinic v
      vsum=numpy.zeros(v.shape[-2:])
      dsum=numpy.zeros(v.shape[-2:])
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
      #iu = nemo_mesh.periodic_i_shift_right(iu,1)   # u: nemo in cell i is hycom in cell i+1
      #iv = nemo_mesh.arctic_patch_shift_up(iu,1)    # v: nemo in cell j is hycom in cell j+1
      #ip = nemo_mesh.sliced(ip)
      #iu = nemo_mesh.sliced(iu)
      #iv = nemo_mesh.sliced(iv)
      #raise NameError,"test"

      # 2D data
      ncid2d=netCDF4.Dataset(file2d,"r")
      ssh          = nemo_mesh.sliced(ncid2d.variables["sossheig"][0,:,:])
      ssh = numpy.where(ssh==ncid2d.variables["sossheig"]._FillValue,0.,ssh)
      ssh = numpy.where(ssh>1e30,0.,ssh) # Hmmmmm
      #bar_height   = nemo_mesh.sliced(ncid2d.variables["sobarhei"][0,:,:] )
      #dyn_height   = nemo_mesh.sliced(ncid2d.variables["sodynhei"][0,:,:] 
      montg1       = ssh * 9.81  #* 1e-3  # Approx
      logger.warning("TODO:montg pot calculation must be checked...")

      # Write to abfile
      outfile = abfile.ABFileArchv(mydt.strftime(fnametemplate),"w",iexpt=10,iversn=22,yrflag=3,)
      logger.info("Writing 2D variables")
      outfile.write_field(montg1,                ip,"montg1"  ,0,0,1,0)
      outfile.write_field(ssh,                   ip,"srfhgt"  ,0,0,0,0)
      outfile.write_field(numpy.zeros(ssh.shape),ip,"surflx"  ,0,0,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),ip,"salflx"  ,0,0,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),ip,"bl_dpth" ,0,0,0,0) # Not used
      outfile.write_field(numpy.zeros(ssh.shape),ip,"mix_dpth",0,0,0,0) # Not used
      outfile.write_field(ubaro                 ,iu,"u_btrop" ,0,0,0,0) # u: nemo in cell i is hycom in cell i+1
      outfile.write_field(vbaro                 ,iv,"v_btrop" ,0,0,0,0) # v: nemo in cell j is hycom in cell j+1
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
         outfile.write_field(ul      ,iu,"u-vel.",0,0,k+1,0) # u: nemo in cell i is hycom in cell i+1
         outfile.write_field(vl      ,iv,"v-vel.",0,0,k+1,0) # v: nemo in cell j is hycom in cell j+1
         outfile.write_field(dtl*onem,ip,"thknss",0,0,k+1,0)
         outfile.write_field(sl      ,ip,"salin" ,0,0,k+1,0)
         outfile.write_field(tl      ,ip,"temp"  ,0,0,k+1,0)


      # TODO: Process ice data
      ncid2d.close()
      ncids.close()
      ncidt.close()
      ncidu.close()
      ncidv.close()
      outfile.close()

      logger.info("Finished writing %s.[ab] "% mydt.strftime(fnametemplate))
   nemo_mesh = []





if __name__ == "__main__" :

   parser = argparse.ArgumentParser(
         description='This tool will convert NEMO netcdf files to hycom archive files. It will also create grid and topo files for hycom.'
#         epilog="""Example:
#   %s GLORYS2V3_mesh_mask.nc GLORYS2V4_1dAV_20130101_20130102_grid2D_R20130102.nc 
#   
#   %s GLORYS2V3_mesh_mask.nc --first_j=100  GLORYS2V4_1dAV_20130101_20130102_grid2D_R20130102.nc 
#   """ %(os.path.basename(sys.argv[0]),os.path.basename(sys.argv[0]))
         )
   parser.add_argument('--first_j',   type=int,default=0,help="first j-index to process. Defaults to 0")
   parser.add_argument('--mean',   action="store_true",default=False,help="if mean flag is set, a mean archive will be created")
   parser.add_argument('meshfile',   type=str,help="NEMO mesh file in netcdf format")
   parser.add_argument('grid2dfile', type=str, nargs="+",help="NEMO 2D data file in netcdf format")

   args = parser.parse_args()
   main(args.meshfile,args.grid2dfile,first_j=args.first_j,mean_file=args.mean)


