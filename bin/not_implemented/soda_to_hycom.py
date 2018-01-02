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
import nemo_mesh_to_hycom

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


def soda_to_regional_grid(fid) :
   # Assume this is for t-cell centers
   plon = fid["longitude"][:]
   plat = fid["latitude"][:]
   dlon = plon[1]-plon[0]
   dlat = plat[1]-plat[0]
   ulon = plon - dlon
   ulat = numpy.copy(plat)
   vlon = numpy.copy(plon )
   vlat = plat-dlat
   qlon = plon-dlon
   qlat = plat-dlat

   plon,plat=numpy.meshgrid(plon,plat)
   ulon,ulat=numpy.meshgrid(ulon,ulat)
   vlon,vlat=numpy.meshgrid(vlon,vlat)
   qlon,qlat=numpy.meshgrid(qlon,qlat)

   # Rough estimate of grid sizes
   scpx=dlon*(110574.2727)*numpy.cos(numpy.radians(plat))
   scpy=dlat*(110574.2727)*numpy.ones(plat.shape)
   scux=dlon*(110574.2727)*numpy.cos(numpy.radians(ulat))
   scuy=dlat*(110574.2727)*numpy.ones(plat.shape)
   scvx=dlon*(110574.2727)*numpy.cos(numpy.radians(vlat))
   scvy=dlat*(110574.2727)*numpy.ones(plat.shape)
   scqx=dlon*(110574.2727)*numpy.cos(numpy.radians(qlat))
   scqy=dlat*(110574.2727)*numpy.ones(plat.shape)
   corio = numpy.sin(numpy.radians(qlat)) * 4. * numpy.pi / 86164.0 
   asp = numpy.where(scpy==0.,99.0,scpx/scpy)
   pang=numpy.zeros(plat.shape)

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





   

def main(startdate,enddate,first_j=0) :

   
   soda_template="/work/shared/nersc/msc/SODA/3.3.1/monthly/soda3.3.1_mn_ocean_reg_%Y.nc"

   # open blkdat.input. Get nesting frequency
   bp=modeltools.hycom.BlkdatParser("blkdat.input")
   nestfq=bp["nestfq"]
   bnstfq=bp["bnstfq"]

   # Read soda-grid and topo from first file
   fid = netCDF4.Dataset(startdate.strftime(soda_template),"r")
   depth=fid["depth"][:]
   soda_to_regional_grid(fid) 


   # Get bathymetry. Set to 0 wherever salinit in top layer is undefined.
   bathy=numpy.zeros(fid["salt"].shape[-2:])
   for k in range(depth.size) :
      bathy[~numpy.squeeze(fid["salt"][0,k,:,:].mask)] = depth[k]
   abfile.write_bathymetry("SODA",22,bathy,0.)
   fid.close()
   
   # TODO:
   ip = bathy == 0.
   iu = bathy == 0.
   iv = bathy == 0.

   onem=9806

   #if mean_file :
   #   fnametemplate_out="archm.%Y_%j_%H"
   #else :
   hycom_template="archv.%Y_%j_%H"



   # Loop over nestfq, bnstfq.
   deltat=enddate-startdate
   dsec=deltat.days*86400+deltat.seconds
   baroclinic_nest_times=[startdate+datetime.timedelta(seconds=s) for s in numpy.arange(0,dsec,nestfq*86400)]
   barotropic_nest_times=[startdate+datetime.timedelta(seconds=s) for s in numpy.arange(0,dsec,bnstfq*86400)]
   tmp=sorted(set(barotropic_nest_times+baroclinic_nest_times))
   for dt in tmp :

      logger.info("Processing time %s"%str(dt))

      # Get "mid-month" dates
      if dt.day >= 15 :
         nm=1+dt.month%12
         ny=dt.year + dt.month/12
         mm0 = datetime.datetime(dt.year,dt.month,15,0,0,0)
         mm1 = datetime.datetime(ny,nm,15,0,0,0) 
      else :
         lm=1+(12 + dt.month-2)%12
         ly=dt.year-lm/12
         print dt.month,lm,ly
         mm0 = datetime.datetime(ly,lm,15,0,0,0) 
         mm1 = datetime.datetime(dt.year,dt.month,15,0,0,0)
      
      # Linear interpolation weights
      deltat=mm1-mm0
      deltat=deltat.days+deltat.seconds/86400.
      w1    = dt - mm0
      w1    = w1.days+w1.seconds/86400.
      w1    = w1/deltat
      w0    = 1.-w1
      flnm0=mm0.strftime(soda_template)
      flnm1=mm1.strftime(soda_template)
      logger.info("Time %s, file %s at %s(w=%.4f) , file %s at %s(w=%.4f)"%(str(dt),flnm0,str(mm0),w0,flnm1,str(mm1),w1))

      
      # Open files
      # TODO: reuse pointers/fields
      fid0=netCDF4.Dataset(flnm0,"r")
      fid1=netCDF4.Dataset(flnm1,"r")

      # Calculate temperature, velocity
      temp= w0*fid0["temp"][0,:,:,:]  + w1*fid1["temp"][0,:,:,:]
      salt= w0*fid0["salt"][0,:,:,:]  + w1*fid1["salt"][0,:,:,:]
      utot= w0*fid0["u"]   [0,:,:,:]  + w1*fid1["u"]   [0,:,:,:]
      vtot= w0*fid0["v"]   [0,:,:,:]  + w1*fid1["v"]   [0,:,:,:]


      #NB: No checks for missing values yet !
      ubaro = numpy.sum(utot,0)
      vbaro = numpy.sum(vtot,0)
      u     = utot - ubaro
      v     = vtot - vbaro

      # 2D vars
      anompb = w0*fid0["anompb"         ][0,:,:]  + w1*fid1["anompb"         ][0,:,:]
      ssh    = w0*fid0["ssh"            ][0,:,:]  + w1*fid1["ssh"            ][0,:,:]
      salflx = w0*fid0["salt_flux_total"][0,:,:]  + w1*fid1["salt_flux_total"][0,:,:]
      surflx = w0*fid0["net_heating"    ][0,:,:]  + w1*fid1["salt_flux_total"][0,:,:]
      montg1 = numpy.zeros(ssh.shape)


      # Write to abfile
      outfile = abfile.ABFileArchv(dt.strftime(hycom_template),"w",iexpt=10,iversn=22,yrflag=3)
      logger.info("Writing 2D variables")
      outfile.write_field(montg1,                ip,"montg1"  ,0,0,1,0)
      outfile.write_field(ssh,                   ip,"srfhgt"  ,0,0,0,0)
      outfile.write_field(surflx                ,ip,"surflx"  ,0,0,0,0) 
      outfile.write_field(salflx                ,ip,"salflx"  ,0,0,0,0) 
      outfile.write_field(numpy.zeros(ssh.shape),ip,"bl_dpth" ,0,0,0,0) 
      outfile.write_field(numpy.zeros(ssh.shape),ip,"mix_dpth",0,0,0,0) 
      outfile.write_field(ubaro                 ,iu,"u_btrop" ,0,0,0,0) 
      outfile.write_field(vbaro                 ,iv,"v_btrop" ,0,0,0,0) 
      #outfile.close() ; raise NameError,"test"
      for k in numpy.arange(u.shape[0]) :
         if k%10==0 : logger.info("Writing 3D variables, level %d of %d"%(k+1,u.shape[0]))

         if k==0 :
            dtl = depth[0]*numpy.where(bathy>=depth[k],1,0)
         else :
            dtl = (depth[k]-depth[k-1])*numpy.where(bathy>=depth[k],1,0)
         print dtl.min(),dtl.max()

         templ = temp[k,:,:]
         saltl = salt[k,:,:]

         # Set to layer above if undefined
         templ[dtl<=0.] = 20.
         saltl[dtl<=0.] = 35.

         outfile.write_field(u[k,:,:]   ,iu,"u-vel.",0,0,k+1,0) 
         outfile.write_field(v[k,:,:]   ,iv,"v-vel.",0,0,k+1,0) 
         outfile.write_field(dtl*onem   ,ip,"thknss",0,0,k+1,0)
         outfile.write_field(saltl      ,ip,"salin" ,0,0,k+1,0)
         outfile.write_field(templ      ,ip,"temp"  ,0,0,k+1,0)

         oldsaltl=saltl
         oldtempl=templ

      # TODO: reuse pointers/fields
      outfile.close()
      fid0.close()
      fid1.close()
      raise NameError,"check vals"




















   raise NameError,"test"

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
      tmp=cfunits.Units(tunit)
      refy, refm, refd=(1958,1,1)
      tmp2=cfunits.Units("seconds since %d-%d-%d 00:00:00"%(refy,refm,refd))            # Units from CF convention
      tmp3=cfunits.Units.conform(time,tmp,tmp2)                                         # Transform to new new unit 
      mydt = datetime.datetime(refy,refm,refd,0,0,0) + datetime.timedelta(seconds=tmp3) # Then calculate dt. Phew!


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
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='From blkdat file and start/stop times, create hycom archv files from soda data ')
   parser.add_argument('startdate',  action=DateTimeParseAction)
   parser.add_argument('enddate',    action=DateTimeParseAction)
   parser.add_argument('--first_j',   type=int,default=0)
   args = parser.parse_args()
   main(args.startdate,args.enddate,first_j=args.first_j)


