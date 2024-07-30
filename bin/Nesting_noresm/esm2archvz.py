#!/usr/bin/env python
#import modeltools.nemo
print("loading 0")
import argparse
import abfile.abfile as abf
import logging
import re
import os.path
import glob
from netCDF4 import Dataset, num2date,date2num
import numpy as np

# Hycom-ified NORCPM files. Approach:
# 1) Create hycom archv files and topo/region files using this routine
# 2) Interpolate NORESM topo file to target region
# 3) Merge target region topography and new region. Set up experiment to use new topo
# 4) Interpolate archive file to new region / experiment
# 5) Remap archive files vertically
# 6) Add montgomery potential in 1st layer to file
# 
# History:
# Mostafa Bakhoday-Paskyabi, December 2018
# Mostafa Bakhoday-Paskyabi, March 2019, bio-nesting (solving interpolation issue)

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



idm=360
jdm=385
latmin=-81
latmax=90
lonmin=0
lonmax=360
mask_method=2
spval=2.0**100.
timeavg_method   = 0             # (1) time average of two consecutive netcdf files; and other values no temporal averaging.


def maplev(a,lpp=1):
    # gapfilling method
    jm,im=a.shape
    J,I=np.where(~np.isnan(a))
    with np.errstate(invalid='ignore'):
        av=np.nansum(a[J,I])/(len(I)*len(J))
    J,I=np.where(np.isnan(a))
    a[J,I]=av
    b=a
    #lpp=1  # it is better to be set to 100, but for practial reasnon, we keep it very small for now
    i=list(range(1,im-1))
    j=list(range(1,jm-1))
    ip1=list(range(2,im))
    jp1=list(range(2,jm))
    im1=list(range(0,im-2))
    jm1=list(range(0,jm-2))
    
    cc=np.zeros(a.shape)
    for k in range(lpp):
        cc[1:-2,1:-2]=b[1:-2,1:-2]+.5/4*( b[1:-2,2:-1]+b[0:-3,1:-2]+b[1:-2,0:-3]+b[2:-1,1:-2]-4.*b[1:-2,1:-2] )
        cc[:,0]=cc[:,1]
        cc[:,im-1]=cc[:,im-2]
        cc[jm-1,:]=cc[jm-2,:]
        b[J,I]=cc[J,I]

    a[J,I]=cc[J,I]
    J,I=np.where(np.isnan(a))
    a[J,I]=0.
    return a


def u_to_hycom_u(field2d)  :
   return np.roll(field2d,1,axis=1)


# nemo V at i,j -> HYCOM U at i,j+1. 
def v_to_hycom_v(field2d,extrapolate="none")  :
   myfield=np.copy(field2d)
   # For now the bottom row is just replicated
   myfield[1:,:] = myfield[:-1,:]
   
   return myfield

def periodic_i_shift_right(field,istep) :
   # shift field left by istep steps
   field2  = np.roll(field,istep,axis=1)
   return field2

def arctic_patch_shift_down(field,jstep) :
   # shift field down
   if jstep != 1 :
      raise NameError("Arctic_patch_shift only with jstep=1 for now")
   field2 = np.copy(field)
   field2[0:-1,:] = field2[1:,:] # Shift down
   tmp=field2[-1,:]              # Top row as top ...
   field2[-1,:] = tmp[::-1]      # .. but reversed direction
   return field2

def depth_u_points(depth) :
   depthip1  = periodic_i_shift_right(depth ,-1)    # noresm values at cell i+1
   depthu  =np.mean(np.array([depth ,depthip1]),axis=0)
   return depthu
      
def depth_v_points(depth) :
   depthjp1  = arctic_patch_shift_down(depth ,1)    # nemo values at cell j+1
   depthv  =np.mean(np.array([depth ,depthjp1]),axis=0)
   return depthv


def read_mesh(filemesh):
   # Read longitude, latitude and depth form the original NORESM-files
   meshfile_area=filemesh
   filepath=os.path.dirname(filemesh)
   filename=os.path.basename(filemesh)
   meshfile_depth=filepath+"/deptho"+filename[9:]

   ncida=Dataset(meshfile_area,"r")
   ncidd=Dataset(meshfile_depth,"r")

   depth = ncidd.variables["deptho" ][:,:] 
   plon = ncida.variables["longitude"][:,:]     # bathy index
   plat  = ncida.variables["latitude"][:,:]     # Total depth of w points

   ncida.close()
   ncidd.close()
  
   return depth,plon,plat

def search_biofile(bio_path,dt):

   logger.info("BIO")
   # filename format for NORESM: var_id+"_Omon_NorESM2-MM_historical_r1i1p1f1_gr_195408_extrap.nc"
   lst=glob.glob(bio_path+"no3_"+esm_id_+"_gr_195408_extrap%s_extrap.nc"%str(dt[:-2]))
   df = np.zeros(len(lst))*np.nan 
   val, idx = min((val, idx) for (idx, val) in enumerate(np.abs(df)))
   return idx,lst[idx]

def check_inputs(x, y, Z, points, mode, bounds_error):
    """Check inputs for interpolate2d function
    """

    msg = 'Only mode "linear" and "constant" are implemented. I got %s' % mode
    if mode not in ['linear', 'constant']:
        raise RuntimeError(msg)

    try:
        x = np.array(x)
    except Exception as e:
        msg = ('Input vector x could not be converted to np array: '
               '%s' % str(e))
        raise Exception(msg)

    try:
        y = np.array(y)
    except Exception as e:
        msg = ('Input vector y could not be converted to np array: '
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
        Z = np.array(Z)
        m, n = Z.shape
    except Exception as e:
        msg = 'Z must be a 2D np array: %s' % str(e)
        raise Exception(msg)

    Nx = len(x)
    Ny = len(y)
    msg = ('Input array Z must have dimensions %i x %i corresponding to the '
           'lengths of the input coordinates x and y. However, '
           'Z has dimensions %i x %i.' % (Nx, Ny, m, n))
    if not (Nx == m and Ny == n):
        raise RuntimeError(msg)

    # Get interpolation points
    points = np.array(points)
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


def main(filemesh,grid2dfiles,first_j=0,mean_file=True,iexpt=10,iversn=22,yrflag=3,bio_path=None,esm_id=None) :

   if mean_file :
      fnametemplate="archm.%Y_%j_%H"
   else :
      fnametemplate="archv.%Y_%j_%H"
   itest=1
   jtest=200
   depth,plon,plat=read_mesh(filemesh)

   depthu=depth_u_points(depth)
   depthv=depth_v_points(depth)

   # Loop over input files. All must be in same directory
   for file2d in grid2dfiles : 

      # See if actually a grid2D file
      dirname=os.path.dirname(file2d)
      m=re.match("(thetao)(_.*\.nc)",os.path.basename(file2d))
      if not m :
         msg="File %s is not a grid2D file, aborting"%file2d
         logger.error(msg)  
         raise ValueError(msg)

      filepath=os.path.dirname(file2d)
      filename=os.path.basename(file2d)
      print(filepath)
      print(filename)
      # Construct remaining files
      filet  =file2d
      files  =filepath + "/so" + filename[6:]
      fileu  =filepath + "/uo" + filename[6:]
      filev  =filepath + "/vo" + filename[6:]
      if bio_path:
         file_no3 = bio_path + "no3"  + filename[6:]
         file_po4 = bio_path + "po4"  + filename[6:]
         file_si  = bio_path + "si"  + filename[6:]
         file_o2  = bio_path + "o2"  + filename[6:]
 
      lenstr=len(filename); bsubstr=lenstr-18; esubstr=lenstr-17;
      print(lenstr,bsubstr,esubstr)
      filessh  =filepath + "/zos" + filename[6:bsubstr]+"n"+filename[esubstr:]
      logger.info("grid2D file: %s"%filessh)

      # P-points
      logger.info("gridS  file: %s"%files)
      logger.info("gridT  file: %s"%filet)
      ncids=Dataset(files,"r")
      ncidt=Dataset(filet,"r")
      ncidu=Dataset(fileu,"r")
      ncidv=Dataset(filev,"r")
      if bio_path:
         ncidno3=Dataset(file_no3,"r")
         ncidsi =Dataset(file_si,"r")
         ncidpo4=Dataset(file_po4,"r")
         ncido2 =Dataset(file_o2,"r")

      # time from gridT file. 
      time = ncidt.variables["time"][0]
      tunit = ncidt.variables["time"].units
      t_cal=ncidt.variables["time"].calendar
      date=num2date(time,units = tunit,calendar = t_cal)

      # Read and calculculate U in hycom U-points. 
      logger.info("gridU, gridV, gridT & gridS  file")
      u=np.squeeze(ncidu.variables["uo"][:,:,:,:])
      v=np.squeeze(ncidv.variables["vo"][:,:,:,:])
      t=np.squeeze(ncidt.variables["thetao"][:,:,:,:])
#      t_fill=ncidt.variables["thetao"]._FillValue
      s=np.squeeze(ncids.variables["so"][:,:,:,:])
#      s_fill=ncids.variables["so"]._FillValue

      if bio_path:
          no3=np.squeeze(ncidno3.variables["no3"][:,:,:,:])
          no3=no3 * 6.625 * 12.01 * 1000.0 # convert from mol/m3 to mg C/m3
          si=np.squeeze(ncidsi.variables["si"][:,:,:,:])
          si=si * 6.625 * 12.01 * 1000.0 # convert from mol/m3 to mg C/m3 
          po4=np.squeeze(ncidpo4.variables["po4"][:,:,:,:])
          po4=po4 * 106.0 * 12.01 * 1000.0 # convert from mol/m3 to mg C/m3 
          o2=np.squeeze(ncido2.variables["o2"][:,:,:,:])
          o2=o2 * 1000.0 # convert from mol/m3 to mmol/m3 

      lev_bnds=ncidu.variables["lev_bnds"][:,:]
      lev=ncidu.variables["lev"][:]

      nlev=np.size(lev)
      dz=lev_bnds[:,1]-lev_bnds[:,0]
   #   uu=np.zeros(np.shape(u))
   #   vv=np.zeros(np.shape(v))
   #   for k in range(nlev) : 
   #      uu[k,:,:] = np.squeeze(u_to_hycom_u(u[k,:,:] ))
   #      vv[k,:,:] = np.squeeze(v_to_hycom_v(v[k,:,:] ))

      u = np.where(np.abs(u)<1e10,u,0.)
      v = np.where(np.abs(v)<1e10,v,0.)

      logger.info("Calculate barotropic velocities ...")

      #U and V are aligned N-S/E-W so need to rotate to be along the grid
      uu_x=np.zeros(np.shape(u)) # U along x-grid
      vv_y=np.zeros(np.shape(v)) # V along y-grif
      for k in range(nlev) :
         uu_x[k,:,:] = np.squeeze(u_to_hycom_u(u[k,:,:] ))  
         vv_y[k,:,:] = np.squeeze(u_to_hycom_u(v[k,:,:] ))  

      del u,v

      #Calculate barotropic and baroclinic velocities
      mbathy=np.zeros(np.shape(depth))
      ii,jj=np.shape(depth)
      mbathy[np.where(depth>=1.0e15)]=0.0
      mbathy[np.where((depth<1.0e15) & (depth >= np.max(lev_bnds[:,0])))]=np.size(lev)-1
      for i in range(ii):
          for j in range(jj):
              if depth[i,j] < np.max(lev_bnds[:,0]):
                  mbathy[i,j]=np.min(np.where(lev_bnds[:,0]>depth[i,j]))-1

      mbathy_u=np.zeros(np.shape(depthu))
      mbathy_u[np.where(depthu>=1.0e15)]=0.0
      mbathy_u[np.where((depthu<1.0e15) & (depthu >= np.max(lev_bnds[:,0])))]=np.size(lev)-1
      for i in range(ii):
          for j in range(jj):
              if depthu[i,j] < np.max(lev_bnds[:,0]):
                  mbathy_u[i,j]=np.min(np.where(lev_bnds[:,0]>depthu[i,j]))-1

      mbathy_v=np.zeros(np.shape(depthv))
      mbathy_v[np.where(depthv>=1.0e15)]=0.0
      mbathy_v[np.where((depthv<1.0e15) & (depthv >= np.max(lev_bnds[:,0])))]=np.size(lev)-1
      for i in range(ii):
          for j in range(jj):
              if depthv[i,j] < np.max(lev_bnds[:,0]):
                  mbathy_v[i,j]=np.min(np.where(lev_bnds[:,0]>depthv[i,j]))-1

      mbathy_u=mbathy_u.astype(int)      
      mbathy_v=mbathy_v.astype(int)

      print("mbathy_u",mbathy_u[200,200],depthu[200,200])#,lev[mbathy_u[200,200]])
      print("mbathy_v",mbathy_v[200,200],np.max(mbathy_v),depthv[200,200])#,lev[mbathy_v[200,200]])

      usum=np.zeros(uu_x.shape[-2:])
      dsumu=np.zeros(uu_x.shape[-2:])
      vsum=np.zeros(vv_y.shape[-2:])
      dsumv=np.zeros(vv_y.shape[-2:])

      for k in range(uu_x.shape[0]-1) : # Dont include lowest layer
         J,I = np.where(mbathy_u>k) 
         usum[J,I] = usum[J,I] + uu_x[k,J,I]*dz[k]
         dsumu[J,I] = dsumu[J,I] + dz[k]
         J,I = np.where(mbathy_v>k)
         vsum[J,I] = vsum[J,I] + vv_y[k,J,I]*dz[k]
         dsumv[J,I] = dsumv[J,I] + dz[k]

      print(dsumu[200,200],depthu[200,200],np.shape(uu_x))
      J,I = np.where(mbathy_u>=0)
      usum[J,I] = usum[J,I] + uu_x[mbathy_u[J,I],J,I]*(depthu[J,I]-dsumu[J,I])
      print("Bottom layer",depthu[200,200]-dsumu[200,200])
      dsumu[J,I] = dsumu[J,I] + depthu[J,I]-dsumu[J,I]
      print(dsumu[200,200],depthu[200,200])
      dsumu=np.where(abs(dsumu)<1e-2,0.05,dsumu)
      ubaro=np.where(dsumu>0.1,usum/dsumu,0.)
      J,I = np.where(mbathy_v>=0)
      vsum[J,I] = vsum[J,I] + vv_y[mbathy_v[J,I],J,I]*(depthv[J,I]-dsumv[J,I])
      dsumv[J,I] = dsumv[J,I] + depthv[J,I]-dsumv[J,I]
      dsumv=np.where(abs(dsumv)<1e-2,0.05,dsumv)
      vbaro=np.where(dsumv>.1,vsum/dsumv,0.)
      fnametemplate="archv.%Y_%j"
      oname=date.strftime(fnametemplate)+"_00"

      # model day
      logger.info("Model day in HYCOM is %s"%str(date.dayofyr))
      print(date,date.day,dir(date))
      model_day=time
      print(model_day)

      # Masks (land:True)
      if mask_method == 1 :
         ip = mbathy   == -1
         iu = mbathy_u == -1
         iv = mbathy_v == -1
      else :
         ip = depth   == 0
         iu = depthu  == 0
         iv = depthv  == 0

      # 2D data
      ncid2d=Dataset(filessh,"r")
      ssh = ncid2d.variables["zos"][0,:,:]
      ssh = np.where(ssh==ncid2d.variables["zos"]._FillValue,0.,ssh)
      ssh = np.where(ssh>1e10,0.,ssh*9.81) # NB: HYCOM srfhgt is in geopotential ...
      montg1=np.zeros(ssh.shape)

      # Write to abfile
      outfile = abf.ABFileArchv("./data/"+oname,"w",iexpt=iexpt,iversn=iversn,yrflag=yrflag,)

      logger.info("Writing 2D variables")
      outfile.write_field(montg1,                ip,"montg1"  ,0,model_day,1,0)
      outfile.write_field(ssh,                   ip,"srfhgt"  ,0,model_day,0,0)
      outfile.write_field(np.zeros(ssh.shape),ip,"surflx"  ,0,model_day,0,0) # Not used
      outfile.write_field(np.zeros(ssh.shape),ip,"salflx"  ,0,model_day,0,0) # Not used
      outfile.write_field(np.zeros(ssh.shape),ip,"bl_dpth" ,0,model_day,0,0) # Not used
      outfile.write_field(np.zeros(ssh.shape),ip,"mix_dpth",0,model_day,0,0) # Not used
      outfile.write_field(ubaro                 ,iu,"u_btrop" ,0,model_day,0,0) # u: nemo in cell i is hycom in cell i+1
      outfile.write_field(vbaro                 ,iv,"v_btrop" ,0,model_day,0,0) # v: nemo in cell j is hycom in cell j+1
      ny=mbathy.shape[0];nx=mbathy.shape[1]
      error=np.zeros((ny,nx))
      for k in np.arange(uu_x.shape[0]) :
         if bio_path:
            no3l=np.squeeze(no3[k,:,:])
            no3l=maplev(no3l)
            no3l = np.where(no3l<1e8,no3l,np.nan)
            no3l = np.minimum(np.maximum(maplev(no3l),0),1.0e8)
            po4l=np.squeeze(po4[k,:,:])
            po4l=maplev(po4l)
            po4l = np.where(po4l<1e8,po4l,np.nan)
            po4l = np.minimum(np.maximum(maplev(po4l),0),1.0e8)
            sil=np.squeeze(si[k,:,:])
            sil=maplev(sil)
            sil = np.where(sil<1e8,sil,np.nan)
            sil = np.minimum(np.maximum(maplev(sil),0),1.0e8)
            o2l=np.squeeze(o2[k,:,:])
            o2l=maplev(o2l)
            o2l = np.where(o2l<1e8,o2l,np.nan)
            o2l = np.minimum(np.maximum(maplev(o2l),0),1.0e8)
            if k%10==0 : logger.info("Writing 3D variables including BIO, level %d of %d"%(k+1,uu_x.shape[0]))
         else:
            if k%10==0 : logger.info("Writing 3D variables, level %d of %d"%(k+1,uu_x.shape[0]))
         #
         ul = np.squeeze(uu_x[k,:,:]) - ubaro # Baroclinic velocity
         vl = np.squeeze(vv_y[k,:,:]) - vbaro # Baroclinic velocity

         # Layer thickness
         dzl=np.zeros(ul.shape)
         if k < uu_x.shape[0]-1 :
            J,I = np.where(mbathy>k)
            dzl[J,I] = dz[k]
            J,I = np.where(mbathy==k)
            dzl[J,I] = depth[J,I]-lev_bnds[k,0]
         else:
            J,I = np.where(mbathy==k)
            dzl[J,I] = depth[J,I]-lev_bnds[k,0]
         
         sl = np.squeeze(s[k,:,:])
         tl = np.squeeze(t[k,:,:])
         sl = np.where(sl<1e2,sl,np.nan)
         sl = np.minimum(np.maximum(maplev(sl),25),80.)
         tl = np.where(tl<=5e2,tl,np.nan)
         tl = np.minimum(np.maximum(maplev(tl),-5.),50.)

         # Fill empty layers with values from above
         if k > 0 :
            K= np.where(dzl < 1e-4)

            tl[K] = tl_above[K]
            sl[K] = sl_above[K]
            if bio_path:
               no3l[K] = no3_above[K]
               po4l[K] = po4_above[K]
               sil[K] = si_above[K]
               o2l[K] = o2_above[K]
   

         onem=9806.
         outfile.write_field(ul      ,iu,"u-vel.",0,model_day,k+1,0) # u: nemo in cell i is hycom in cell i+1
         outfile.write_field(vl      ,iv,"v-vel.",0,model_day,k+1,0) # v: nemo in cell j is hycom in cell j+1
         outfile.write_field(dzl*onem,ip,"thknss",0,model_day,k+1,0)
         outfile.write_field(tl      ,ip,"temp"  ,0,model_day,k+1,0)
         outfile.write_field(sl      ,ip,"salin" ,0,model_day,k+1,0)
         if bio_path :
            outfile.write_field(no3l      ,ip,"ECO_no3" ,0,model_day,k+1,0)
            outfile.write_field(po4l      ,ip,"ECO_pho" ,0,model_day,k+1,0)
            outfile.write_field(sil       ,ip,"ECO_sil" ,0,model_day,k+1,0)
            outfile.write_field(o2l       ,ip,"ECO_oxy" ,0,model_day,k+1,0)

         tl_above=np.copy(tl)
         sl_above=np.copy(sl)
         if bio_path:
            no3_above=np.copy(no3l)
            po4_above=np.copy(po4l)
            si_above=np.copy(sil)
            o2_above=np.copy(o2l)
         

      # TODO: Process ice data
      ncid2d.close()
      outfile.close()
      ncidt.close()
      ncids.close()
      ncidu.close()
      ncidv.close()
      if bio_path:
         ncidno3.close()
         ncidsi.close()
         ncidpo4.close()
         ncido2.close()
   nemo_mesh = []





if __name__ == "__main__" :

   parser = argparse.ArgumentParser(
         description='This tool will convert NEMO netcdf files to hycom archive files. It will also create grid and topo files for hycom.'
         )
   parser.add_argument('--first-j',   type=int,default=0,help="first j-index to process. Defaults to 0")
   parser.add_argument('--mean',   action="store_true",default=False,help="if mean flag is set, a mean archive will be created")
   parser.add_argument('meshfile',   type=str,help="ESM mesh file in netcdf format")
   parser.add_argument('grid2dfile', type=str, nargs="+",help="ESM 2D data file in netcdf format")
   parser.add_argument('--iexpt',    type=int,default=10,  help="    ")
   parser.add_argument('--iversn',   type=int,default=22,  help="    ")
   parser.add_argument('--yrflag',   type=int,default=3,   help="    ")
   parser.add_argument('--bio_path',   type=str,   help="    ")
   parser.add_argument('--interp_method',   type=int,default=3,   help="    ")
   parser.add_argument('--esm_id',   type=str,default="thetao_Omon_NorESM2-MM_historical_r1i1p1f1",   help="    ")

   args = parser.parse_args()
   main(args.meshfile,args.grid2dfile,first_j=args.first_j,mean_file=args.mean,iexpt=args.iexpt,iversn=args.iversn,yrflag=args.yrflag,bio_path=args.bio_path,esm_id=args.esm_id)
