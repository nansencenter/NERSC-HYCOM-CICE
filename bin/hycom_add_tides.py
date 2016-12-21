#!/usr/bin/env python
import modeltools.hycom
import modeltools.tools
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
import netCDF4
import logging
import re
import cfunits
import os
import os.path

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


def check_grids(plon,plon2,plat,plat2) :
   # Check grids match
   maxdlon = numpy.amax(numpy.abs(plon -plon2 ))
   maxdlat = numpy.amax(numpy.abs(plat -plat2 ))
   if maxdlon > 1e-4 or maxdlat > 1e-4 :
      msg="grid file mismatch max lon diff =%g , max lat diff = %g"%(maxdlon,maxdlat)
      logger.error(msg)
      raise ValueError,msg

def check_depths(depth,depth2):
   # Check depths match. NB: Since central region can be filled, we only check
   # where depth > 0
   tmp1=depth>.1
   tmp2=depth2>.2
   tmp1=numpy.logical_and(tmp1,tmp2)
   maxddep = numpy.amax(numpy.abs(depth-depth2)[tmp1])
   if maxddep > 1e-4 :
      msg="depth file mismatch max diff =%g , max lat diff = %g"%(maxddep)
      logger.error(msg)
      raise ValueError,msg

def cf_time_to_datetime(times,time_units) :
   # Time processing
   tmp=cfunits.Units(time_units)
   refy, refm, refd=(1950,1,1)                                              # Reference time for this routine
   tmp2=cfunits.Units("seconds since %d-%d-%d 00:00:00"%(refy,refm,refd))   # Units from CF convention
   tmp3=cfunits.Units.conform(times,tmp,tmp2)                               # Transform to new new unit  (known to this routine)
   # Then calculate dt. Phew!
   mydt = [ datetime.datetime(refy,refm,refd,0,0,0) +
         datetime.timedelta(seconds=int(elem)) for elem in tmp3]
   return mydt

def diff_in_seconds(deltat) :
   return deltat.days*86400. + deltat.seconds


def main(tide_file,archv_files,include_uv=False):

   # 1) If this routine is called without any archive files (empty list), then 
   # Files suitable for barotropic nesting only are created. The new archive files are then 
   # chosen to match times in tide file.

   # 2) If routines are called with archive files, then times matching the archive file times are
   # sought from the tide file. It they are found, srfhgt and montg1 are adjusted 
   # to match the new tidal data.


   # Read plon,plat and depth from regional files. Mainly used to check that
   # grid is ok ...
   logger.info("Opening regional.grid.[ab]")
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   pang=gfile.read_field("pang") # For rotation of tidal current
   gfile.close()

   logger.info("Opening regional.depth.[ab]")
   bathyfile=abfile.ABFileBathy("regional.depth","r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   depth=bathyfile.read_field("depth")
   bathyfile.close()
   depth = depth.filled(0.)
   ip=depth>0.0
   iu=numpy.copy(ip)
   iu[:,1:] = numpy.logical_and(iu[:,1:],iu[:,0:-1])
   iv=numpy.copy(ip)
   iv[1:,:] = numpy.logical_and(iv[1:,:],iv[0:-1,:])

   # Open netcdf file, get time variable and some basic stuff
   print os.getcwd(),tide_file
   logger.info("Opening %s"%tide_file)
   nc_h = netCDF4.Dataset(tide_file,"r")
   plon_h=nc_h.variables["longitude"][:]
   plat_h=nc_h.variables["latitude"][:]
   depth_h=nc_h.variables["depth"][:]
   check_grids(plon,plon_h,plat,plat_h)
   check_depths(depth,depth_h)

   # Time processing for tidal elevations 
   time_h=nc_h.variables["time"][:]
   tunit = nc_h.variables["time"].units
   mydt_h = cf_time_to_datetime(time_h,tunit) 

   if include_uv :

      m=re.match("^(.*)_h.nc$",tide_file)
      if m :
         tide_file_u = m.group(1)+"_u.nc"
      else :
         msg="Unable to guesstimate tidal u component from tidsl heights file %s "%tide_file_h
         logger.error(msg)
         raise ValueError,msg

      m=re.match("^(.*)_h.nc$",tide_file)
      if m :
         tide_file_v = m.group(1)+"_v.nc"
      else :
         msg="Unable to guesstimate tidal u component from tidsl heights file %s "%tide_file_h
         logger.error(msg)
         raise ValueError,msg

      logger.info("Opening %s"%tide_file_u)
      nc_u = netCDF4.Dataset(tide_file_u,"r")
      plon_u=nc_u.variables["longitude"][:]
      plat_u=nc_u.variables["latitude"][:]
      depth_u=nc_u.variables["depth"][:]
      check_grids(plon,plon_u,plat,plat_u)
      check_depths(depth,depth_u)

      # Time processing for tidal elevations 
      time_u=nc_u.variables["time"][:]
      tunit = nc_u.variables["time"].units
      mydt_u = cf_time_to_datetime(time_u,tunit) 

      logger.info("Opening %s"%tide_file_v)
      nc_v = netCDF4.Dataset(tide_file_v,"r")
      plon_v=nc_v.variables["longitude"][:]
      plat_v=nc_v.variables["latitude"][:]
      depth_v=nc_v.variables["depth"][:]
      check_grids(plon,plon_v,plat,plat_v)
      check_depths(depth,depth_v)

      # Time processing for tidal elevations 
      time_v=nc_v.variables["time"][:]
      tunit = nc_v.variables["time"].units
      mydt_v = cf_time_to_datetime(time_v,tunit) 

      # restriction for now, u and v must have same time steps as h
      # TODO: Loosen restriction
      try :
         difftu=[abs(diff_in_seconds(elem[0]-elem[1])) for elem in zip(mydt_h,mydt_u)]
         difftv=[abs(diff_in_seconds(elem[0]-elem[1])) for elem in zip(mydt_h,mydt_v)]
      except:
         # Probably due to size mismatch, but could be more descriptive. 
         # TODO: Add more descriptive error message
         msg="Error when subtracting times from u/v from h. Check your data"
         logger.error(msg)
         raise ValueError,msg

      #print difftu
      #print difftv
      if any([ elem > 10. for elem in difftu]) or any([ elem > 10. for elem in difftv]):
         msg="Times in tidal u/v vs tidal h mismatch. Time series must be estimated at the same times"
         logger.error(msg)
         raise ValueError,msg


   # Create output dir.
   path0=os.path.join(".","archv_with_tide")
   if os.path.exists(path0) and os.path.isdir(path0) :
      pass
   else :
      os.mkdir(path0)

   # Open blkdat files. Get some properties
   bp=modeltools.hycom.BlkdatParser("blkdat.input")
   idm    = bp["idm"]
   jdm    = bp["jdm"]
   kdm    = bp["kdm"]
   thflag = bp["thflag"]
   thbase = bp["thbase"]
   kapref = bp["kapref"]
   iversn = bp["iversn"]
   iexpt  = bp["iexpt"]
   yrflag = bp["yrflag"]
   thref=1e-3
   if kapref == -1 : 
      kapnum = 2
      msg="Only kapref>=0 is implemented for now"
      logger.error(msg)
      raise ValueError,msg
   else :
      kapnum = 1 

   if kapnum > 1 :
      msg="Only kapnum=1 is implemented for now"
      logger.error(msg)
      raise ValueError,msg


   # hycom sigma and kappa, written in python. NB: sigver is not used here.
   # Modify to use other equations of state. For now we assume sigver is:
   #    1 (7-term eqs referenced to    0 bar)
   #    2 (7-term eqs referenced to 2000 bar)
   if thflag == 0 :
      sigver=1
   else :
      sigver=2
   sig  = modeltools.hycom.Sigma(thflag)
   if kapref > 0  : kappa = modeltools.hycom.Kappa(kapref,thflag*1000.0e4) # 



   # Now loop through tide_times
   for rec,tide_time in enumerate(mydt_h) :

      # Construct archive file name to create
      iy = tide_time.year
      id,ih,isec = modeltools.hycom.datetime_to_ordinal(tide_time,yrflag)
      archv_file_in_string = "archv.%04d_%03d_%02d"%(iy,id,ih)

      # Is there match for this file name in list of archive files?
      I=[elem for elem in archv_files if os.path.basename(elem)[:17] == archv_file_in_string ]
      state_from_archv=len(I)>0
      if state_from_archv : archv_file_in =I[0]

      # Output file name
      fnameout = os.path.join(path0,os.path.basename(archv_file_in_string))
      arcfile_out=abfile.ABFileArchv(fnameout,"w",
            iversn=iversn,
            yrflag=yrflag,
            iexpt=iexpt,mask=False,
            cline1="TIDAL data has been added")

      tide_h=numpy.copy(nc_h.variables["h"][rec,:,:])
      tide_h=numpy.where(tide_h==nc_h.variables["h"]._FillValue,0.,tide_h)
      #print tide_h.min(),tide_h.max()
      if include_uv :
         tide_u=numpy.copy(nc_u.variables["u"][rec,:,:])
         tide_v=numpy.copy(nc_v.variables["v"][rec,:,:])
         #print tide_u.min(),tide_u.max()
         #print tide_v.min(),tide_u.max()

         tide_u=numpy.where(tide_u==nc_u.variables["u"]._FillValue,0.,tide_u)
         tide_v=numpy.where(tide_v==nc_v.variables["v"]._FillValue,0.,tide_v)

         # Rotate vectors to align with grid
         tide_u= tide_u*numpy.cos(pang) + tide_v*numpy.sin(pang)
         tide_v=-tide_u*numpy.sin(pang) + tide_v*numpy.cos(pang) #tide_v=tide_u*numpy.cos(pang+.5*numpy.pi) + tide_v*numpy.sin(pang+.5*numpy.pi)

         # From P-point to u. 2nd dim in python = 1st dim in Fortran
         tide_u[:,1:] =.5*(tide_u[:,1:] + tide_u[:,0:-1])
         tide_u=numpy.where(iu,tide_u,0.)

         # From P-point to v. 1st dim in python = 2nd dim in Fortran
         tide_v[1:,:] =.5*(tide_v[1:,:] + tide_v[0:-1,:])
         tide_v=numpy.where(iv,tide_v,0.)



      if state_from_archv :

         logger.info("Adding tidal values to existing state:%s"%arcfile_out.basename)
         arcfile=abfile.ABFileArchv(archv_file_in,"r")
         if arcfile.idm <> plon.shape[1] or  arcfile.jdm <> plon.shape[0] :
            msg="Grid size mismatch between %s and %s "%(tide_file,archv_file_in)

         # Read all layers .. (TODO: If there are memory problems, read and estimate sequentially)
         temp    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
         saln    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
         th3d  =numpy.ma.zeros((kdm,jdm,idm))
         thstar=numpy.ma.zeros((kdm,jdm,idm))
         dp    =numpy.ma.zeros((jdm,idm))
         p     =numpy.ma.zeros((kdm+1,jdm,idm))
         logger.info("Reading layers to get thstar and p")
         for k in range(kdm) :
            logger.debug("Reading layer %d from %s"%(k,archv_file_in))
            temp  =arcfile.read_field("temp",k+1)
            saln  =arcfile.read_field("salin",k+1)
            #dp    [k  ,:,:]=arcfile.read_field("thknss",k+1)
            dp    [:,:]=arcfile.read_field("thknss",k+1)
            th3d  [k  ,:,:]=sig.sig(temp,saln) - thbase
            p     [k+1,:,:]= p[k,:,:] + dp[:,:]
            thstar[k  ,:,:]=numpy.ma.copy(th3d  [k  ,:,:])
            if kapref > 0 :
               thstar[k  ,:,:]=thstar  [k  ,:,:] + kappa.kappaf(
                     temp[:,:], saln[:,:], th3d[k,:,:]+thbase, p[k,:,:])
            elif kapref < 0 :
               msg="Only kapref>=0 is implemented for now"
               logger.error(msg)
               raise ValueError,msg


         # Read montg1 and srfhgt, and set new values
         # ... we have ...
         # montg1 = montgc + montgpb * pbavg
         # srfhgt = montg1 + thref*pbavg
         # ...
         montg1  = arcfile.read_field("montg1",thflag)
         srfhgt  = arcfile.read_field("srfhgt",0)

         # New surface height - 
         montg1pb=modeltools.hycom.montg1_pb(thstar,p)
         montg1  = montg1 + montg1pb * modeltools.hycom.onem * tide_h
         srfhgt  = montg1 + thref*tide_h*modeltools.hycom.onem

         # Barotrpic velocities 
         if include_uv :
            ubavg  = arcfile.read_field("u_btrop",0)
            vbavg  = arcfile.read_field("v_btrop",0)
            ubavg  = ubavg + tide_u
            vbavg  = vbavg + tide_v

         # Loop through original fields and write
         for key in sorted(arcfile.fields.keys()) :
            fieldname = arcfile.fields[key]["field"]
            time_step = arcfile.fields[key]["step"]
            model_day = arcfile.fields[key]["day"]
            k         = arcfile.fields[key]["k"]
            dens      = arcfile.fields[key]["dens"]
            fld       =arcfile.read_field(fieldname,k)

            if fieldname == "montg1" :
               logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,fnameout))
               arcfile_out.write_field(montg1,None,fieldname,time_step,model_day,sigver,thbase) 
            elif fieldname == "srfhgt" :
               logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,fnameout))
               arcfile_out.write_field(srfhgt,None,fieldname,time_step,model_day,sigver,thbase) 
            elif fieldname == "u_btrop" and include_uv :
               logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,fnameout))
               arcfile_out.write_field(ubavg,None,fieldname,time_step,model_day,sigver,thbase) 
            elif fieldname == "v_btrop" and include_uv :
               logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,fnameout))
               arcfile_out.write_field(vbavg,None,fieldname,time_step,model_day,sigver,thbase) 
            else :
               arcfile_out.write_field(fld   ,None,fieldname,time_step,model_day,k,dens) 
               #logger.info("Writing field %10s at level %3d to %s (copy from original)"%(fieldname,k,fnameout))

         arcfile.close()


      else : 
         logger.info("Crating archv file with tidal data   :%s"%arcfile_out.basename)

         montg1=numpy.zeros((jdm,idm,))
         srfhgt=tide_h*modeltools.hycom.onem*thref
         arcfile_out.write_field(montg1,None,"montg1",0,0.,sigver,thbase) 
         arcfile_out.write_field(srfhgt,None,"srfhgt",0,0.,0,0.0) 
         if include_uv :
            ubavg  = tide_u
            vbavg  = tide_v
            arcfile_out.write_field(ubavg ,None,"ubavg" ,0,0.,0,0.0) 
            arcfile_out.write_field(vbavg ,None,"vbavg" ,0,0.,0,0.0) 



      logger.info("Finished writing to %s"%fnameout)
      arcfile_out.close()

   logger.info("Files containing tidal data in directory %s"%path0)
   logger.warning("Sigver assumed to be those of 7 term eqs")
   logger.warning("    1 for sigma-0/thflag=0, 2 for sigma-2/thflag=2")





if __name__ == "__main__" :
   class PointParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values[0].split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       tmp1= getattr(args, self.dest)
       tmp1.append(tmp)
       setattr(args, self.dest, tmp1)

   parser = argparse.ArgumentParser(description="""This routine will add
   previously generated tidal data to a archv file. The resulting file can be
   used as input to a nested hycom simulation""")
   parser.add_argument('--include-uv'  , action="store_true", default=False,
         help="Also add tidal u and v components")
   parser.add_argument('tide_file', type=str)
   parser.add_argument('archv', type=str,nargs="*")
   args = parser.parse_args()
   
   main(args.tide_file,args.archv,include_uv=args.include_uv)
