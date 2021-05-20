#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import netCDF4
import os.path
import logging
import abfile
import numpy 
import scipy.interpolate
import struct
import modeltools.tools

'''
CAGLAR - Sep11,2020
Reads in annually averaged TCO2 and TAlk, extrapolates to fill the land points
and saves as monthly .d files.
USAGE: python hycom_glodap_zfiles.py path
e.g. : python hycom_glodap_zfiles.py /cluster/projects/nn2993k/ModelInput/GLODAPV2/ 
'''
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


def main(path,months) :
   logger.info("Path=%s"%path)
   dump_netcdf=True
   if months :
      months = [elem-1 for elem in months]
   else :
      months = range(12) 

   # Open annual files
   ncid={}
   ncid["annual"]={"TAlk" :[],"TCO2" : []}
   
   
   fname = path + "GLODAPv2.2016b.TAlk.nc"
   logger.debug("Opening %s"%fname)
   ncid["annual"]["TAlk"].append(netCDF4.Dataset(fname,"r"))

   fname = path + "GLODAPv2.2016b.TCO2.nc"
   logger.debug("Opening %s"%fname)
   ncid["annual"]["TCO2"].append(netCDF4.Dataset(fname,"r"))

   # Get lon_bnds from one of the files
   lon=ncid["annual"]["TAlk"][0].variables["lon"][:]
   lat=ncid["annual"]["TAlk"][0].variables["lat"][:]
   depths=ncid["annual"]["TAlk"][0].variables["Depth"][:]
   nlon=lon.size
   nlat=lat.size

   dlon=1.
   dlat=1.

   firstmonth = True

   kkannual = ncid["annual"]["TAlk"][0].variables["Depth"].shape[0]
   logger.info("kkannual=%3d"%(kkannual))


   def fortran_write_title(fid,title,endianness=">") :
      if len(title) <> 40 : 
         raise NameError,"title must be 40 chars"
      nbytes=len(title)
      fid.write(struct.pack("%si"%endianness,nbytes))
      fid.write(struct.pack("%s%dc"%(endianness,nbytes),*title))
      fid.write(struct.pack("%si"%endianness,nbytes))
      #print nbytes

   def fortran_write_header(fid,nlon,nlat,flon,flat,dlon,dlat,levs,endianness=">") :
      nbytes=(7+levs.size)*4
      fid.write(struct.pack("%si"%endianness,nbytes))
      fid.write(struct.pack("%s3i4f"%endianness,nlon,nlat,levs.size,flon,flat,dlon,dlat))
      fid.write(struct.pack("%s%df"%(endianness,levs.size),*levs))
      fid.write(struct.pack("%si"%endianness,nbytes))
      #print nbytes

   def fortran_write_data(fid,data,endianness=">"):
      nbytes=data.size*4
      fid.write(struct.pack("%si"%endianness,nbytes))
      data.astype("%sf"%endianness).tofile(fid)
      fid.write(struct.pack("%si"%endianness,nbytes))


   # Loop over months
   for month in months :
   
      # Open output files
      f_talk   = open("TAlk%02d.d"%(month+1),"w")
      f_co2 = open("TCO2%02d.d"%(month+1),"w")

      if firstmonth :
       if dump_netcdf :
         fname_out_nc = "extrapolated_GlODAPV2_annual.nc"
         ds_out       = netCDF4.Dataset(fname_out_nc, "w", format="NETCDF4")
         ds_out.createDimension("depth", kkannual)
         ds_out.createDimension("lon", nlon)
         ds_out.createDimension("lat", nlat)
         ds_out.createVariable("depth", "f4",("depth",))
         ds_out.createVariable("lat", "f4",("lat",))
         ds_out.createVariable("lon","f4",("lon",))
         ds_out.createVariable("TAlk","f4",("depth","lat","lon",))
         ds_out.createVariable("TCO2"   ,"f4",("depth","lat","lon",))
         ds_out.createVariable("TAlk_e","f4",("depth","lat","lon",))
         ds_out.createVariable("TCO2_e"   ,"f4",("depth","lat","lon",))
         ds_out.variables["lat"][:]=lat
         ds_out.variables["lon"][:]=lon
         ds_out.variables["depth"][:]=ncid["annual"]["TAlk"][0].variables["Depth"][:]

      # season index and weights. Hardcoded, but possible to estimate from clim_bnds - file 1 is Jan, Feb, March, File 2 is April, MAy, June, etc...
#      i0,i1,w0,w1 = month_weights(month+1) 
#      this is commented out since deeper levels for nutrients only exist in annual files, thus no month weighing required

      # Write headers to output file
      monthname= datetime.date(1900, month+1, 1).strftime('%B')
      ttitle="GloadV2.2 gridded annual TAlk"
      ctitle="GloadV2.2 gridded annual TCO2"
      ttitle = ttitle + " " * (40-len(ttitle))
      ctitle = ctitle + " " * (40-len(ctitle))

      fortran_write_title(f_talk,ttitle)
      fortran_write_title(f_co2,ctitle)

      # note that lon is written to header as lon - 20. so the longitudes are compatible with other variables
      # tracer data is rearranged below
      fortran_write_header(f_talk,nlon,nlat,lon[0]-20.,lat[0],dlon,dlat,depths)
      fortran_write_header(f_co2,nlon,nlat,lon[0]-20.,lat[0],dlon,dlat,depths)

      #raise NameError,"tes"


      first=True
      t_over = numpy.zeros((nlat,nlon,))
      c_over = numpy.zeros((nlat,nlon,))
      for k in range(kkannual) :

         logger.debug("Reading netcdf data")
         myfile="annual"
         depth= ncid["annual"]["TAlk"][0].variables["Depth"][k]

         t_out = ncid["annual"]["TAlk"][0].variables["TAlk"][k,:,:]
         c_out = ncid["annual"]["TCO2"][0].variables["TCO2"][k,:,:]
         # rearrange longitudes, so that these fit other hycom variable longitudes
         t_copy = numpy.zeros((nlat,nlon,))
         t_copy[:,0:20] = t_out[:,-20:]
         t_copy[:,-340:] = t_out[:,0:340]
         t_copy = numpy.ma.masked_where(t_copy<-990.,t_copy)
         t_out = t_copy

         c_copy = numpy.zeros((nlat,nlon,))
         c_copy[:,0:20] = c_out[:,-20:]
         c_copy[:,-340:] = c_out[:,0:340]
         c_copy = numpy.ma.masked_where(c_copy<-990.,c_copy)
         c_out = c_copy


         # write nc file
         if firstmonth:
          if dump_netcdf :
            logger.info("Writing to netcdf file %s" % fname_out_nc)
            ds_out.variables["TAlk"][k,:,:] = t_out
            ds_out.variables["TCO2"][k,:,:] = c_out


         logger.info("%s file, level %3d/%3d, depth=%5.0f"%(myfile, k+1,kkannual,depth))
         # unmask using nearest neighbour for surface, linear deeper (faster). Since this is a one time job we do the slow stuff
         if first :
            method="nearest"
         else :
            #method="linear"
            method="nearest"
         logger.debug("Unmasking silicate using %s"%method)
         t_out = modeltools.tools.extrapolate_data(t_out,method)
         logger.debug("Unmasking nitrate using %s"%method)
         c_out = modeltools.tools.extrapolate_data(c_out,method)


         # This will fill the remainder with values from above. This will happen if 
         # we move to climatology layers deeper than we have in the region
         if numpy.count_nonzero(t_out.mask)>0 : #or  numpy.count_nonzero(t_out.mask)>0 :
            logger.debug("Unmasking TAlk using layer above")
            t_out[t_out.mask] = t_over[t_out.mask]
            logger.debug("Unmasking TCO2 using layer above")
            c_out[c_out.mask] = c_over[c_out.mask]

 
         t_out = numpy.ma.filled(t_out)
         c_out = numpy.ma.filled(c_out)

         fortran_write_data(f_talk,t_out)
         fortran_write_data(f_co2,c_out)

         # write nc file
         if firstmonth:
          if dump_netcdf :
            logger.info("Writing to netcdf file %s" % fname_out_nc)
            ds_out.variables["TAlk_e"][k,:,:] = t_out
            ds_out.variables["TCO2_e"][k,:,:] = c_out
            ds_out.sync()

         # Some stats
         logger.info("**TAlk    min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (t_out.min(),t_out.max(),t_out.mean(),numpy.sqrt(numpy.mean(t_out**2))))
         logger.info("**TCO2 min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (c_out.min(),c_out.max(),c_out.mean(),numpy.sqrt(numpy.mean(c_out**2))))

         t_over=t_out
         c_over=c_out
         first=False
         if firstmonth:
            if dump_netcdf :ds_out.close()
         firstmonth=False

      f_talk.close()
      f_co2.close()





if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Routine prepares .d files, suitable to use as input to hycom relaxation setup')
   parser.add_argument('path',help="Directory where GLODAPV2 netcdf files are present")
   parser.add_argument('months',type=int,nargs="*",help="months to process for (1 to 12), all months processed if not present")
   args = parser.parse_args()

   main(args.path,args.months)
