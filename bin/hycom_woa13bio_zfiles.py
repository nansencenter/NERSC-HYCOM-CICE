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
Reads in annually and monthly averaged WOA2013 nutrients, 
extrapolates to fill the land points and saves as monthly .d files.
Note that monthly averages are limited to mesopelagic, thus deep layers are 
taken from annual files.
USAGE: python hycom_woa13bio_zfiles.py path
e.g. : python hycom_woa13bio_zfiles.py /cluster/projects/nn2993k/ModelInput/WOA2013/Nutrients/ 
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

   fnametemplate="woa13_all_%s%02d_01.nc" # 1.00 degrees
   tracer = ['i','n','o','p']

   # Open annual files
   ncid={}
   ncid["month"]={"i" :[],"n" : [],"p" :[],"o" : []}
   ncid["annual"]={"i" :[],"n" : [],"p" :[],"o" : []}
   
   for item in tracer[:]:
       fname=os.path.join(path,fnametemplate%(item,0))
       logger.debug("Opening %s"%fname)
       ncid["annual"][item].append(netCDF4.Dataset(fname,"r"))

       # Open monthly files
       for i in range(12) :
           fname=os.path.join(path,fnametemplate%(item,i+1))
           ncid["month"][item].append(netCDF4.Dataset(fname,"r"))
           logger.debug("Opening %s"%fname)

   # Get lon_bnds from one of the files
   lon_bnds = ncid["month"]["i"][0].variables["lon_bnds"][:]
   lat_bnds = ncid["month"]["i"][0].variables["lat_bnds"][:]
   lon=ncid["month"]["i"][0].variables["lon"][:]
   lat=ncid["month"]["i"][0].variables["lat"][:]
   depths=ncid["annual"]["i"][0].variables["depth"][:]
   nlon=lon.size
   nlat=lat.size

   # Get grid spacing from bnds
   dlon=lon_bnds[:,1]-lon_bnds[:,0]
   dlat=lat_bnds[:,1]-lat_bnds[:,0]
   if numpy.any(numpy.abs(dlon - dlon[0] )> numpy.abs(dlon[0] *1e-7)) :
      logger.error("longitude spacing not uniform")
   if numpy.any(numpy.abs(dlat - dlat[0] )> numpy.abs(dlat[0] *1e-7)) :
      logger.error("latgitude spacing not uniform")
   dlon=dlon[0]
   dlat=dlat[0]

   # Consistency checks for bnds
   for mkey in ncid.keys() :
      for vkey in ncid[mkey].keys() :
         for (i,tmpnc) in enumerate(ncid[mkey][vkey]) :
            #print mkey,vkey,i,tmpnc.filepath()
            tmplo = tmpnc.variables["lon_bnds"]
            tmpla = tmpnc.variables["lat_bnds"]

            if any([elem[1]-elem[0] <> 0 for elem in zip(lon_bnds.shape,tmplo.shape)]) :
               loggear.error("longitude shapes differ between files")

            if any([elem[1]-elem[0] <> 0 for elem in zip(lat_bnds.shape,tmpla.shape)]) :
               logger.error("latitude shapes differ between files")

            if numpy.any(numpy.abs(lon_bnds - tmplo )> 1e-7) :
               logger.error("longitudes differ between files")

            if numpy.any(numpy.abs(lat_bnds - tmpla )> 1e-7) :
               logger.error("latitudes differ between files")



   kkannual = ncid["annual"]["i"][0].variables["i_an"].shape[1]
   kkmonth  = ncid["month"]["i"][0].variables["i_an"].shape[1]
   logger.info("kkannual=%3d, kkmonth=%3d"%(kkannual,kkmonth))


   #NB: Not general

   def month_weights(month) :
         i0=((month+1-2)/3)%4
         i1=(i0+1)%4
         month0=i0*3+2
         w1=((month+1-month0)%12)/3.
         w0=1.-w1
         return i0,i1,w0,w1



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
      f_nitrate   = open("n%02d.d"%(month+1),"w")
      f_phosphate = open("p%02d.d"%(month+1),"w")
      f_silicate  = open("i%02d.d"%(month+1),"w")
      f_oxygen    = open("o%02d.d"%(month+1),"w")

      if dump_netcdf :
         fname_out_nc = "extrapolated_WOA2013_m%02d.nc"%(month + 1)
         ds_out       = netCDF4.Dataset(fname_out_nc, "w", format="NETCDF4")
         ds_out.createDimension("depth", kkannual)
         ds_out.createDimension("lon", nlon)
         ds_out.createDimension("lat", nlat)
         ds_out.createVariable("depth", "f4",("depth",))
         ds_out.createVariable("lat", "f4",("lat",))
         ds_out.createVariable("lon","f4",("lon",))
         ds_out.createVariable("nitrate","f4",("depth","lat","lon",))
         ds_out.createVariable("silicate"   ,"f4",("depth","lat","lon",))
         ds_out.createVariable("phosphate","f4",("depth","lat","lon",))
         ds_out.createVariable("oxygen"   ,"f4",("depth","lat","lon",))
         ds_out.createVariable("nitrate_e","f4",("depth","lat","lon",))
         ds_out.createVariable("silicate_e"   ,"f4",("depth","lat","lon",))
         ds_out.createVariable("phosphate_e","f4",("depth","lat","lon",))
         ds_out.createVariable("oxygen_e"   ,"f4",("depth","lat","lon",))
         ds_out.variables["lat"][:]=lat
         ds_out.variables["lon"][:]=lon
         ds_out.variables["depth"][:]=ncid["annual"]["i"][0].variables["depth"][:]

      # season index and weights. Hardcoded, but possible to estimate from clim_bnds - file 1 is Jan, Feb, March, File 2 is April, MAy, June, etc...
#      i0,i1,w0,w1 = month_weights(month+1) 
#      this is commented out since deeper levels for nutrients only exist in annual files, thus no month weighing required

      # Write headers to output file
      monthname= datetime.date(1900, month+1, 1).strftime('%B')
      ititle="WOA 2013" + " %s "%(monthname) + "Silicate"
      ntitle="WOA 2013" + " %s "%(monthname) + "Nitrate"
      ptitle="WOA 2013" + " %s "%(monthname) + "Phosphate"
      otitle="WOA 2013" + " %s "%(monthname) + "Oxygen"
      ititle = ititle + " " * (40-len(ititle))
      ntitle = ntitle + " " * (40-len(ntitle))
      ptitle = ptitle + " " * (40-len(ptitle))
      otitle = otitle + " " * (40-len(otitle))

      fortran_write_title(f_silicate,ititle)
      fortran_write_title(f_nitrate,ntitle)
      fortran_write_title(f_phosphate,ptitle)
      fortran_write_title(f_oxygen,otitle)

      fortran_write_header(f_silicate,nlon,nlat,lon[0],lat[0],dlon,dlat,depths)
      fortran_write_header(f_nitrate,nlon,nlat,lon[0],lat[0],dlon,dlat,depths)
      fortran_write_header(f_phosphate,nlon,nlat,lon[0],lat[0],dlon,dlat,depths)
      fortran_write_header(f_oxygen,nlon,nlat,lon[0],lat[0],dlon,dlat,depths)

      #raise NameError,"tes"


      first=True
      i_over = numpy.zeros((nlat,nlon,))
      n_over = numpy.zeros((nlat,nlon,))
      p_over = numpy.zeros((nlat,nlon,))
      o_over = numpy.zeros((nlat,nlon,))
      for k in range(kkannual) :

         logger.debug("Reading netcdf data")
         if k < kkmonth : 
            myfile="month"
            depth= ncid["month"]["i"][month].variables["depth"][k]

            i_out = ncid["month"]["i"][month].variables["i_an"][0,k,:]
            n_out = ncid["month"]["n"][month].variables["n_an"][0,k,:]
            p_out = ncid["month"]["p"][month].variables["p_an"][0,k,:]
            o_out = ncid["month"]["o"][month].variables["o_an"][0,k,:]
         else  :
            myfile="annual"
            depth= ncid["annual"]["i"][0].variables["depth"][k]

            i_out = ncid["annual"]["i"][0].variables["i_an"][0,k,:]
            n_out = ncid["annual"]["n"][0].variables["n_an"][0,k,:]
            p_out = ncid["annual"]["p"][0].variables["p_an"][0,k,:]
            o_out = ncid["annual"]["o"][0].variables["o_an"][0,k,:]


         # write nc file
         if dump_netcdf :
            logger.info("Writing to netcdf file %s" % fname_out_nc)
            ds_out.variables["silicate"][k,:,:] = i_out
            ds_out.variables["nitrate"][k,:,:] = n_out
            ds_out.variables["phosphate"][k,:,:] = p_out
            ds_out.variables["oxygen"][k,:,:] = o_out


         logger.info("%s file, level %3d/%3d, depth=%5.0f"%(myfile, k+1,kkannual,depth))
         # unmask using nearest neighbour for surface, linear deeper (faster). Since this is a one time job we do the slow stuff
         if first :
            method="nearest"
         else :
            #method="linear"
            method="nearest"
         logger.debug("Unmasking silicate using %s"%method)
         i_out = modeltools.tools.extrapolate_data(i_out,method)
         logger.debug("Unmasking nitrate using %s"%method)
         n_out = modeltools.tools.extrapolate_data(n_out,method)
         logger.debug("Unmasking phosphate using %s"%method)
         p_out = modeltools.tools.extrapolate_data(p_out,method)
         logger.debug("Unmasking oxygen using %s"%method)
         o_out = modeltools.tools.extrapolate_data(o_out,method)


         # This will fill the remainder with values from above. This will happen if 
         # we move to climatology layers deeper than we have in the region
         if numpy.count_nonzero(i_out.mask)>0 : #or  numpy.count_nonzero(t_out.mask)>0 :
            logger.debug("Unmasking silicate using layer above")
            i_out[i_out.mask] = i_over[i_out.mask]
            logger.debug("Unmasking nitrate using layer above")
            n_out[n_out.mask] = n_over[n_out.mask]
            logger.debug("Unmasking phosphate using layer above")
            p_out[p_out.mask] = p_over[p_out.mask]
            logger.debug("Unmasking oxygen using layer above")
            o_out[o_out.mask] = o_over[o_out.mask]

 
         i_out = numpy.ma.filled(i_out)
         n_out = numpy.ma.filled(n_out)
         p_out = numpy.ma.filled(p_out)
         o_out = numpy.ma.filled(o_out)

         fortran_write_data(f_silicate,i_out)
         fortran_write_data(f_nitrate,n_out)
         fortran_write_data(f_phosphate,p_out)
         fortran_write_data(f_oxygen,o_out)

         # write nc file
         if dump_netcdf :
            logger.info("Writing to netcdf file %s" % fname_out_nc)
            ds_out.variables["silicate_e"]   [k,:,:] = i_out
            ds_out.variables["nitrate_e"][k,:,:] = n_out
            ds_out.variables["phosphate_e"]    [k,:,:] = p_out
            ds_out.variables["oxygen_e"]    [k,:,:] = o_out
            ds_out.sync()

         # Some stats
         logger.info("**Silicate    min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (i_out.min(),i_out.max(),i_out.mean(),numpy.sqrt(numpy.mean(i_out**2))))
         logger.info("**Nitrate min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (n_out.min(),n_out.max(),n_out.mean(),numpy.sqrt(numpy.mean(n_out**2))))
         logger.info("**Phosphate     min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (p_out.min(),p_out.max(),p_out.mean(),numpy.sqrt(numpy.mean(p_out**2))))
         logger.info("**Oxygen     min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (o_out.min(),o_out.max(),o_out.mean(),numpy.sqrt(numpy.mean(o_out**2))))

         i_over=i_out
         n_over=n_out
         p_over=p_out
         o_over=o_out
         first=False

      f_nitrate.close()
      f_silicate.close()
      f_phosphate.close()
      f_oxygen.close()
      if dump_netcdf :ds_out.close()


if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Routine prepares .d files, suitable to use as input to hycom relaxation setup')
   parser.add_argument('path',help="Directory where WOA13 netcdf files are present")
   parser.add_argument('months',type=int,nargs="*",help="months to process for (1 to 12), all months processed if not present")
   args = parser.parse_args()

   main(args.path,args.months)
