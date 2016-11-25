#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime
import netCDF4
import os.path
import logging
import abfile
import numpy 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import scipy.interpolate
import struct
import matplotlib.mlab
import modeltools.tools

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


def plot_test(fld,filename) :
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   #P=ax.pcolormesh(fld)
   #P=ax.pcolormesh(fld[2200:2800,3500:4500],cmap=cmap)
   P=ax.pcolormesh(fld,cmap="jet")
   ax.figure.colorbar(P)
   figure.canvas.print_figure(filename)




def main(path,sigma,months,resolution) :
   logger.info("Path=%s"%path)
   dump_netcdf=True
   if months :
      months = [elem-1 for elem in months]
   else :
      months = range(12) 

   if resolution == 0.25 :
      fnametemplate="woa13_decav_%s%02d_04v2.nc" # 0.25 degrees
   else :
      fnametemplate="woa13_decav_%s%02d_01v2.nc" # 1.00 degrees


   # Open seasonal files
   ncid={}
   ncid["month"]={"t" :[],"s" : []}
   ncid["season"]={"t" :[],"s" : []}
   for i in range(4) :
      fname=os.path.join(path,fnametemplate%("s",i+13))
      logger.debug("Opening %s"%fname)
      ncid["season"]["s"].append(netCDF4.Dataset(fname,"r"))

      fname=os.path.join(path,fnametemplate%("t",i+13))
      logger.debug("Opening %s"%fname)
      ncid["season"]["t"].append(netCDF4.Dataset(fname,"r"))

   # Open monthly files
   seasonal_ncid={"t" :[],"s" : []}
   for i in range(12) :
      fname=os.path.join(path,fnametemplate%("s",i+1))
      ncid["month"]["s"].append(netCDF4.Dataset(fname,"r"))
      logger.debug("Opening %s"%fname)

      fname=os.path.join(path,fnametemplate%("t",i+1))
      ncid["month"]["t"].append(netCDF4.Dataset(fname,"r"))
      logger.debug("Opening %s"%fname)


   # Get lon_bnds from one of the files
   lon_bnds = ncid["month"]["s"][0].variables["lon_bnds"][:]
   lat_bnds = ncid["month"]["s"][0].variables["lat_bnds"][:]
   lon=ncid["month"]["s"][0].variables["lon"][:]
   lat=ncid["month"]["s"][0].variables["lat"][:]
   depths=ncid["season"]["s"][0].variables["depth"][:]
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


   sig = modeltools.hycom.Sigma(sigma)

   kkseason = ncid["season"]["s"][0].variables["s_an"].shape[1]
   kkmonth  = ncid["month"]["s"][0].variables["s_an"].shape[1]
   logger.info("kkseason=%3d, kkmonth=%3d"%(kkseason,kkmonth))


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
      f_saln = open("s_m%02d.d"%(month+1),"w")
      f_temp = open("t_m%02d.d"%(month+1),"w")
      f_dens = open("r_m%02d.d"%(month+1),"w")

      if dump_netcdf :
         fname_out_nc = "extrapolated_WOA2013_m%02d.nc"%(month + 1)
         ds_out       = netCDF4.Dataset(fname_out_nc, "w", format="NETCDF4")
         ds_out.createDimension("depth", kkseason)
         ds_out.createDimension("lon", nlon)
         ds_out.createDimension("lat", nlat)
         ds_out.createVariable("depth", "f4",("depth",))
         ds_out.createVariable("lat", "f4",("lat",))
         ds_out.createVariable("lon","f4",("lon",))
         ds_out.createVariable("temperature","f4",("depth","lat","lon",))
         ds_out.createVariable("salinity"   ,"f4",("depth","lat","lon",))
         ds_out.createVariable("temperature_e","f4",("depth","lat","lon",))
         ds_out.createVariable("salinity_e"   ,"f4",("depth","lat","lon",))
         ds_out.createVariable("density_e"    ,"f4",("depth","lat","lon",))
         ds_out.variables["lat"][:]=lat
         ds_out.variables["lon"][:]=lon
         ds_out.variables["depth"][:]=ncid["season"]["s"][0].variables["depth"][:]

      # season index and weights. Hardcoded, but possible to estimate from clim_bnds - file 1 is Jan, Feb, March, File 2 is April, MAy, June, etc...
      i0,i1,w0,w1 = month_weights(month+1) 

      # Write headers to output file
      monthname= datetime.date(1900, month+1, 1).strftime('%B')
      stitle="WOA 2013 Stable Salinity, %s"%monthname
      ttitle="WOA 2013 Stable Pot. Temp, %s"%monthname
      rtitle="WOA 2013 Stable Density, %s"%monthname
      stitle = stitle + " " * (40-len(stitle))
      ttitle = ttitle + " " * (40-len(ttitle))
      rtitle = rtitle + " " * (40-len(rtitle))

      fortran_write_title(f_saln,stitle)
      fortran_write_title(f_temp,ttitle)
      fortran_write_title(f_dens,rtitle)

      fortran_write_header(f_saln,nlon,nlat,lon[0],lat[0],dlon,dlat,depths)
      fortran_write_header(f_temp,nlon,nlat,lon[0],lat[0],dlon,dlat,depths)
      fortran_write_header(f_dens,nlon,nlat,lon[0],lat[0],dlon,dlat,depths)

      #raise NameError,"tes"


      logger.info("month = %2d/%2d, i0=%d, i1=%d, w0=%6.2f, w1=%6.2f"%(month+1,12,i0,i1,w0,w1))
      first=True
      d_over = numpy.zeros((nlat,nlon,))
      t_over = numpy.zeros((nlat,nlon,))
      s_over = numpy.zeros((nlat,nlon,))
      for k in range(kkseason) :
      #for k in range(1,kkseason,10) :
      #for k in range(kkseason-10,kkseason) :

         # Read 10 levels at a time

         logger.debug("Reading netcdf data")
         if k < kkmonth : 
            myfile="month"
            depth= ncid["month"]["s"][month].variables["depth"][k]

            s_out = ncid["month"]["s"][month].variables["s_an"][0,k,:]
            t_out = ncid["month"]["t"][month].variables["t_an"][0,k,:]
         else  :
            myfile="season"
            depth= ncid["season"]["s"][i0].variables["depth"][k]

            s_out_0 = ncid["season"]["s"][i0].variables["s_an"][0,k,:]
            s_out_1 = ncid["season"]["s"][i1].variables["s_an"][0,k,:]

            t_out_0 = ncid["season"]["t"][i0].variables["t_an"][0,k,:]
            t_out_1 = ncid["season"]["t"][i1].variables["t_an"][0,k,:]
            s_out=s_out_0*w0 + s_out_1*w1
            t_out=t_out_0*w0 + t_out_1*w1

         # write nc file
         if dump_netcdf :
            logger.info("Writing to netcdf file %s" % fname_out_nc)
            ds_out.variables["salinity"]   [k,:,:] = s_out
            ds_out.variables["temperature"][k,:,:] = t_out


         logger.info("%s file, level %3d/%3d, depth=%5.0f"%(myfile, k+1,kkseason,depth))
         # unmask using nearest neighbour for surface, linear deeper (faster). Since this is a one time job we do the slow stuff
         if first :
            method="nearest"
         else :
            #method="linear"
            method="nearest"
         logger.debug("Unmasking salinity using %s"%method)
         s_out = modeltools.tools.extrapolate_data(s_out,method)
         logger.debug("Unmasking temperature using %s"%method)
         t_out = modeltools.tools.extrapolate_data(t_out,method)

         # This will fill the remainder with values from above. This will happen if 
         # we move to climatology layers deeper than we have in the region
         if numpy.count_nonzero(s_out.mask)>0 or  numpy.count_nonzero(t_out.mask)>0 :
            logger.debug("Unmasking salinity using layer above")
            s_out[s_out.mask] = s_over[s_out.mask]
            logger.debug("Unmasking temperature using layer above")
            t_out[t_out.mask] = t_over[t_out.mask]

 
         # Density
         d_out = sig.sig(t_out,s_out)


         I=numpy.where(d_out<d_over)
         if I and len(I[0]) <> 0:
            d_out[I] = d_over[I]
            s_out[I] = sig.SOFSIG( d_out[I], t_out[I])
            logger.info("Corrected %d points for density"%len(I[0]))
         if first :
            plot_test(s_out,"s_out.png")
            plot_test(t_out,"t_out.png")
            plot_test(d_out,"d_out.png")

         d_out = numpy.ma.filled(d_out)
         t_out = numpy.ma.filled(t_out)
         s_out = numpy.ma.filled(s_out)

         fortran_write_data(f_saln,s_out)
         fortran_write_data(f_temp,t_out)
         fortran_write_data(f_dens,d_out)


         # write nc file
         if dump_netcdf :
            logger.info("Writing to netcdf file %s" % fname_out_nc)
            ds_out.variables["salinity_e"]   [k,:,:] = s_out
            ds_out.variables["temperature_e"][k,:,:] = t_out
            ds_out.variables["density_e"]    [k,:,:] = d_out
            ds_out.sync()

         # Some stats
         logger.info("**Salinity    min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (s_out.min(),s_out.max(),s_out.mean(),numpy.sqrt(numpy.mean(s_out**2))))
         logger.info("**Temperature min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (t_out.min(),t_out.max(),t_out.mean(),numpy.sqrt(numpy.mean(t_out**2))))
         logger.info("**Density     min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (d_out.min(),d_out.max(),d_out.mean(),numpy.sqrt(numpy.mean(d_out**2))))

         d_over=d_out
         s_over=s_out
         t_over=t_out
         first=False

      f_saln.close()
      f_temp.close()
      f_dens.close()
      if dump_netcdf :ds_out.close()












if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Routine prepares .d files, suitable to use as input to hycom relaxation setup')
   parser.add_argument('--resolution',type=float,help="Resolution of netcdf files",default=0.25)
   parser.add_argument('path',help="Directory where WOA13 netcdf files are present")
   parser.add_argument('sigma',type=int, help="Value of sigma to use (hycom thflag)")
   parser.add_argument('months',type=int,nargs="*",help="months to process for (1 to 12), all months processed if not present")
   args = parser.parse_args()

   # Set up AtmosphericForcing object, which keeps track of data to be read
   main(args.path,args.sigma,args.months,args.resolution)
