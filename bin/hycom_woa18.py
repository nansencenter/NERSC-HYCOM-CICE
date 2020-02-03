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
import modeltools.hycom
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





def main(path,resolution=0.25,write_netcdf=False) :

   if resolution == 0.25 :
      fnametemplate="woa13_decav_%s%02d_04v2.nc" # 0.25 degrees
   elif resolution == 1.00 :
      fnametemplate="woa13_decav_%s%02d_01v2.nc" # 1.00 degrees
   else :
      msg="Resolution must be 0.25 or 1.00"
      logger.error(msg)
      raise ValueError(msg)
 

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


   # Open regional grid file
   gridfile = abfile.ABFileGrid("regional.grid","r")
   plat=gridfile.read_field("plat")
   plon=gridfile.read_field("plon")
   gridfile.close()


   # Get lon_bnds from one of the files
   lon_bnds = ncid["month"]["s"][0].variables["lon_bnds"][:]
   lat_bnds = ncid["month"]["s"][0].variables["lat_bnds"][:]
   #print lon_bnds
   #print lat_bnds

   # Get grid spacing from bnds
   #print lon_bnds.shape
   dlon=lon_bnds[:,1]-lon_bnds[:,0]
   dlat=lat_bnds[:,1]-lat_bnds[:,0]
   if numpy.any(numpy.abs(dlon - dlon[0] )> numpy.abs(dlon[0] *1e-7)) :
      logger.error("longitude spacing not uniform")
   if numpy.any(numpy.abs(dlat - dlat[0] )> numpy.abs(dlat[0] *1e-7)) :
      logger.error("latgitude spacing not uniform")

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


   # Bounds ok, Find pivot points for hycom on woa grid
   lon0=lon_bnds[0,0]
   lat0=lat_bnds[0,0]
   dlon=dlon[0]
   dlat=dlat[0]
   tmp = numpy.fmod(720.+plon-lon0,360.)
   ipiv = tmp / dlon
   jpiv = (plat-lat0)/dlat
   plot_test(ipiv,"ipiv.png")
   plot_test(jpiv,"jpiv.png")

   # Find bilinear interpolation weights and cornerpoints
   s=ipiv-numpy.floor(ipiv)
   t=jpiv-numpy.floor(jpiv)
   wll=(1.-s)*(1.-t)
   wlr=    s *(1.-t)
   wur=    s *    t
   wul=(1.-s)*    t
   tmpsum= wll+wlr+wur+wul
   #print tmpsum.min(),tmpsum.max()
   plot_test(wll,"wll.png")
   ipiv=numpy.floor(ipiv).astype(int)
   jpiv=numpy.floor(jpiv).astype(int)
   print lat_bnds[-1,:]
   ipib=numpy.mod(ipiv+1,lon_bnds.shape[0])
   jpib=numpy.minimum(jpiv+1,lat_bnds.shape[0]-1)

   ## Test on actual field
   #fld = ncid["month"]["s"][0].variables["s_an"][0,0,:]
   #print fld.shape
   #newfld = fld[jpiv,ipiv]*wll + fld[jpiv,ipib]*wlr + fld[jpib,ipib]*wur + fld[jpib,ipiv]*wul
   #plot_test(newfld,"newfld.png")

   sig = modeltools.hycom.Sigma(0)

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




   # Loop over months
   for month in range(12) :
      # season index and weights. Hardcoded, but possible to estimate from clim_bnds - file 1 is Jan, Feb, March, File 2 is April, MAy, June, etc...
      i0,i1,w0,w1 = month_weights(month+1) 

      if  write_netcdf :
         fname_out_nc = "extrapolated_WOA2018_modelgrid_m%02d.nc"%(month+1)
         ds_out       = netCDF4.Dataset(fname_out_nc, "w", format="NETCDF4")
         ds_out.createDimension("depth", kkseason)
         ds_out.createDimension("idm", ipiv.shape[1])
         ds_out.createDimension("jdm", ipiv.shape[0])
         ds_out.createVariable("depth", "f4",("depth",))
         ds_out.createVariable("latitude", "f4",("jdm","idm",))
         ds_out.createVariable("longitude","f4",("jdm","idm",))
         ds_out.createVariable("temperature","f4",("depth","jdm","idm",))
         ds_out.createVariable("salinity"   ,"f4",("depth","jdm","idm",))
         ds_out.createVariable("temperature_e","f4",("depth","jdm","idm",))
         ds_out.createVariable("salinity_e"   ,"f4",("depth","jdm","idm",))
         ds_out.createVariable("density_e"    ,"f4",("depth","jdm","idm",))
         ds_out.variables["latitude"]=plat
         ds_out.variables["longitude"]=plon
         ds_out.variables["depth"][:]=ncid["season"]["s"][0].variables["depth"][:]


      # Open HYCOM .a files
      t_abfile= abfile.ABFileRelaxZ("temp_sig%d_m%02d"%(sig.sigma,month+1),"w",cline1="WOA2018 monthly",cline2="Potential temperature")
      s_abfile= abfile.ABFileRelaxZ("saln_sig%d_m%02d"%(sig.sigma,month+1),"w",cline1="WOA2018 monthly",cline2="Salinity")
      d_abfile= abfile.ABFileRelaxZ("dens_sig%d_m%02d"%(sig.sigma,month+1),"w",cline1="WOA2018 monthly",cline2="Potential density(Sigma-%d)"%sig.sigma)


      logger.info("month = %2d/%2d, i0=%d, i1=%d, w0=%6.2f, w1=%6.2f"%(month,12,i0,i1,w0,w1))
      d_over = numpy.zeros(ipiv.shape)
      t_over = numpy.zeros(ipiv.shape)
      s_over = numpy.zeros(ipiv.shape)
      first=True
      for k in range(kkseason) :
      #for k in range(1,kkseason,10) :
      #for k in range(kkseason-10,kkseason) :

         # Read 10 levels at a time

         logger.debug("Reading netcdf data")
         if k < kkmonth : 
            myfile="month"
            depth= ncid["month"]["s"][month].variables["depth"][k]

            s_in = ncid["month"]["s"][month].variables["s_an"][0,k,:]
            t_in = ncid["month"]["t"][month].variables["t_an"][0,k,:]
         else  :
            myfile="season"
            depth= ncid["season"]["s"][i0].variables["depth"][k]

            s_in_0 = ncid["season"]["s"][i0].variables["s_an"][0,k,:]
            s_in_1 = ncid["season"]["s"][i1].variables["s_an"][0,k,:]

            t_in_0 = ncid["season"]["t"][i0].variables["t_an"][0,k,:]
            t_in_1 = ncid["season"]["t"][i1].variables["t_an"][0,k,:]
            s_in=s_in_0*w0 + s_in_1*w1
            t_in=t_in_0*w0 + t_in_1*w1



         logger.info("%s file, level %3d/%3d, depth=%5.0f"%(myfile, k,kkseason,depth))
         s_out = s_in[jpiv,ipiv]*wll + s_in[jpiv,ipib]*wlr + s_in[jpib,ipib]*wur + s_in[jpib,ipiv]*wul
         t_out = t_in[jpiv,ipiv]*wll + t_in[jpiv,ipib]*wlr + t_in[jpib,ipib]*wur + t_in[jpib,ipiv]*wul
         t_out = numpy.maximum(-0.055*s_out,t_out)

         # write nc file
         if  write_netcdf :
            logger.info("Writing to netcdf file %s"%fname_out_nc)
            ds_out.variables["salinity"]   [k,:,:] = s_out
            ds_out.variables["temperature"][k,:,:] = t_out
            ds_out.sync()

         # unmask using nearest neighbour
         if first : 
            method="nearest"
         else :
            method="linear"

         # KAL - for low number of points this crashes with qhull error- figure out why
         # For now this should be relatively safe ... May catch exception in stead or may use other hull method...
         if numpy.count_nonzero(~s_out.mask) > 20 :
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
            first=False

         d_over = d_out
         t_over = t_out
         s_over = s_out


         # Write hycom abfile
         logger.info("Writing to hycom files")
         s_abfile.write_field(s_out,s_out,"salinity             ",depth)
         t_abfile.write_field(t_out,t_out,"potential temperature",depth)
         d_abfile.write_field(d_out,d_out,"sigma-%d"%sig.sigma   ,depth)

         # write nc file
         if  write_netcdf :
            logger.info("Writing to netcdf file %s"%fname_out_nc)
            ds_out.variables["salinity_e"]   [k,:,:] = s_out
            ds_out.variables["temperature_e"][k,:,:] = t_out
            ds_out.variables["density_e"]    [k,:,:] = d_out
            ds_out.sync()

         # Some stats
         logger.info("**Salinity    min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (s_out.min(),s_out.max(),s_out.mean(),numpy.sqrt(numpy.mean(s_out**2))))
         logger.info("**Temperature min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (t_out.min(),t_out.max(),t_out.mean(),numpy.sqrt(numpy.mean(t_out**2))))
         logger.info("**Density     min=%10.3f, max=%10.3f, ave=%14.3f, rms=%10.3f" % (d_out.min(),d_out.max(),d_out.mean(),numpy.sqrt(numpy.mean(d_out**2))))

         first = False

      s_abfile.close()
      d_abfile.close()
      t_abfile.close()


      if  write_netcdf : ds_out.close()











if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Create a z-level climatology interpolated to a model region from WOA2018')
   parser.add_argument('path')
   parser.add_argument('sigma',type=int)
   parser.add_argument('--resolution',type=float,default=0.25,help = "Resolution of data to use (0.25 or 1.00). Default is 0.25")
   parser.add_argument('--write-netcdf',action="store_true",default=False,help = "write to netcdf files as well as hycom .ab files")
   args = parser.parse_args()

   # Set up AtmosphericForcing object, which keeps track of data to be read
   main(args.path,
         resolution=args.resolution,
         write_netcdf=args.write_netcdf)
