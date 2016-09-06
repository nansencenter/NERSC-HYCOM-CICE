#!/usr/bin/env python
##!/usr/bin/python -E
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import modeltools.hycom
import logging
import argparse
import datetime
import numpy
import os
import netCDF4
import re
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

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


def main(saltfile,lon1,lon2,lat1,lat2,layerskip=1,minlayer=0,maxlayer=None):

   #TODO: Allow mupltiple salinity files ?

   logger.info("Salinity file:%s"%saltfile)
   m=re.match("(.*)_s([0-9]{2})_([0-9a-z]{4}\.nc)",saltfile)
   if m :
      tempfile = m.group(1)+"_t"+m.group(2)+"_"+m.group(3)
   else :
      msg="Cant deduce temp file from salt file name"
      logger.error(msg)
      raise ValueError,msg
   logger.info("Temp     file:%s"%tempfile)

   sig0 = modeltools.hycom.Sigma(0)
   sig2 = modeltools.hycom.Sigma(2)

   # output filename info from season value
   if int(m.group(2)) > 0 and int(m.group(2)) < 13  :
      sinfo=m.group(2)
   elif int(m.group(2)) == 0 :
      sinfo="year"
   elif int(m.group(2)) == 13 :
      sinfo="jfm"
   elif int(m.group(2)) == 14 :
      sinfo="amj"
   elif int(m.group(2)) == 15 :
      sinfo="jas"
   elif int(m.group(2)) == 16 :
      sinfo="ond"


   s_ncid = netCDF4.Dataset(saltfile,"r")
   t_ncid = netCDF4.Dataset(tempfile,"r")

   #Selected area:
   s_lon=s_ncid.variables["lon"][:]
   s_lat=s_ncid.variables["lat"][:]
   dprof = s_ncid["depth"][:]
   nz=dprof.size

   I=numpy.where(numpy.logical_and(s_lon>lon1,s_lon<lon2))
   J=numpy.where(numpy.logical_and(s_lat>lat1,s_lat<lat2))
   #print J
   #print I
   I,J=numpy.meshgrid(I,J)
   #print J
   #print I
   la = s_lon[I]
   lo = s_lat[J]



   if maxlayer is None :
      maxlayer=nz 
   else :
      maxlayer=min(maxlayer,nz)

   # Loop over vertical layers.  Skip layers if requested.
   nzs = maxlayer/layerskip +1 
   s0_max = numpy.zeros((nzs))
   s0_plo = numpy.zeros((nzs))
   s0_p10 = numpy.zeros((nzs))
   s0_p50 = numpy.zeros((nzs))
   s0_p90 = numpy.zeros((nzs))
   s0_phi = numpy.zeros((nzs))
   s0_min = numpy.zeros((nzs))

   s2_min = numpy.zeros((nzs))
   s2_plo = numpy.zeros((nzs))
   s2_p10 = numpy.zeros((nzs))
   s2_p50 = numpy.zeros((nzs))
   s2_p90 = numpy.zeros((nzs))
   s2_phi = numpy.zeros((nzs))
   s2_max = numpy.zeros((nzs))
   dprofs= numpy.zeros(nzs)

   plo=5
   phi=95


   t_bins=numpy.linspace(-2,30,50)
   s_bins=numpy.linspace(0,39,50)
   T_BINS,S_BINS = numpy.meshgrid(t_bins,s_bins)
   tsbin = numpy.zeros((t_bins.size,s_bins.size))
   for ks,k in enumerate(range(minlayer,maxlayer,layerskip)) :
      logger.info("Reading layer %d"%k)
      s = s_ncid["s_an"][0,k,:,:]
      s = s[J,I]
      t = t_ncid["t_an"][0,k,:,:]
      t = t[J,I]

      # For geo stats
      mylo=lo[~s.mask]
      myla=la[~s.mask]

      # Valid points
      s = s[~s.mask]
      t = t[~t.mask]
      #print s.shape

      # Densities:
      s0 =sig0.SIG(t,s)
      s2 =sig2.SIG(t,s)

      # Percentiles
      s0_plo[ks],s0_p10[ks],s0_p50[ks],s0_p90[ks],s0_phi[ks] = numpy.percentile(s0,[plo,10,50,90,phi])
      s0_max[ks] = numpy.max(s0)
      s0_min[ks] = numpy.min(s0)
      s2_plo[ks],s2_p10[ks],s2_p50[ks],s2_p90[ks],s2_phi[ks] = numpy.percentile(s2,[plo,10,50,90,phi])
      s2_max[ks] = numpy.max(s2)
      s2_min[ks] = numpy.min(s2)
      dprofs[ks] = dprof[k]

      # Bin salt and t. Find ways to speed up ....
      ind_s_bin=numpy.digitize(s,s_bins)
      ind_t_bin=numpy.digitize(t,t_bins)
      for i,j in zip(ind_s_bin,ind_t_bin) :
         tsbin[j,i] = tsbin[j,i] +1
   tsbin = numpy.ma.masked_where(tsbin<1,tsbin)

   # Scatter plot
   f1,ax1 = plt.subplots(1)
   #P=ax1.pcolormesh(s_bins,t_bins,tsbin,cmap=cmap)#,vmin=1,vmax=1000)
   P=ax1.pcolormesh(s_bins,t_bins,tsbin,
         norm=matplotlib.colors.LogNorm(vmin=tsbin.min(), vmax=tsbin.max()),
         cmap="summer")
   plt.colorbar(P,ax=ax1)
   cs0 = ax1.contour(S_BINS,T_BINS,sig0.SIG(T_BINS,S_BINS),20,colors="r",lw=.005,label="Sigma-0")
   plt.clabel(cs0, fontsize=6, inline=1,fmt="%1.1f")
   cs2=ax1.contour(S_BINS,T_BINS,sig2.SIG(T_BINS,S_BINS),20,colors="b",lw=.005,label="Sigma-2")
   plt.clabel(cs2, fontsize=6, inline=1,fmt="%1.1f")
   ax1.set_xlabel("Salinity")
   ax1.set_ylabel("Temperature")
   ax1.legend(loc="lower left")
   ax1.set_title("Densities from %d to %d m"%(dprof[minlayer],dprof[maxlayer-1]))
   f1.savefig("density_dist_bin.png",dpi=180)

   # Plot density distribution
   f1,ax1 = plt.subplots(1)
   ax1.fill_betweenx(-dprofs,s0_p10,s0_p90,-dprofs,color="g",alpha=".5")
   #ax1.plot(s0_min,-dprofs,linewidth=2,color="r")
   #ax1.plot(s0_max,-dprofs,linewidth=2,color="r")
   ax1.plot(s0_plo,-dprofs,linewidth=2,color="r",label="p%02d/p%02d"%(plo,phi))
   ax1.plot(s0_phi,-dprofs,linewidth=2,color="r")
   ax1.plot(s0_p10,-dprofs,linewidth=2,color="b",label="p10/p90")
   ax1.plot(s0_p50,-dprofs,linewidth=2,color="k",label="p50")
   ax1.plot(s0_p90,-dprofs,linewidth=2,color="b")
   ax1.grid(True)
   ax1.legend(loc="lower left")
   f1.savefig("vert_density_dist_sigma0.png",dpi=180)
   #
   f1,ax1 = plt.subplots(1)
   ax1.fill_betweenx(-dprofs,s2_p10,s2_p90,-dprofs,color="g",alpha=".5")
   #ax1.plot(s2_min,-dprofs,linewidth=2,color="r")
   #ax1.plot(s2_max,-dprofs,linewidth=2,color="r")
   ax1.plot(s2_plo,-dprofs,linewidth=2,color="r",label="p%02d/p%02d"%(plo,phi))
   ax1.plot(s2_phi,-dprofs,linewidth=2,color="r")
   ax1.plot(s2_p10,-dprofs,linewidth=2,color="b",label="p10/p90")
   ax1.plot(s2_p50,-dprofs,linewidth=2,color="k",label="p50")
   ax1.plot(s2_p90,-dprofs,linewidth=2,color="b")
   ax1.grid(True)
   ax1.legend(loc="lower left")
   f1.savefig("vert_density_dist_sigma2.png",dpi=180)





if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--layerskip'     , type=int,default=1) 
   parser.add_argument('--maxlayer'      , type=int,default=None) 
   parser.add_argument('--minlayer'      , type=int,default=0) 
   parser.add_argument('lon1'      , type=float) 
   parser.add_argument('lat1'      , type=float) 
   parser.add_argument('lon2'      , type=float) 
   parser.add_argument('lat2'      , type=float) 
   parser.add_argument('saltfile' , help="")

   args = parser.parse_args()
   main(args.saltfile,args.lon1,args.lon2,args.lat1,args.lat2,layerskip=args.layerskip,maxlayer=args.maxlayer,minlayer=args.minlayer)



