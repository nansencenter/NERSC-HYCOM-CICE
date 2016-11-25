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


def main(saltfile,lon1,lon2,lat1,lat2,
      layerskip=1,
      minlayer=0,
      maxlayer=None,
      blkdatfile="",
      dpi=180,
      showminmax=False,
      extreme_percentile=2):

   if blkdatfile :
      bp=modeltools.hycom.BlkdatParser(blkdatfile)

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
   lo = s_lon[I]
   la = s_lat[J]



   if maxlayer is None :
      maxlayer=nz 
   else :
      maxlayer=min(maxlayer,nz)


   # Loop over vertical layers.  Skip layers if requested.
   nzs = maxlayer/layerskip
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
   #print maxlayer,nzs,nz

   plo=extreme_percentile
   phi=100-plo


   s0_min_map=numpy.ones(I.shape)*1e8
   s2_min_map=numpy.ones(I.shape)*1e8
   s0_max_map=numpy.ones(I.shape)*-1e8
   s2_max_map=numpy.ones(I.shape)*-1e8

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
      mask0=~s.mask

      # For geo stats
      mylo=lo[~s.mask]
      myla=la[~s.mask]

      # Valid points
      s = s[~s.mask]
      t = t[~t.mask]

      if s.size == 0 : 
         break

      # Densities:
      s0 =sig0.sig(t,s)
      s2 =sig2.sig(t,s)


      s0_min_map[mask0]=numpy.minimum(s0_min_map[mask0],s0)
      s0_max_map[mask0]=numpy.maximum(s0_max_map[mask0],s0)
      s2_min_map[mask0]=numpy.minimum(s2_min_map[mask0],s2)
      s2_max_map[mask0]=numpy.maximum(s2_max_map[mask0],s2)


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

      maxk=k
   tsbin = numpy.ma.masked_where(tsbin<1,tsbin)

   s0_min_map=numpy.ma.masked_where(numpy.abs(s0_min_map)>1e6,s0_min_map)
   s0_max_map=numpy.ma.masked_where(numpy.abs(s0_max_map)>1e6,s0_max_map)
   s2_min_map=numpy.ma.masked_where(numpy.abs(s0_min_map)>1e6,s2_min_map)
   s2_max_map=numpy.ma.masked_where(numpy.abs(s2_max_map)>1e6,s2_max_map)


   dprofs=dprofs[0:maxk]
   s0_min=s0_min[0:maxk]
   s0_max=s0_max[0:maxk]
   s0_plo=s0_plo[0:maxk]
   s0_phi=s0_phi[0:maxk]
   s0_p50=s0_p50[0:maxk]
   s0_p10=s0_p10[0:maxk]
   s0_p90=s0_p90[0:maxk]
   #
   s2_min=s2_min[0:maxk]
   s2_max=s2_max[0:maxk]
   s2_plo=s2_plo[0:maxk]
   s2_phi=s2_phi[0:maxk]
   s2_p50=s2_p50[0:maxk]
   s2_p10=s2_p10[0:maxk]
   s2_p90=s2_p90[0:maxk]


   # Scatter plot
   f1,ax1 = plt.subplots(1)
   #P=ax1.pcolormesh(s_bins,t_bins,tsbin,cmap=cmap)#,vmin=1,vmax=1000)
   P=ax1.pcolormesh(s_bins,t_bins,tsbin,
         norm=matplotlib.colors.LogNorm(vmin=tsbin.min(), vmax=tsbin.max()),
         cmap="summer")
   plt.colorbar(P,ax=ax1)
   cs0 = ax1.contour(S_BINS,T_BINS,sig0.sig(T_BINS,S_BINS),20,colors="r",lw=.005,label="Sigma-0")
   plt.clabel(cs0, fontsize=6, inline=1,fmt="%1.1f")
   cs2=ax1.contour(S_BINS,T_BINS,sig2.sig(T_BINS,S_BINS),20,colors="b",lw=.005,label="Sigma-2")
   plt.clabel(cs2, fontsize=6, inline=1,fmt="%1.1f")
   #

   if blkdatfile and bp["thflag"]== 0 :
      cs3=ax1.contour(S_BINS,T_BINS,sig0.sig(T_BINS,S_BINS),bp["sigma"],colors="k",lw=.010,label="Sigma-2")
   else :
      cs3=ax1.contour(S_BINS,T_BINS,sig2.sig(T_BINS,S_BINS),bp["sigma"],colors="k",lw=.010,label="Sigma-2")
   #
   ax1.set_xlabel("Salinity")
   ax1.set_ylabel("Temperature")
   ax1.legend(loc="lower left")
   ax1.set_title("Binned T,S,Sigma-0, and Sigma-2 from %d to %d m [$kg/m^3$]"%(dprof[minlayer],dprof[maxlayer-1]))
   f1.canvas.print_figure("density_dist_bin.png",dpi=dpi)


   # Plot density distribution
   f1,ax1 = plt.subplots(1)
   if blkdatfile and bp["thflag"]== 0 :
      for k in range(bp["kdm"]) :
         if k==0 :
            label="blkdat"
         else :
            label=None
         ax1.plot(bp["sigma"][k]*numpy.ones(2),[dprofs[0],dprofs[-1]],color="m",label=label,lw=.5)
   ax1.fill_betweenx(dprofs,s0_p10,s0_p90,color="g",alpha=".5")
   if showminmax :ax1.plot(s0_min,dprofs,"k+",linewidth=2,label="min/max")
   if showminmax :ax1.plot(s0_max,dprofs,"k+",linewidth=2)
   ax1.plot(s0_plo,dprofs,linewidth=2,color="r",label="p%02d/p%02d"%(plo,phi))
   ax1.plot(s0_phi,dprofs,linewidth=2,color="r")
   ax1.plot(s0_p10,dprofs,linewidth=2,color="b",label="p10/p90")
   ax1.plot(s0_p50,dprofs,linewidth=2,color="k",label="p50")
   ax1.plot(s0_p90,dprofs,linewidth=2,color="b")
   ax1.text(0.05,0.95,"min: %5.2f"%(s0_min.min()),transform=ax1.transAxes)
   ax1.text(0.05,0.91,"max: %5.2f"%(s0_max.max()),transform=ax1.transAxes)
   ax1.grid(True)
   ax1.set_yscale("log")
   ax1.invert_yaxis()
   ax1.legend(loc="lower left")
   ax1.set_ylabel("Depth[m]")
   ax1.set_xlabel("Sigma-0")
   xlim=ax1.get_xlim()
   ax1.xaxis.set_ticks(numpy.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/30))
   for t in ax1.get_xticklabels() :
      t.set_size(8)
      t.set_rotation(-90)
   for line in ax1.get_xgridlines() +ax1.get_ygridlines() :
      line.set_linewidth(.1)
      line.set_linestyle("-")
      line.set_color(".5")
   f1.canvas.print_figure("sigma0_vdist.png",dpi=dpi)

   #
   f1,ax1 = plt.subplots(1)
   if blkdatfile and bp["thflag"]== 2 :
      for k in range(bp["kdm"]) :
         if k==0 :
            label="blkdat"
         else :
            label=None
         ax1.plot(bp["sigma"][k]*numpy.ones(2),[dprofs[0],dprofs[-1]],color="m",label=label,lw=.5)
   ax1.fill_betweenx(dprofs,s2_p10,s2_p90,color="g",alpha=".5")
   if showminmax :ax1.plot(s2_min,dprofs,"k+",linewidth=2,label="min/max")
   if showminmax :ax1.plot(s2_max,dprofs,"k+",linewidth=2)
   ax1.plot(s2_plo,dprofs,linewidth=2,color="r",label="p%02d/p%02d"%(plo,phi))
   ax1.plot(s2_phi,dprofs,linewidth=2,color="r")
   ax1.plot(s2_p10,dprofs,linewidth=2,color="b",label="p10/p90")
   ax1.plot(s2_p50,dprofs,linewidth=2,color="k",label="p50")
   ax1.plot(s2_p90,dprofs,linewidth=2,color="b")
   ax1.text(0.05,0.95,"min: %5.2f"%(s2_min.min()),transform=ax1.transAxes)
   ax1.text(0.05,0.91,"max: %5.2f"%(s2_max.max()),transform=ax1.transAxes)
   ax1.grid(True)
   ax1.set_yscale("log")
   ax1.invert_yaxis()
   ax1.legend(loc="lower left")
   ax1.set_ylabel("Depth[m]")
   ax1.set_xlabel("Sigma-2")
   ax1.legend(loc="lower left")
   #ax1.set_xlim(s2_plo.min(),s2_phi.max())
   xlim=ax1.get_xlim()
   ax1.xaxis.set_ticks(numpy.arange(xlim[0],xlim[1],(xlim[1]-xlim[0])/30))
   for t in ax1.get_xticklabels() :
      t.set_size(8)
      t.set_rotation(-90)
   for line in ax1.get_xgridlines() +ax1.get_ygridlines() :
      line.set_linewidth(.1)
      line.set_linestyle("-")
      line.set_color(".5")
   f1.canvas.print_figure("sigma2_vdist.png",dpi=dpi)


   #Plot of maps
   cmap="Paired"
   f1,ax1 = plt.subplots(1)
   vmin=s0_min_map.min()
   vmax=s0_min_map.max()
   vmax=vmin+(vmax-vmin)*.20
   P=ax1.pcolormesh(lo,la,s0_min_map,
         #norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax),
         vmin=vmin,vmax=vmax,
         cmap=cmap)
   f1.colorbar(P,ax=ax1,extend="both")
   ax1.set_xlabel("Longitude")
   ax1.set_ylabel("Latitude")
   ax1.set_title("Minimum sigma-0 density in water column")
   f1.canvas.print_figure("sigma0_min.png",dpi=dpi)
   #
   f1,ax1 = plt.subplots(1)
   vmin=s0_max_map.min()
   vmax=s0_max_map.max()
   vmin=vmax-(vmax-vmin)*.20
   P=ax1.pcolormesh(lo,la,s0_max_map,
         #norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax),
         vmin=vmin,vmax=vmax,
         cmap=cmap)
   f1.colorbar(P,ax=ax1,extend="both")
   ax1.set_xlabel("Longitude")
   ax1.set_ylabel("Latitude")
   ax1.set_title("Maximum sigma-0 density in water column")
   f1.canvas.print_figure("sigma0_max.png",dpi=dpi)



   #Plot of maps
   cmap="Paired"
   f1,ax1 = plt.subplots(1)
   vmin=s2_min_map.min()
   vmax=s2_min_map.max()
   vmax=vmin+(vmax-vmin)*.20
   P=ax1.pcolormesh(lo,la,s2_min_map,
         #norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax),
         vmin=vmin,vmax=vmax,
         cmap=cmap)
   f1.colorbar(P,ax=ax1,extend="both")
   ax1.set_xlabel("Longitude")
   ax1.set_ylabel("Latitude")
   ax1.set_title("Minimum sigma-2 density in water column")
   f1.canvas.print_figure("sigma2_min.png",dpi=dpi)
   #
   f1,ax1 = plt.subplots(1)
   vmin=s2_max_map.min()
   vmax=s2_max_map.max()
   vmin=vmax-(vmax-vmin)*.20
   P=ax1.pcolormesh(lo,la,s2_max_map,
         #norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax),
         vmin=vmin,vmax=vmax,
         cmap=cmap)
   f1.colorbar(P,ax=ax1,extend="both")
   ax1.set_xlabel("Longitude")
   ax1.set_ylabel("Latitude")
   ax1.set_title("Maximum sigma-2 density in water column")
   f1.canvas.print_figure("sigma2_max.png",dpi=dpi)


if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--dpi'     , type=int,default=180,help="dpi of output plots") 
   parser.add_argument('--layerskip'     , type=int,default=1,help="number of layers to jump when reading data to be analyzed") 
   parser.add_argument('--maxlayer'      , type=int,default=None,help="max layer to extract data for") 
   parser.add_argument('--minlayer'      , type=int,default=0,help="min layer to extract data for") 
   parser.add_argument('--blkdatfile'      , type=str,default="",help="If specified, some plots will be populated by data from blkdat.input") 
   parser.add_argument('lon1'      , type=float,help="Lower left longitude of region") 
   parser.add_argument('lat1'      , type=float,help="Lower left latitude  of region") 
   parser.add_argument('lon2'      , type=float,help="Upper left longitude of region") 
   parser.add_argument('lat2'      , type=float,help="Upper left latitude  of region") 
   parser.add_argument('saltfile' , help="WOA2013 salinity file")

   args = parser.parse_args()
   main(args.saltfile,args.lon1,args.lon2,args.lat1,args.lat2,
         layerskip=args.layerskip,
         maxlayer=args.maxlayer,
         minlayer=args.minlayer,
         blkdatfile=args.blkdatfile,
         dpi=args.dpi)



