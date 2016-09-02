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


def main(blkdatfile,saltfile,lon,lat):

   logger.info("Salinity file:%s"%saltfile)
   m=re.match("(.*)_s([0-9]{2})_([0-9a-z]{4}\.nc)",saltfile)
   if m :
      tempfile = m.group(1)+"_t"+m.group(2)+"_"+m.group(3)
   else :
      msg="Cant deduce temp file from salt file name"
      logger.error(msg)
      raise ValueError,msg
   logger.info("Temp     file:%s"%saltfile)

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

   #Closest point
   s_lon=s_ncid.variables["lon"][:]
   s_lat=s_ncid.variables["lat"][:]
   j=numpy.argmin(numpy.abs(numpy.mod(s_lon-lon+360.,360.)))
   i=numpy.argmin(numpy.abs(s_lat-lat))
   print lon,s_lon[j]
   print lat,s_lat[i]

   dprof = s_ncid["depth"][:]
   salprof = s_ncid["s_an"][0,:,i,j]
   temprof = t_ncid["t_an"][0,:,i,j]
   try :
      nz=salprof[~salprof.mask].size
   except :
      nz=salprof.size


   # Interface values dprofi
   dprofi=numpy.zeros(dprof.size+1)
   for k in range(dprof.size) :
      if k < dprof.size -1 :
         dprofi[k] = 0.5 * (dprof[k]+dprof[k+1])
      else :
         dprofi[k] = dprof[k]

   
   cols="b"
   colt="r"
   cold="m"


   # MAx depth in data
   try :
      maxd = numpy.max(dprof[~salprof.mask])
   except :
      maxd = numpy.max(dprof)
   logger.info("MAx depth in data is %d m"%maxd)


   
   # Plot vertical profile of t and s in point
   f,axes = plt.subplots(2,2,sharey=True)
   ax1=axes[0,0]
   ax2=axes[0,1]
   ax3=axes[1,0]
   ax4=axes[1,1]
   ax1.set_title("Profiles at lon=%6.2f, lat=%6.2f"%(lon,lat))
   ax1.plot(temprof,-dprof,lw=2,color=colt)
   ax1.set_ylabel("Temperature[C]",color=colt)
   ax1.grid(True)
   for t in ax1.get_xticklabels() :
      t.set_color(colt)
      t.set_size(8)
      t.set_rotation(-45)
   ax2.plot(salprof,-dprof,lw=2,color=cols)
   ax2.set_ylabel("Salinity[psu]",color=cols)
   ax2.grid(True)
   for t in ax2.get_xticklabels() :
      t.set_color(cols)
      t.set_size(8)
      t.set_rotation(-45)
   ax3.plot(sig0.SIG(temprof,salprof),-dprof,lw=2,color=cold)
   ax3.set_ylabel("Density[sigma-0]",color=cold)
   ax3.grid(True)
   for t in ax3.get_xticklabels() :
      t.set_color(cold)
      t.set_size(8)
      t.set_rotation(-45)
   ax4.plot(sig2.SIG(temprof,salprof),-dprof,lw=2,color=cold)
   ax4.set_ylabel("Density[sigma-2]",color=cold)
   ax4.grid(True)
   for t in ax4.get_xticklabels() :
      t.set_color(cold)
      t.set_size(8)
      t.set_rotation(-45)
   plt.gcf().subplots_adjust()
   plt.gcf().savefig("ts_%s.png"%(sinfo),dpi=180)
   


   # Open blkdat.input
   bp=modeltools.hycom.BlkdatParser(blkdatfile)
   dp0k = bp["dp0k"]
   ds0k = bp["ds0k"]
   sigma = bp["sigma"]
   kdm=bp["kdm"]
   nhybrd=bp["nhybrd"]
   nsigma=bp["nsigma"]
   thflag=bp["thflag"]
   mysig   = modeltools.hycom.Sigma(thflag)
   sigprof = mysig.SIG(temprof,salprof)
   print dp0k

   # Interface from dp
   intf0k=numpy.zeros((kdm+1))
   intf0s=numpy.zeros((kdm+1))
   intf  =numpy.zeros((kdm+1))
   for i in range(nsigma) :
      intf0s[i+1] = intf0s[i] + ds0k[i]
   for i in range(nhybrd) :
      intf0k[i+1] = intf0k[i] + dp0k[i]

   # Shallow  level
   if   intf0s[nsigma] > maxd :
      intf=numpy.copy(intf0s)
      logger.info("Shallow z level")
   else :
      # Sigma level (stretched w depth)
      if intf0k[nsigma] > maxd:
         logger.info("Sigma level")
         fac = maxd/intf0k[nsigma] 
         intf=numpy.copy(intf0k) * fac
         print intf
      # Deep z level
      else :
         logger.info("Deep z level")
         intf=numpy.copy(intf0k)
   intf=numpy.minimum(intf,maxd)
   print "intf0s:",intf0s
   print "intf0k:",intf0k
   print "intf:",intf
   dp0=intf[1:]-intf[:-1]
   intfmid=(intf[1:]+intf[:-1])*.5

   # Plot vertical profile of sig-2 vs depth amd index ?
   f,(ax1,ax2) = plt.subplots(1,2,sharey=False)
   ax1.plot(sigma,-intfmid,lw=2,color=colt)
   ax1.set_ylabel("Min interface depth",color=colt)
   ax1.set_xlabel("sigma-%d"%thflag,color=colt)
   ax2.plot(sigma,-numpy.arange(intfmid.size),lw=2,color=colt)
   ax2.set_ylabel("Layer index",color=colt)
   ax2.set_xlabel("sigma-%d"%thflag,color=colt)
   #ax1.set_ylabel("Temperature[C]",color=colt)
   #ax1.grid(True)
   #for t in ax1.get_xticklabels() :
   #   t.set_color(colt)
   #   t.set_size(8)
   #   t.set_rotation(-45)
   plt.gcf().savefig("blkdat_sigma.png",dpi=180)


   # We have interface values, now use designated layer sigma values and actual sigma
   # values to find an approximate vertical coordinate setup. 

   # Loop over output layers
   newintf=numpy.zeros(intf.shape)
   intsig=numpy.zeros(dp0.shape)
   for k in range(kdm) :

      # Target layer upper interface
      upint=newintf[k]

      # Mix water over integration range
      dpsum=0.
      sgsum=0.
      for k2 in range(nz) :

         # Range of this layer
         upint2=dprofi[k2]
         lwint2=dprofi[k2+1]

         # Part of this layer in target layer
         u = max(upint2,upint)
         l = lwint2
         dp=max(0.,l-u)
         
         if dpsum > 0. : 

            # Integrated value of sigma up to this point
            sg = sgsum/dpsum

            # target layer lighter than integrated value. 
            # No need to add layers as they will only become heavier.
            if sigma[k] <= sg :
               pass

            # target layer heavier than integrated value. Ok to add more layers
            elif sigma[k] > sg :
               
               # Sum of layers at this point
               dpsum2=dpsum+dp
               sgsum2=sgsum+sigprof[k2]*dp
               sg2 = sgsum2/dpsum2

               # new layer > sigma[k]. Add only a part of the new layer
               if sigma[k] < sg2 :
                  dpfrac = (dpsum*sg - dpsum*sigma[k]) / (sigma[k] - sigprof[k2])
                  dpsum=dpsum+dpfrac
                  sgsum=sgsum+sigprof[k2]*dpfrac
               # new layer still < sigma[k]. Add full layer
               else :
                  dpsum=dpsum+dp
                  sgsum=sgsum+sigprof[k2]*dp

         # dpsum = 0.
         else :
            sg = sigprof[k2]
            # target layer lighter than layer value. 
            # No need to add layers as they will only become heavier.
            if sigma[k] <= sg :
               pass
            # target layer heavier than layer value. Ok to add more layers to mix downward
            else :
               #print "hei",sg,sigma[k2]
               dpsum=dpsum+dp
               sgsum=sgsum+sigprof[k2]*dp
            #print sg,sigma[k],dpsum,sgsum



      # End k2 loop

      # Make sure dpsum adheres to minimum layer thickness
      #print "0:",k,dpsum,sgsum
      dpsum=max(dpsum,dp0[k])
      #print "1:",k,dpsum
      
      # Adjust layer interface
      newintf[k+1] = newintf[k] + dpsum

      # MAke sure lowest interface is above sea floor
      newintf[k+1] = min(newintf[k+1],maxd)
      
      # Effective layer thickness
      dpsum = newintf[k+1] - newintf[k]

      # Estimated sigma
      intsig[k]=sgsum/dpsum
   # End k loop

   # Make sure bottom layer reaches sea floor
   newintf[kdm] = max(newintf[kdm],maxd)



   
   # Plot vertical profile and layers
   f,ax = plt.subplots(1)
   ax.set_title("Sigma-%d profile and vertical layers at lon=%6.2f, lat=%6.2f"%(thflag,lon,lat))
   ax.plot(sigprof,-dprof,lw=2,color=cold,label="Sigma-%d"%thflag)
   #ax.semilogy(sigprof,dprof,lw=2,color=cold)
   xlim=ax.get_xlim()
   ylim=ax.get_ylim()
   for k in range(kdm+1) :


      if k<kdm :
         ax.plot([-100,100],-newintf[k]*numpy.ones((2)),lw=.1,color=".5")
         textx = numpy.sum(xlim)/2.
         texty = -0.5*(newintf[k] + newintf[k+1])
         ax.text(textx,texty,"layer %d: %4.2f"%(k+1,sigma[k]),verticalalignment="center",horizontalalignment="center",fontsize=6)

         # Check if layer is isopycnal
         isopyc = numpy.abs(intsig[k]-sigma[k])<1e-4

         if isopyc :
            ax.fill_between([-100,100],-newintf[k+1]*numpy.ones((2)),-newintf[k]*numpy.ones((2)),color="g",alpha=".5")

      else :
         ax.plot([-100,100],-newintf[k]*numpy.ones((2)),lw=2,color="k",label="Sea floor")


   # Some dummy plots for the legend
   ax.plot([-100,100],[1e5,1e5],lw=.1,color=".5",label="Layer interface")
   ax.fill_between([-100,100],[1e5,1e5],[2e5,2e5],color="g",alpha=".5",label="Isopycnal layer") #NB: Label doesnt work


   
   ax.set_xlim(xlim)
   ax.set_ylim(ylim)
   ax.set_ylabel("Density[sigma-%d]"%thflag,color=cold)
   for t in ax.get_xticklabels() :
      t.set_color(cold)
      t.set_size(8)
      t.set_rotation(-45)
   #ax.legend(bbox_to_anchor=(1.2, 0.95))
   ax.legend(loc="lower left")
   plt.gcf().savefig("vert_%s.png"%sinfo,dpi=180)
   ax.set_ylim(-750,0)
   plt.gcf().savefig("vert750_%s.png"%sinfo,dpi=180)









   raise NameError,"test"




if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('blkdatfile' , help="")
   parser.add_argument('saltfile' , help="")
   parser.add_argument('lon'      , type=float) 
   parser.add_argument('lat'      , type=float) 

   args = parser.parse_args()
   main(args.blkdatfile,args.saltfile,args.lon,args.lat)



