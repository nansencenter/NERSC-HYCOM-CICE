#!/usr/bin/env python
##!/usr/bin/python -E
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import modeltools.hycom
import gridxsec
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


def main(blkdatfile,saltfile,lon,lat,lon2=None,lat2=None,sectionid=None,dpi=180):

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

   if sectionid : sinfo=sinfo+"_%s"%sectionid


   s_ncid = netCDF4.Dataset(saltfile,"r")
   t_ncid = netCDF4.Dataset(tempfile,"r")
   s_lon=s_ncid.variables["lon"][:]
   s_lat=s_ncid.variables["lat"][:]
   dprof = s_ncid["depth"][:]
   nz=dprof.size

   #Closest point
   if lon2 is None or lat2 is None :
      section=False
      i=numpy.argmin(numpy.abs(numpy.mod(s_lon-lon+360.,360.)))
      j=numpy.argmin(numpy.abs(s_lat-lat))

      salprof = s_ncid["s_an"][0,:,j,i]
      temprof = t_ncid["t_an"][0,:,j,i]
      salprof=numpy.expand_dims(salprof,axis=1)
      temprof=numpy.expand_dims(temprof,axis=1)
   
   # Section
   else  :
      section=True
      lon2d,lat2d=numpy.meshgrid(s_lon,s_lat)
      sec=gridxsec.Section([lon,lon2],[lat,lat2],lon2d,lat2d)
      i,j=sec.grid_indexes
      salprof = s_ncid["s_an"][0,:,:,:]
      temprof = t_ncid["t_an"][0,:,:,:]
      salprof = salprof[:,j,i]
      temprof = temprof[:,j,i]

   # MAsk salinity and temperature profiles
   salprof=numpy.ma.masked_where(salprof ==  s_ncid["s_an"]._FillValue,salprof)
   temprof=numpy.ma.masked_where(temprof ==  t_ncid["t_an"]._FillValue,temprof)

   # Interface values dprofi
   dprofi=numpy.zeros(nz+1)
   for k in range(dprof.size) :
      if k < dprof.size -1 :
         dprofi[k] = 0.5 * (dprof[k]+dprof[k+1])
      else :
         dprofi[k] = dprof[k]


   # MAx depth in data. TODO. Get from profiles
   maxd=numpy.zeros(salprof[0,:].shape)
   for k in range(nz) :
      maxd[~salprof[k,:].mask] = dprof[k]
   logger.info("MAx depth in data is %d m"%numpy.max(maxd))

   # Open blkdat.input
   bp=modeltools.hycom.BlkdatParser(blkdatfile)
   dp0k = bp["dp0k"]
   ds0k = bp["ds0k"]
   sigma = bp["sigma"]
   kdm=bp["kdm"]
   nhybrd=bp["nhybrd"]
   nsigma=bp["nsigma"]
   thflag=bp["thflag"]
   kapref=bp["kapref"]
   isotop=bp["isotop"]
   eqstate   = modeltools.hycom.Sigma(thflag)
   #eqstate   = modeltools.hycom.Sigma12Term(thflag)  #Test 12 term sigma
   #eqstate   = modeltools.hycom.Sigma17Term(thflag)  #Test 17 term sigma
   sigprof = eqstate.sig(temprof,salprof)
   sigprof=numpy.ma.masked_where(salprof.mask,sigprof)
  

   # Thermobaricity 
   if kapref <> 0 :
      if thflag == 2 :
         mykappa=modeltools.hycom.Kappa(kapref,2000.0e4)
      elif thflag == 0  :
         mykappa=modeltools.hycom.Kappa(kapref,0.)
      else :
         raise ValueError,"Unknown value of thflag"
      tmp = mykappa.kappaf(salprof.transpose(),temprof.transpose(),sigprof.transpose(),dprof*9806)
      #print tmp.shape,sigprof.shape
      sigprof_tbar = sigprof - tmp.transpose()



   # Min thickness interface values
   intf,masks=bp.intf_min_profile(numpy.array(maxd))
   dp0=intf[:,1:]-intf[:,:-1]
   intfmid=(intf[:,1:]+intf[:,:-1])*.5
   #print dp0[numpy.argmax(maxd),:],dp0.shape


   # We have interface values, now use designated layer sigma values and actual sigma
   # values to find an approximate vertical coordinate setup. 
   newintf,intsig=modeltools.tools.isopycnal_coordinate_layers(dprofi,sigprof,dp0,numpy.array(sigma),isotop)


   # TODO: move into function
#   # Loop over output layers
#   newintf=numpy.zeros(intf.shape)
#   intsig=numpy.zeros(dp0.shape)
#   for k in range(kdm) :
#
#      # Target layer upper interface
#      upint=newintf[:,k]
#      #print upint
#
#      # Mix water over integration range
#      dpsum=numpy.zeros(upint.shape)
#      sgsum=numpy.zeros(upint.shape)
#      sg   =numpy.zeros(upint.shape)
#      for k2 in range(nz) :
#
#
#         # Range of this layer
#         upint2=dprofi[k2]
#         lwint2=dprofi[k2+1]
#
#         # Part of this layer in target layer
#         u = numpy.maximum(upint2,upint)
#         l = lwint2
#         dp=numpy.maximum(0.,l-u)
#
#         #Mask where dpsum > 0  
#         Imask = dpsum > 0.
#
#         ####################  dpsum > 0  and summed layer density < target density ###########
#
#         # Integrated value of sigma up to this point
#         sg[Imask]=sgsum[Imask]/dpsum[Imask]
#
#         # Mask where target layer heavier than integrated value => Ok to add more layers
#         Jmask = numpy.logical_and(Imask,sigma[k] > sg)
#         Jmask = numpy.logical_and(Jmask,~salprof.mask[k2,:])
#
#         # Fraction of layer to be added
#         dpfrac = (dpsum[Jmask]*sg[Jmask] - dpsum[Jmask]*sigma[k]) / (sigma[k] - sigprof[k2,Jmask])
#
#         # Can not add more than dp!
#         dpfrac = numpy.minimum(dpfrac,dp[Jmask]) #
#
#         # Dpfrac can be negative if sigprof[k2,Jmask] < sigma[k] ( must mix "negatively"). 
#         # In that case add entire layer. Perhaps next heavy layer can tip the scale ...
#         dpfrac=numpy.where(dpfrac < 0.,dp[Jmask],dpfrac)
#
#         # Update values 
#         dpsum[Jmask]=dpsum[Jmask]+dpfrac
#         sgsum[Jmask]=sgsum[Jmask]+sigprof[k2,Jmask]*dpfrac
#
#         # NB: No need to treat summed layer density > target density, since mixing wont 
#         # enable us to reach target density. (sg can only increase as we sum deeper)
#
#
#         ####################  dpsum == 0  and sigma[k] > sigprof[k2,:] ######################
#
#         #sg = sigprof[k2,k2]
#
#         # target layer heavier than layer value. Ok to add more layers to mix downward
#         Jmask = numpy.logical_and(~Imask,sigma[k] >  sigprof[k2,:])
#         Jmask = numpy.logical_and(Jmask,~salprof.mask[k2,:])
#         dpsum[Jmask]=dpsum[Jmask]+dp[Jmask]
#         sgsum[Jmask]=sgsum[Jmask]+sigprof[k2,Jmask]*dp[Jmask]
#
#         # No need to treat new layer density > target layer density, since mixing wont enable 
#         # us to reach target density
#
#         #################### ######################################### ######################
#
#
#
#
#      # end k2 loop
#      #print dpsum
#
#      # Make sure dpsum adheres to minimum layer thickness
#      dpsum=numpy.maximum(dpsum,dp0[:,k])
#      
#      # Adjust layer interface
#      newintf[:,k+1] = newintf[:,k] + dpsum
#
#      # MAke sure lowest interface is above sea floor
#      newintf[:,k+1] = numpy.minimum(newintf[:,k+1],maxd)
#      
#      # Effective layer thickness
#      dpsum = newintf[:,k+1] - newintf[:,k]
#
#      # Estimated sigma
#      msk=dpsum>0.
#      intsig[msk,k]=sgsum[msk]/dpsum[msk]
#
#      # TODO: Keep track of final layer density for all layers (not just isopycnals)
#   # End k loop
#
#   # Make sure bottom layer reaches sea floor
#   newintf[:,kdm] = numpy.maximum(newintf[:,kdm],maxd)
#   newintfmid=(newintf[:,1:]+newintf[:,:-1])*.5



#   # Smooth eine bitchen...
#   if maxd.size > 5 :
#      w=numpy.ones(5)
#      for k in range(kdm+1) :
#         newintf[:,k]=numpy.convolve(w/w.sum(),newintf[:,k],mode="same")
#         # But make sure we adhere to minimum layer thicknesses and max depth
#         if k>0 :
#            newintf[:,k]=numpy.maximum( newintf[:,k-1]+dp0[:,k-1], newintf[:,k])
#         newintf[:,k]=numpy.minimum( newintf[:,k],maxd[:])


   # TODO: Integrate over temp and sal profiles for plotting these





   # Find segments where its ok to place layer legend
   import scipy.ndimage.measurements
   lab,num_features =scipy.ndimage.measurements.label(maxd > numpy.max(maxd)*.25)
   feat_count = {}
   for i in range(num_features): feat_count[i+1] = numpy.count_nonzero(lab==i+1)
   tmp = sorted(feat_count.items(), key=lambda x:x[1],reverse=True) # Order by feature count
   main_feature=tmp[0][0]
   ft_ind=[]
   for i_enum,i in enumerate(tmp) :
      #logger.info( "Feature %03d: %d cells" %(i[0],i[1]))
      tmp=numpy.where(lab==i[0])
      ft_ind.append(int(tmp[0].mean()))
   #print ft_ind
   #tmp=numpy.where(lab==main_feature)
   #print tmp
   #sec_xpind = int (tmp[0].mean())
   #print sec_xpind
   
   cols="b"
   colt="r"
   cold="m"
   
   # Plot vertical profile and layers. Two cases: section or not
   if section :
      f,ax = plt.subplots(1,figsize=(10,5))
      ax.set_title("sigma-%d layers from lon=%6.2f, lat=%6.2f to lon=%6.2f, lat=%6.2f"%(thflag,lon,lat,lon2,lat2))
      x = sec.distance / 1000.
      ymult=1.
      xlim=[x.min(),x.max()]
   else :
      f,ax = plt.subplots(1)
      ax.set_title("Sigma-%d profile and layers at lon=%6.2f, lat=%6.2f"%(thflag,lon,lat))
      x = numpy.array([-100,100])
      ymult=numpy.ones((2))
      ax.plot(sigprof,-dprof,lw=2,color=cold,label="Sigma-%d"%thflag)
      xlim=ax.get_xlim()
   #print ax.get_position()
   ax.set_position([.1,.1,.7,.7])


   ylim=[-numpy.max(maxd)+10,0.]
   for k in range(kdm+1) :

      if k<kdm :
         ax.plot(x,-newintf[:,k]*ymult,lw=.1,color=".5",label="layer %d: %6.3f"%(k+1,sigma[k]))
         if section :
            ##ind = x.size/2
            #ind = numpy.argmax(maxd)
            #xp = x[ind]
            #yp = -0.5*(newintf[ind,k] + newintf[ind,k+1])
            ##ax.text(xp,yp,"layer %d: %4.2f"%(k+1,sigma[k]),verticalalignment="center",horizontalalignment="left",fontsize=6)
            #ax.text(xp,yp,"%d"%(k+1),verticalalignment="center",horizontalalignment="left",fontsize=6)
            for ind in ft_ind :
               xp = x[ind]
               yp = -0.5*(newintf[ind,k] + newintf[ind,k+1])
               #ax.text(xp,yp,"layer %d: %4.2f"%(k+1,sigma[k]),verticalalignment="center",horizontalalignment="left",fontsize=6)
               ax.text(xp,yp,"%d"%(k+1),verticalalignment="center",horizontalalignment="left",fontsize=6)
            # TODO
            #ind = x.size-x.size/6
            #xp = x[ind]
            #yp = -0.5*(newintf[ind,k] + newintf[ind,k+1])
            #ax.annotate("layer %d: %4.2f"%(k+1,sigma[k]),
            #    xy=(xp, yp), xycoords='data',
            #    xytext=(1.01,0.0-k*.03), textcoords='axes fraction',
            #    arrowprops=dict(facecolor='b', shrink=0.01,alpha=.5,ec="none",width=1,headwidth=3,frac=.03),
            #    horizontalalignment='left', verticalalignment='center',
            #    fontsize=6
            #    )
         else :
            xp = numpy.sum(xlim)/2.
            yp = -0.5*(newintf[0,k] + newintf[0,k+1])
            ax.text(xp,yp,"layer %d: %4.2f"%(k+1,sigma[k]),verticalalignment="center",horizontalalignment="left",fontsize=6)

         isopyc = numpy.abs(intsig[:,k]-sigma[k])<1e-3
         ax.fill_between(x,-newintf[:,k+1]*ymult,-newintf[:,k]*ymult,color="g",alpha=".5",where=isopyc*ymult)
      else :
         # Plot sea floor
         #ax.plot(x,-newintf[:,k]*ymult,lw=2,color="k",label="Sea floor")
         ax.fill_between(x,-newintf[:,k]*ymult,-10*numpy.ones(x.shape)*maxd.max(),color=".5")

   # Some dummy plots for the legend
   ax.plot([-100,100],[1e5,1e5],lw=.1,color=".5",label="Layer interface")
#   ax.fill_between([-100,100],[1e5,1e5],[2e5,2e5],color="g",alpha=".5",label="Isopycnal layer") #NB: Label doesnt work



   ax.grid(False)
   ax.set_xlim(xlim)
   ax.set_ylim(ylim)
   ax.set_ylabel("Depth")
   if section : 
      ax.set_xlabel("Distance along section[km]")
   else :
      ax.set_xlabel("Density[sigma-%d]"%thflag,color=cold)
   for t in ax.get_xticklabels() :
      if not section : t.set_color(cold)
      t.set_size(8)
      t.set_rotation(-45)
   #ax.legend(loc="lower left")
   leg=ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,labelspacing=-.5)
   ltext=leg.get_texts()
   plt.setp(ltext,fontsize=4)


   
   fname="layers_%s.png"%sinfo
   logger.info("Layers in %s"%fname)
   plt.gcf().savefig(fname,dpi=dpi)
   #
   ax.set_ylim(-750,0)
   fname="layers750_%s.png"%sinfo
   logger.info("Layers for top 750m in %s"%fname)
   plt.gcf().savefig("layers750_%s.png"%sinfo,dpi=dpi)
   #
   ax.set_ylim(-100,0)
   fname="layers100_%s.png"%sinfo
   logger.info("Layers for top 100m in %s"%fname)
   plt.gcf().savefig("layers100_%s.png"%sinfo,dpi=dpi)


   # Plot density of original data
   f,ax = plt.subplots(1,figsize=(10,5))
   P=ax.pcolormesh(x,-dprof,sigprof,cmap="Paired")
   CS=ax.contour(x,-dprof,sigprof,sigma)
   ax.clabel(CS, inline=1, fontsize=4)
   f.colorbar(P,ax=ax)
   f.savefig("dens_data_%s.png"%sinfo,dpi=dpi)

   # Plot density of original data
   if kapref > 0 :
      f,ax = plt.subplots(1,figsize=(10,5))
      P=ax.pcolormesh(x,-dprof,sigprof_tbar,cmap="Paired")
      CS=ax.contour(x,-dprof,sigprof_tbar,sigma)
      ax.clabel(CS, inline=1, fontsize=4)
      f.colorbar(P,ax=ax)
      f.savefig("dens_data_kappa%d_%s.png"%(kapref,sinfo),dpi=dpi)

   # Plot vertical density gradient of original data
   #print sigprof.shape,dprof.shape
   dsigprof = numpy.zeros(sigprof.shape)
   for k in range(nz) :
      #print k,nz
      if k<nz-1 :
         dsigprof[k,:] = (sigprof[k+1,:] - sigprof[k,:])/dprof[k]
      else :
         dsigprof[k,:] = dsigprof[k-1,:] 
   dsigprof = numpy.ma.masked_where(numpy.abs(dsigprof)>1e10,dsigprof)
   dsigprof=dsigprof * 100 # DS per 100 meters
   f,ax = plt.subplots(1,figsize=(10,5))
   P=ax.pcolormesh(x,-dprof,dsigprof,
         norm=matplotlib.colors.LogNorm(vmin=0.001,vmax=10),
         cmap="Paired")
   CS=ax.contour(x,-dprof,sigprof,sigma)
   ax.clabel(CS, inline=1, fontsize=4)
   ax.set_title("Delta rho per 100 meter in the vertical")
   f.colorbar(P,ax=ax)
   f.savefig("densgrad_data_%s.png"%sinfo,dpi=dpi)

   # Plot Buoyancy frequency
   Nsq=dsigprof/100. * 9.81/1000.
   f,ax = plt.subplots(1,figsize=(10,5))
   P=ax.pcolormesh(x,-dprof,Nsq,
         norm=matplotlib.colors.LogNorm(vmin=1e-9,vmax=1e-4),
         cmap="Paired")
   #CS=ax.contour(x,-dprof,Nsq,sigma)
   ax.clabel(CS, inline=1, fontsize=4)
   f.colorbar(P,ax=ax)
   ax.set_title("N^2")
   f.savefig("Nsq_data_%s.png"%sinfo,dpi=dpi)




   # TODO: Plot Temperature and Salinity sections





if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('blkdatfile' , help="")
   parser.add_argument('saltfile' , help="")
   parser.add_argument('lon'      , type=float) 
   parser.add_argument('lat'      , type=float) 
   parser.add_argument('--lon2'      , type=float) 
   parser.add_argument('--lat2'      , type=float) 
   parser.add_argument('--sectionid'      , type=str,default="") 



   args = parser.parse_args()
   #print args
   main(args.blkdatfile,args.saltfile,args.lon,args.lat,lon2=args.lon2,lat2=args.lat2,sectionid=args.sectionid)



