import modeltools.hycom
import matplotlib
import cmocean
import abfile.abfile as abf
import numpy as np
import datetime
import re
import scipy.interpolate
import os.path
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, DayLocator 


"""
This module contains useful functions, and always we can add more here

it can be imported as: import mod_hyc2plot
"""



def replaceZeroes(data):
 min_nonzero = np.min(data[np.nonzero(data)])
 data[data == 0] = min_nonzero
 return data

### def sph2gr(long,latt,xpole,ypole,dx,dy,ylong):
###     """
###     Convert from longitude,latitude (given as vectors) to grid coordinates
###     defined by xpole,ypole,dx,dy,ylong. (dx,dy has unit km)
###     """
###     long=nx.array(long)
###     latt=nx.array(latt)
###     rad = 3.14159265 / 180.
###     rearth = 6370.
###     phinul = 60. * rad
### 
###     #  convert degrees to radians.
###     lamb = long  * rad
###     phi  = latt  * rad
###     alpha= ylong * rad
### 
###     #  the computations.
###     r = rearth * nx.cos(phi) * (1 + nx.sin(phinul)) \
###         /(1 + nx.sin(phi))
###     x = xpole + r * nx.sin(lamb-alpha) / dx
###     y = ypole - r * nx.cos(lamb-alpha) / dy
### 
###     return x,y
### 
### #function [tke,mke,eke]=calc_kinetic(uin,vin)

def calc_kinetic(uin,vin):

   """
   % function [tke,mke,eke]=calc_kinetic(uin,vin)
   % 
   % routine to calculate total kinetic energy (tke), kinetic energy of the mean flow (mke),
   % and eddy kinetic energy (eke).
   %
   % tke can be split into mke and eke as follows:
   % tke = mke + eke
   %
   % mke reflects the energy from the mean flow and eke reflects the energy from the
   % mesoscale variability (eddies).
   %
   % applying Reynold's averaging: U = Ubar + U'
   % where Ubar is the time mean of U, and U' the mesoscale deviations (U' = U - Ubar)
   % the same applies for V
   % 
   % mke = (Ubar*Ubar + Vbar*Vbar)/2
   %
   % eke = (U'U'bar + V'V'bar)/2
   % because U'U'bar (the time mean of the velocity correlations) is difficult to calculate
   % first calculate Ubar (time average U) and UUbar (time average U*U) and then
   % U'U'bar = UUbar - Ubar*Ubar (same for V)
   % 
   """
   print(">>>uin.shape, vin.shape="),print(np.shape(uin)),print(np.shape(vin)) 
   # calculate Ubar / Vbar
   ubar=np.nanmean(uin,axis=0)
   vbar=np.nanmean(vin,axis=0)
   # calculate mke = (Ubar*Ubar + Vbar*Vbar)/2
   mke=(ubar*ubar+vbar*vbar)/2
   
   # calculate UUbar and VVbar
   uubar=np.nanmean(uin*uin,axis=0)
   vvbar=np.nanmean(vin*vin,axis=0)
   
   # calculate eke =  (U'U'bar + V'V'bar)/2
   # where U'U'bar = UUbar - Ubar*Ubar
   # and V'V'bar = VVbar - Vbar*Vbar
   uprimeuprimebar=uubar-ubar*ubar
   vprimevprimebar=vvbar-vbar*vbar
   eke=(uprimeuprimebar+vprimevprimebar)/2
   
   # calculate tke
   tke=mke+eke
   
   return tke,mke,eke


def runningMean(x, N):
    #running mean
    #y = np.zeros((len(x),))
    #y = np.zeros(shape=x.shape)
    print("Calculating running mean over"), print(x.shape), print("mean window="), print(N) 
    #for ctr in range(x.shape[0]):
    #     y[ctr] = np.nansum(x[ctr:(ctr+N)],axis=0)
    #return y/N
    #assert N%2==1
    y=np.zeros(shape=x.shape)
    y[:]=np.nan
    for ii in range(x.shape[1]):
       for jj in range(x.shape[2]):
          #y[:,ii,jj]=np.convolve(x[:,ii,jj], np.ones((N,dtype='float'))/N, mode='same')
          if ~np.isnan(x[0,ii,jj]):
             y[:,ii,jj]=np.convolve(x[:,ii,jj], np.ones(N,dtype='float'), mode='same')/ \
                np.convolve(np.ones(x.shape[0]), np.ones((N)), mode='same')
             #print "Numer =", np.convolve(np.ones(x.shape[0]), np.ones((N)), mode='same')
             #sci.convolve(a,np.ones(n,dtype='float'), 'same')/sci.convolve(np.ones(len(a)),np.ones(n), 'same')
             #y[:,ii,jj]=pd.Series(x).rolling(window=N,axis=0).mean().iloc[N-1:].values
    #y=np.convolve(x[:], np.ones((N,))/N, mode='same')
    return y


def spatiomean(fldin,regi_mask):
  #  fldin=[]
    ab = abf.ABFileGrid("regional.grid","r")
    pplon=ab.read_field("plon")
    pplat=ab.read_field("plat")
    scppx=ab.read_field("scpx")
    scppy=ab.read_field("scpy")
    abdpth = abf.ABFileBathy('regional.depth',"r",idm=ab.idm,jdm=ab.jdm)
    mdpth=abdpth.read_field('depth')
    maskdd=mdpth.data
    maskdd[maskdd>1e29]=np.nan
    #fldin[fldin>1e29]=np.nan
    #scppx[np.isnan(maskdd)]=np.nan
    #scppy[np.isnan(maskdd)]=np.nan
###     # mask for specific region
###     #Nordic=False
###     if Nordmask:
###        print 'Compute for Nordic------>>>>>>-'
###        maskdd[pplat>80]=np.nan
###        maskdd[pplat<55]=np.nan
###        maskdd[pplon>60]=np.nan
###        maskdd[pplon<-60]=np.nan
###        #Norid
###        #fldin[np.isnan(maskdd)]=np.nan
###        #scppx[np.isnan(maskdd)]=np.nan
###        #scppy[np.isnan(maskdd)]=np.nan
###     #--
    numer=fldin*scppx*scppy
    denum=scppx*scppy
    numer[np.isnan(regi_mask)]=np.nan
    denum[np.isnan(regi_mask)]=np.nan
    fldin_avg=np.nansum(numer)/np.nansum(denum)
    print('np.nansum(numer)='), print(np.nansum(numer))
    print('np.nansum(denum)='), print(np.nansum(denum))
    print('fldin_avg='), print(fldin_avg)
    #direct mean
    return fldin_avg



# Sigma function. NB: This is the standard 7-term approach
def sig(T,S) :
   # Sigma function. NB: This is the standard 7-term approach
   # --- coefficients for sigma-0 (based on Brydon & Sun fit)
   C1=-1.36471E-01
   C2= 4.68181E-02
   C3= 8.07004E-01
   C4=-7.45353E-03
   C5=-2.94418E-03
   C6= 3.43570E-05
   C7= 3.48658E-05
   #pref=0.       # reference pressure, Pascals
   # --- sigma-theta as a function of temp (deg c) and salinity (mil)
   # --- (friedrich-levitus 3rd degree polynomial fit)
   dens=(C1+C3*S+np.multiply(T,(C2+C5*S+np.multiply(T,(C4+C7*S+C6*T)))))
   print("dens min, max="),print(dens.min()),print(dens.max())
   #return (C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))
   return dens


def plot_hycom_func1D(tid,dat1,fig_tit,labl1='dat1',dat2=None,\
                  labl2='dat2',unit=None):
   #Plot hycom field
   figure, ax =plt.subplots()
   days = DayLocator()   # every year
   years = YearLocator()   # every year
   months = MonthLocator()  # every month
   yearsFmt = DateFormatter('%Y')
   daysFmt = DateFormatter('%Y-%m-%d')
   #ax=figure.add_subplot(111)
   ax.plot(tid,dat1,'-',color='g',label=labl1)
   if dat2 is not None:
      ax.plot(tid,dat2,'--',color='b',label=labl2)
   #ax.grid(True)
   #ax.xaxis.set_major_locator(years)
   #ax.xaxis.set_major_formatter(yearsFmt)
   ax.xaxis.set_major_locator(days)
   ax.xaxis.set_major_formatter(daysFmt)
   #ax.xaxis.set_minor_locator(months)
   # format the coords message box
   ax.autoscale_view()
   # format the coords message box
   def price(x):
       return '$%1.2f' % x
   ax.fmt_xdata = DateFormatter('%Y-%m-%d')
   #ax.fmt_ydata = price
   figure.autofmt_xdate()
   #ax.locator_params(axis='x', nbins=10)
   print('len(tid)='), print(len(tid))
   if len(tid) > 100:
      every_nth = 30
   else:
      every_nth = 5
   for n, label in enumerate(ax.xaxis.get_ticklabels()):
      if n % every_nth != 0:
         label.set_visible(False)
   legend =plt.legend(loc='upper left',fontsize=8)
   fnamepng='Fig1D_'+fig_tit
   print("output in  %s"%fnamepng)
   figure.canvas.print_figure(fnamepng,bbox_inches='tight',dpi=180)




def plot_hycom_func(lon,lat,fld,fig_tit,uu=None,vv=None, \
                  cmap=None,clim=None, fldname=None,unit=None):
   #Plot hycom field
   #bm=None
   from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
   from mpl_toolkits.basemap import Basemap
   LinDic=cmap_dict('sawtooth_fc100.txt')
   if cmap is None:
      cmap= matplotlib.colors.LinearSegmentedColormap('my_colormap',LinDic)
   else:
      cmap=plt.get_cmap("jet")

   bm = Basemap(width=7400000,height=7400000, \
         resolution='i',projection='stere',\
         lat_ts=70,lat_0=85,lon_0=-40.)
   x,y=bm(lon,lat)   
   J,I=np.meshgrid(np.arange(fld.shape[0]),np.arange(fld.shape[1])) 
   #
   print("---------min,max of data="), print(fld.min()), print(fld.max())
   figure = plt.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P=bm.pcolormesh(x[J,I],y[J,I],fld[J,I],cmap=cmap)
   if 'temp' in fig_tit:
      P1=bm.contour(x[J,I],y[J,I],fld[J,I],levels=[-1.,1,4.0,8], \
             colors=('w',),linestyles=('-',),linewidths=(1.5,))
      plt.clabel(P1, fmt = '%2.1d', colors = 'w', fontsize=10) #contour line labels
   bm.drawcoastlines(linewidth=0.05)
   bm.drawcountries(linewidth=0.05)
   bm.fillcontinents(color='.8',lake_color='white')
   bm.drawparallels(np.arange(-80.,81.,40.),linewidth=0.2)
   bm.drawmeridians(np.arange(-180.,181.,40.),linewidth=0.2)
   bm.drawmapboundary(linewidth=0.2) #fill_color='aqua')
   ##
   if uu is not None:
      skip=10
      print("ploting quiver .......>>> ")
      I2=I[::skip,::skip]
      J2=J[::skip,::skip]
      ax.quiver(x[J2,I2],y[J2,I2],uu[J2,I2],vv[J2,I2])
               

   # Print figure.
   aspect = 40
   pad_fraction = 0.25
   divider = make_axes_locatable(ax)
   width = axes_size.AxesY(ax, aspect=1./aspect)
   pad = axes_size.Fraction(pad_fraction, width)
   cax = divider.append_axes("right", size=width, pad=pad)
   cb=ax.figure.colorbar(P,cax=cax,extend='both')
   #P.set_clim([fld.min(),fld.max()])
   if clim is not None : P.set_clim(clim)
   if fldname is not None:
      ax.set_title(fig_tit+' :'+fldname+'['+unit+']')
   else:
      ax.set_title(fig_tit)
   # Print figure.
   fnamepng='Fig_'+fig_tit
   print("output in  %s"%fnamepng)
   dpi=180
   figure.canvas.print_figure(fnamepng,bbox_inches='tight',dpi=dpi)
   plt.close(figure)


def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    # https://ocefpaf.github.io/python4oceanographers/blog/2014/03/10/gearth/
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    fig = plt.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    ax.set_axis_off()
    return fig, ax


def interpolate_to_latlon(lon,lat,data,res=0.1) :
   lon2 = np.mod(lon+360.,360.)
   # New grid
   minlon=np.floor((np.min(lon2)/res))*res
   minlat=max(-90.,np.floor((np.min(lat)/res))*res)
   maxlon=np.ceil((np.max(lon2)/res))*res
   maxlat=min(90.,np.ceil((np.max(lat)/res))*res)
   #maxlat=90.
   lo1d = np.arange(minlon,maxlon+res,res)
   la1d = np.arange(minlat,maxlat,res)
   lo2d,la2d=np.meshgrid(lo1d,la1d)
   print(minlon,maxlon,minlat,maxlat)

   if os.path.exists("grid.info") :
      import modelgrid

      grid=modelgrid.ConformalGrid.init_from_file("grid.info")
      map=grid.mapping

      # Index into model data, using grid info
      I,J=map.ll2gind(la2d,lo2d)

      # Location of model p-cell corner 
      I=I-0.5
      J=J-0.5

      # Mask out points not on grid
      K=J<data.shape[0]-1
      K=np.logical_and(K,I<data.shape[1]-1)
      K=np.logical_and(K,J>=0)
      K=np.logical_and(K,I>=0)

      # Pivot point 
      Ii=I.astype('i')
      Ji=J.astype('i')

      # Takes into account data mask
      tmp =np.logical_and(K[K],~data.mask[Ji[K],Ii[K]])
      K[K]=tmp
      
      tmp=data[Ji[K],Ii[K]]
      a,b=np.where(K) 
      z=np.zeros(K.shape)
      z[a,b] = tmp
      z=np.ma.masked_where(~K,z)

   # Brute force ...
   else  :
      K=np.where(~data.mask)
      z=scipy.interpolate.griddata( (lon2[K],lat[K]),data[K],(lo2d,la2d),'cubic')
      z=np.ma.masked_invalid(z)
      z2=scipy.interpolate.griddata( (lon2.flatten(),lat.flatten()),data.mask.flatten(),(lo2d,la2d),'nearest')
      z2=z2>.1
      z=np.ma.masked_where(z2,z)

   return lo2d,la2d,z


# Simple routine to create a kml file from field

def to_kml(lon,lat,data,fieldname,cmap,res=0.1,clim=None) :
   import simplekml
   
   lo2d,la2d,z= interpolate_to_latlon(lon,lat,data,res=res)
   urcrnrlon = lo2d[-1,-1]
   urcrnrlat = la2d[-1,-1]
   llcrnrlon = lo2d[0,0]
   llcrnrlat = la2d[0,0]

   fig,ax= gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024)
   figname="overlay.png"
   P=ax.pcolormesh(lo2d,la2d,z,cmap=cmap)
   if clim is not None : P.set_clim(clim)
   fig.canvas.print_figure(figname,Transparent=True)

   kw={}
   kml=simplekml.Kml()
   draworder = 0
   draworder += 1
   ground = kml.newgroundoverlay(name='GroundOverlay')
   ground.draworder = draworder
   ground.visibility = kw.pop('visibility', 1)
   ground.name = kw.pop('name', fieldname)
   ground.color = kw.pop('color', 'ddffffff') # First hex gives transparency
   ground.atomauthor = kw.pop('author', 'NERSC')
   ground.latlonbox.rotation = kw.pop('rotation', 0)
   ground.description = kw.pop('description', 'Matplotlib figure')
   ground.gxaltitudemode = kw.pop('gxaltitudemode', 'clampToGround')
   ground.icon.href = figname
   ground.latlonbox.east = llcrnrlon
   ground.latlonbox.south = llcrnrlat
   ground.latlonbox.north = urcrnrlat
   ground.latlonbox.west = urcrnrlon
   kmzfile="overlay.kmz"
   kml.savekmz(kmzfile)




def cmap_dict(cmapfile):
   #LinL = np.loadtxt('sawtooth_0-1.txt')
   LinL = np.loadtxt(cmapfile)
   #LinL = np.loadtxt('sawtooth_fc100.txt')
   
   b3=LinL[:,2] # value of blue at sample n
   b2=LinL[:,2] # value of blue at sample n
   b1=np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1
   
   # setting up columns for list
   g3=LinL[:,1]
   g2=LinL[:,1]
   g1=np.linspace(0,1,len(g2))
   
   r3=LinL[:,0]
   r2=LinL[:,0]
   r1=np.linspace(0,1,len(r2))
   
   # creating list
   R=zip(r1,r2,r3)
   G=zip(g1,g2,g3)
   B=zip(b1,b2,b3)
   
   # transposing list
   RGB=zip(R,G,B)
   rgb=zip(*RGB)
   # print rgb
   
   # creating dictionary
   k=['red', 'green', 'blue']
   LinearL=dict(zip(k,rgb)) # makes a dictionary from 2 lists

   return LinearL
   

#rgb=np.array([[1.0000, 1.0000, 1.0000],
#        [0.9763, 0.9235, 0.9955], 
#             [0.9567, 0.8603, 0.9918]])

