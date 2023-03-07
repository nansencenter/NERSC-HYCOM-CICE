import numpy as N
import numpy.ma as ma
##import Scientific.IO.NetCDF as S
#from netCDF4 import Dataset
import netCDF4 as nc 
import matplotlib.pyplot as plt
#import mpl_toolkits.basemap as bm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

# reading the observation file
with open('Obs_drift.txt','r') as file:
  Fline=file.readline()
  ds=nc.Dataset(Fline[:-1],'r')
  dx=ds['dX'][0,:,:]; dy=ds['dY'][0,:,:]
  lat1=ds['lat'][:,:]; lon1=ds['lon'][:,:]; derr=ds['uncert_dX_and_dY'][0,:,:]
  dsp=N.sqrt(dx*dx+dy*dy)

# reading the model result
filename0='Mod_drift000.nc'
ds=nc.Dataset(filename0,'r')
dx0=ds['dX-ice'][:,:]; dy0=ds['dY-ice'][:,:]
lat=ds['lat'][:,:]; lon=ds['lon'][:,:]; dflg=ds['ObsFlg'][:,:]
dsp0=N.sqrt(dx0*dx0+dy0*dy0)

Cname='jet'

SMALL_SIZE =10 
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

cmap=plt.get_cmap(Cname)
#str0='SIT_Fore\n 19785 \n '+str(N.mean(tmpave))
#str0="Mod\n drift \n %0.2f (m)"%N.mean(tmpave)
lonT=100
latT=100
for record in range(6):
   if record==1:
      Expnam="Mod"
      fld0=ma.masked_where(dflg<30,dx0)
      str0="dx Mod \n (km/d)"
      Fout=Expnam+"%03d.png"%record
   elif record==2:
      Expnam="Mod"
      fld0=ma.masked_where(dflg<30,dy0)
      str0="dy Mod \n (km/d)"
      Fout=Expnam+"%03d.png"%record
   elif record==0:
      str0="drift \n (km/d)"
      Expnam="Mod_ds"
      fld0=ma.masked_where(dflg<30,N.sqrt(dx0*dx0+dy0*dy0))
      Fout=Expnam+"%03d.png"%record
   elif record==3:
      str0="drift \n (km/d)"
      Expnam="Obs_ds"
      fld0=ma.masked_where(dflg<30,N.sqrt(dx*dx+dy*dy))
      Fout=Expnam+"%03d.png"%(record-3)
   elif record==4:
      Expnam="Obs"
      fld0=ma.masked_where(dflg<30,dx)
      str0="dx Osisaf \n (km/d)"
      Fout=Expnam+"%03d.png"%(record-3)
   elif record==5:
      Expnam="Obs"
      fld0=ma.masked_where(dflg<30,dy)
      str0="dy Osisaf \n (km/d)"
      Fout=Expnam+"%03d.png"%(record-3)
   
   figure = plt.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   if record==0 or record==3:
      P=ax.pcolormesh(fld0,cmap=cmap,vmin=0,vmax=15)
   else:
      P=ax.pcolormesh(fld0,cmap=cmap,vmin=-10,vmax=10)
   clb=ax.figure.colorbar(P)
   plt.text(lonT,latT,str0,fontsize=14,color='w');
   ax.set_facecolor("gray")
   plt.xlabel('J-index',fontsize=14)
   plt.ylabel('I-index',fontsize=14)
   figure.canvas.print_figure(Fout)

for record in range(3):
   if record==0:
      fld0=dsp0   # model
      fld1=dsp    # observation
      str0="Obs drift (km/d)"
      str1="Mod drift (km/d)"
   elif record==1:
      fld0=dx0   # model
      fld1=dx    # observation
      str0="Obs dx (km/d)"
      str1="Mod dx (km/d)"
   elif record==2:
      fld0=dy0   # model
      fld1=dy    # observation
      str0="Obs dy (km/d)"
      str1="Mod dy (km/d)"

   Fout="DriftScatter_"+"%03d.png"%record
   figure = plt.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   P0=ax.scatter(fld0,fld1,c="blue")
   ax.grid()
   plt.xlabel(str0,fontsize=14)
   plt.ylabel(str1,fontsize=14)
   figure.canvas.print_figure(Fout)
  
