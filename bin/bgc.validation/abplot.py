import modeltools.hycom
import abfile
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from netCDF4 import Dataset as NetCDFFile
#import mycolors as mycmaps
import nclcmaps
import cmocean
import warnings
warnings.filterwarnings('ignore')
import gridxsec
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

'''
Usage:

import abplot
ab = abplot.openfile('/cluster/work/users/cagyum/TP5a0.06/expt_01.0/data/archm.2013_014_12.b')
no3 = abplot.getvarib(ab,"no3",1)
'''

def getgrid(*args, **kwargs):
  region = kwargs.get('region',None)
  region = 'TP5a0.06' if region is None else region
  abgrid = abfile.ABFileGrid("/cluster/work/users/cagyum/"+region+"/topo/regional.grid","r")
  plon=abgrid.read_field("plon")
  plat=abgrid.read_field("plat")
  jdm,idm=plon.shape
  return plon,plat

def getdepth(*args, **kwargs):
  region = kwargs.get('region',None)
  region = 'TP5a0.06' if region is None else region
  #print(region+" specific depth is retrieved")
  plon,plat = getgrid(region=region)
  abdepth = abfile.ABFileBathy("/cluster/work/users/cagyum/"+region+"/topo/regional.depth.b", \
            "r",idm=plon.shape[1],jdm=plon.shape[0])
  depthm=abdepth.read_field("depth")
#  print(region+" depth shape is: "+depthm.shape)
  return depthm

def openfile(filename) :
  '''
  ab = abplot.openfile('/cluster/work/users/cagyum/TP5a0.06/expt_01.0/data/archm.2013_014_12.b')
  '''
  i_abfile = abfile.ABFileArchv(filename,"r")
  return i_abfile

def openrestart(filename,*args, **kwargs):

  region = kwargs.get('region',None)
  region = 'TP5a0.06' if region is None else region
  plon,plat = getgrid(region=region)
  res = abfile.ABFileRestart(filename,"r",idm=plon.shape[1],jdm=plon.shape[0])
  return res

def getvarib(abfile,fieldname,fieldlevel) :
  '''
  how to:
  no3 = abplot.getvarib(ab,"no3",1)
  '''
  C2NIT = 6.625
  C2PHO = 106.0
  C2SIL = 6.625
  C2CAR = 12.01

  if str(abfile)[21:26] == 'Archv': abtype = 'archive'
  if str(abfile)[21:28] == 'Restart': abtype = 'restart'
  
  if fieldname == "no3" :
#     print('convert mgC m-3 --> mmolN m-3')
     if abtype == 'archive':
        fld = abfile.read_field('ECO_no3',fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field('ECO_no3',fieldlevel,1)
     fld = fld/C2CAR/C2NIT
  elif fieldname == "pho" :
#     print('convert mgC m-3 --> mmolP m-3')
     if abtype == 'archive':
        fld = abfile.read_field('ECO_pho',fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field('ECO_pho',fieldlevel,1)
     fld = fld/C2CAR/C2PHO
  elif fieldname == "sil" :
#     print('convert mgC m-3 --> mmolSi m-3')
     if abtype == 'archive':
        fld = abfile.read_field('ECO_sil',fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field('ECO_sil',fieldlevel,1)
     fld = fld/C2CAR/C2SIL
  elif fieldname == "det" :
#     print('convert mgC m-3 --> mmolC m-3')
     if abtype == 'archive':
        fld = abfile.read_field('ECO_det',fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field('ECO_det',fieldlevel,1)
     fld = fld/C2CAR
  elif fieldname == "chl" :
#     print('diatom + flagellate chlorophyll mgChl m-3')
     if abtype == 'archive':
        fld1 = abfile.read_field('ECO_diac',fieldlevel)
        fld2 = abfile.read_field('ECO_flac',fieldlevel)
        try:
           fld3 = abfile.read_field('ECO_cclc',fieldlevel)
           fld2 = fld2 + fld3
        except:
           pass
     if abtype == 'restart':
        fld1 = abfile.read_field('ECO_diac',fieldlevel,1)
        fld2 = abfile.read_field('ECO_flac',fieldlevel,1)
        try:
           fld3 = abfile.read_field('ECO_cclc',fieldlevel,1)
           fld2 = fld2 + fld3
        except:
           pass
     fld = fld1 + fld2
  elif fieldname == "poc" :
     if abtype == 'archive':
        fld1 = abfile.read_field('ECO_dia',fieldlevel)
        fld2 = abfile.read_field('ECO_fla',fieldlevel)
        fld3 = abfile.read_field('ECO_det',fieldlevel)
        fld4 = abfile.read_field('ECO_micr',fieldlevel)         
        try:
            fld5 = abfile.read_field('ECO_ccl',fieldlevel)
            fld4 = fld4 + fld5
        except:
           pass
     if abtype == 'restart':
        fld1 = abfile.read_field('ECO_dia',fieldlevel,1)
        fld2 = abfile.read_field('ECO_fla',fieldlevel,1)
        fld3 = abfile.read_field('ECO_det',fieldlevel,1)
        fld4 = abfile.read_field('ECO_micr',fieldlevel,1)         
        try:
            fld5 = abfile.read_field('ECO_ccl',fieldlevel,1)
            fld4 = fld4 + fld5
        except:
           pass
     fld = fld1 + fld2 + fld3 + fld4
  elif fieldname == "pbio" :
#     print('diatom + flagellate convert mgC m-3 --> mmolC m-3')
     if abtype == 'archive':
        fld1 = abfile.read_field('ECO_dia',fieldlevel)
        fld2 = abfile.read_field('ECO_fla',fieldlevel)
        try:
           fld3 = abfile.read_field('ECO_ccl',fieldlevel)
           fld2 = fld2 + fld3
        except:
           pass
     if abtype == 'restart':
        fld1 = abfile.read_field('ECO_dia',fieldlevel,1)
        fld2 = abfile.read_field('ECO_fla',fieldlevel,1)
        try:
           fld3 = abfile.read_field('ECO_ccl',fieldlevel,1)
           fld2 = fld2 + fld3
        except:
           pass
     fld = (fld1 + fld2)/C2CAR
  elif fieldname == "zbio" :
#     print('micro + meso convert mgC m-3 --> mmolC m-3')
     if abtype == 'archive':
        fld1 = abfile.read_field('ECO_micr',fieldlevel)
        fld2 = abfile.read_field('ECO_meso',fieldlevel)
     if abtype == 'restart':
        fld1 = abfile.read_field('ECO_micr',fieldlevel,1)
        fld2 = abfile.read_field('ECO_meso',fieldlevel,1)
     fld = (fld1 + fld2)/C2CAR
  elif fieldname == "sed1" :
#     print('convert mgC m-2 --> mmolC m-2')
     if abtype == 'archive':
        fld = abfile.read_field('ECO_sed1',fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field('ECO_sed1',fieldlevel,1)
     fld = fld/C2CAR
  elif fieldname == "sed2" :
#     print('convert mgC m-2 --> mmolSi m-2')
     if abtype == 'archive':
        fld = abfile.read_field('ECO_sed2',fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field('ECO_sed2',fieldlevel,1)
     fld = fld/C2CAR/C2SIL
  elif fieldname == "sed3" :
#     print('convert mgC m-2 --> mmolP m-2')
     if abtype == 'archive':
        fld = abfile.read_field('ECO_sed3',fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field('ECO_sed3',fieldlevel,1)
     fld = fld/C2CAR/C2PHO
  elif fieldname == "pp" :
     fld = integr(abfile,'ECO_prim',200.)
     fld = fld * 60. * 60. * 24.
  else :
#     print(fieldname)
     if abtype == 'archive':
        fld = abfile.read_field(fieldname,fieldlevel)
     if abtype == 'restart':
        fld = abfile.read_field(fieldname,fieldlevel,1)


  if abtype == 'restart':
     if fld.shape[0] == 760 :
        region='TP5a0.06'
     if fld.shape[0] == 460 :
        region='TZ4a0.10'
     if fld.shape[0] == 380 :
        region='TP2a0.10'      
     depthm = getdepth(region=region) 
     fld = np.ma.masked_where(depthm.mask,fld)

  return fld

def get3d(abfile,fieldname) :
  fld1 = getvarib(abfile,fieldname,1)
  dm=max(abfile.fieldlevels)
  fld = np.zeros((dm,fld1.shape[0],fld1.shape[1]))
  for k in range(dm) :
    fld1 = getvarib(abfile,fieldname,k+1)
    fld[k,:,:] = fld1
  fld = np.ma.masked_where(fld>1e10,fld)
  return fld

def depth3d(abfile):
    thk = get3d(abfile,"thknss")
    thk = thk/modeltools.hycom.onem;
    dm=max(abfile.fieldlevels)
    fld = np.zeros((thk.shape))
    fld[0,:,:] = thk[0,:,:]/2.
    for k in range(dm-1) :
        fld[k+1,:,:] = fld[k,:,:] + thk[k,:,:]/2. + thk[k+1,:,:]/2.
    fld = np.ma.masked_where(fld>1e10,fld)
    return fld

def integr(abfile,fieldname,intdepth):
    dm  = max(abfile.fieldlevels)

    k=1
    thk   = abfile.read_field('thknss',k)/modeltools.hycom.onem;
    deep  = np.copy(thk)
    field = getvarib(abfile,fieldname,k) * thk

    for k in range(2,dm):
        thk = abfile.read_field('thknss',k)/modeltools.hycom.onem;
        thickness = deep + thk
        thickness[thickness>=intdepth] = intdepth
        thickness = thickness - deep
        thickness[thickness <= 0.] = 0.

        field = getvarib(abfile,fieldname,k) * thickness + field

        deep = deep + thk

    return field

def plota(field,*args, **kwargs) :
    '''
    usage: abplot.plota(field)
    options:
            figx=10     (x-size of window, defaults to 7)
            figy=8      (y-size of window, defaults to 6)
            title='TBF' (plot title, defaults to empty)
            cmin=0      (color range minimum, defaults to field minimum)
            cmax=10     (color range maximum, defaults to field maximum)
            nlevels=20  (number of colors to use, defaults to 40)
    '''
    map = kwargs.get('map',None)
    mapNA = kwargs.get('mapNA',None)
    zoom = kwargs.get('zoom',None)
    figx = kwargs.get('figx',None)
    figx = 7 if figx is None else figx
    if zoom: figx=5
    figy = kwargs.get('figy',None)
    figy = 6 if figy is None else figy
    if zoom: figy=7
    title = kwargs.get('title',None)
    title = '' if title is None else title

    cmax = kwargs.get('cmax',None)
    cmax = field.max() if cmax is None else cmax
    cmin = kwargs.get('cmin',None)
    cmin = field.min() if cmin is None else cmin

    if field.shape[0] == 380:
       region = "TP2a0.10"
    if field.shape[0] == 760:
       region = "TP5a0.06"
    if field.shape[0] == 460:
       region = "TZ4a0.10"
    if field.shape[0] == 110:
       region = 'TP0a1.00'

    cmap = kwargs.get('cmap',None)
    if cmap is not None:
       if cmap == 'nclchl':
          cmap = nclcmaps.cmap('WhiteBlueGreenYellowRed')
       if cmap == 'spectral':
          cmap = plt.cm.get_cmap('Spectral_r')
       if cmap == 'diff':
          cmap = plt.get_cmap('RdYlBu_r')
    else:
       cmap = nclcmaps.cmap('GMT_haxby')
    depthm = getdepth(region=region)
    field = np.ma.masked_where(depthm.mask,field)

    fig=plt.figure(figsize=(figx,figy),facecolor='w',dpi=200)

    if map == True:
         import cartopy.crs as ccrs
         import cartopy.feature as cfeature
         import matplotlib.path as mpath
         
         ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-15))
         ax.set_extent([-180, 180, 54, 90], ccrs.PlateCarree())

         ax.gridlines(color='gray', linestyle=':', linewidth=0.5)
         theta = np.linspace(0, 2*np.pi, 100)
         center, radius = [0.5, 0.5], 0.6  # Increase the radius to cover more empty space
         verts = np.vstack([np.sin(theta), np.cos(theta)]).T
         circle = mpath.Path(verts * radius + center)
         ax.set_boundary(circle, transform=ax.transAxes)

         plon,plat = getgrid(region=region)
         pmesh = ax.pcolormesh(plon, plat, field, transform=ccrs.PlateCarree(), cmap=cmap)

         #ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black')
    elif mapNA == True:
         import cartopy.crs as ccrs
         import cartopy.feature as cfeature
         import matplotlib.path as mpath

         ax = plt.axes(projection=ccrs.Miller())
         ax.set_extent([-80, 30, 20, 75], crs=ccrs.PlateCarree())

         ax.gridlines(color='gray', linestyle=':', linewidth=0.5,zorder=3)
         plon,plat = getgrid(region=region)
         pmesh = ax.pcolormesh(plon, plat, field, transform=ccrs.PlateCarree(), cmap=cmap, zorder=0)
         ax.coastlines(zorder=2)
         ax.add_feature(cfeature.LAND, facecolor='lightgray',zorder=1)

    else:
         ax = fig.add_subplot(1,1,1)
         pmesh = plt.pcolormesh(field,cmap=cmap)
         ax.set_facecolor('lightgrey')
    plt.title(title, backgroundcolor='white', bbox=dict(facecolor='white', edgecolor='black'))

    pmesh.set_clim(cmin,cmax)
    if map or mapNA:
      cbar = plt.colorbar(pmesh, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
      cbar.ax.yaxis.set_ticks_position('right')
      cbar.ax.yaxis.set_label_position('right')
      cbar.ax.tick_params(labelsize=14)
      plt.tight_layout()
      cbar.ax.set_position([0.875, 0.25, 0.04, 0.5])
    else:
      cbar=plt.colorbar(pmesh)
      plt.tight_layout()
    
    plt.show()

def pltsec0(field,lon1,lat1,lon2,lat2,*args, **kwargs) :
# field is expected to be in 2D
   print(field.min(),field.max())
   region = kwargs.get('region',None)
   region = 'TP5a0.06' if region is None else region
   print('setting domain to: ', region)
   plon,plat = getgrid(region=region)

   depthm = getdepth(region=region)

   cooINDEX = abs( plat-lat1 ) + abs( plon-lon1 )
   JJ1,II1 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
   cooINDEX = abs( plat-lat2 ) + abs( plon-lon2 )
   JJ2,II2 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)

   sec = gridxsec.SectionIJSpace([II1,II2],[JJ1,JJ2],plon,plat)
   I,J=sec.grid_indexes
   dist=sec.distance
   slon=sec.longitude
   slat=sec.latitude

   xaxis = kwargs.get('xaxis',None)
   if xaxis is not None:
     if xaxis == "lon" :
      x=slon
      xlab="Longitude"
     elif xaxis == "lat" :
      x=slat
      xlab="Latitude"
   else:
      x=dist/1000.
      xlab="Distance along section[km]"     

   data2d = field[:,:]
   data2d=np.ma.filled(data2d,1e30)
   datasec = field[J,I]#data2d[J,I]

   fig=plt.figure(figsize=(13,4),facecolor='w')
   
   ax = fig.add_subplot(2,1,1)
   ax.set_position([0.075,0.12,0.6,0.8])
   ax.set_facecolor('lightgray')
   P=ax.plot(x,datasec)

   ax.set_xlabel(xlab)

   title = kwargs.get('title',None)
   title = '' if title is None else title
   plt.title(title)

   cmax = kwargs.get('cmax',None)
   cmax = field.max() if cmax is None else cmax
   cmin = kwargs.get('cmin',None)
   cmin = field.min() if cmin is None else cmin

   axm=fig.add_subplot(212)
   axm.set_position([0.72,0.12,0.225,0.85])
   data2d = field[:,:]
   data2d = np.ma.masked_where(depthm.mask,data2d)
   axm.set_facecolor('xkcd:gray')
   cmap = plt.cm.get_cmap('Spectral_r')
   pmesh = plt.pcolormesh(field,cmap=cmap)
   section = plt.plot(I,J,'r',lw=2)
   plt.clim(datasec.min(),datasec.max())
   cbaxes = fig.add_axes([0.95, 0.12, 0.02, 0.55])
   cb = ax.figure.colorbar(pmesh, cax = cbaxes)

   plt.show()

def pltsec(field,thk,lon1,lat1,lon2,lat2,*args, **kwargs) :
# field, thk are expected to be in 3D
   region = kwargs.get('region',None)
   region = 'TP5a0.06' if region is None else region
   print('setting domain to: ', region)
   plon,plat = getgrid(region=region)

   depthm = getdepth(region=region)

   cooINDEX = abs( plat-lat1 ) + abs( plon-lon1 )
   JJ1,II1 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
   cooINDEX = abs( plat-lat2 ) + abs( plon-lon2 )
   JJ2,II2 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)

#   sec = gridxsec.SectionIJSpace([lon1,lon2],[lat1,lat2],plon,plat)
   sec = gridxsec.SectionIJSpace([II1,II2],[JJ1,JJ2],plon,plat)
   I,J=sec.grid_indexes
   dist=sec.distance
   slon=sec.longitude
   slat=sec.latitude

   xaxis = kwargs.get('xaxis',None)
   if xaxis is not None:
     if xaxis == "lon" :
      x=slon
      xlab="Longitude"
     elif xaxis == "lat" :
      x=slat
      xlab="Latitude"
   else:
      x=dist/1000.
      xlab="Distance along section[km]"     

   fig=plt.figure(figsize=(13,4),facecolor='w')

   cmap = plt.cm.get_cmap('Spectral_r')


   kdm=field.shape[0]
   intfsec=np.zeros((kdm+1,I.size))
   datasec=np.zeros((kdm+1,I.size))

   for k in range(kdm) :
    dp2d = thk[k,:,:]
    data2d = field[k,:,:]
    dp2d=np.ma.filled(dp2d,0.)
    dp2d = dp2d/9860.
    data2d=np.ma.filled(data2d,1e30)
      
    intfsec[k+1,:] = intfsec[k,:] + dp2d[J,I]
    if k==0 : datasec[k,:] = data2d[J,I]
    datasec[k+1,:] = data2d[J,I]
   datasec = np.ma.masked_where(datasec>0.5*1e30,datasec)
   datasec = np.ma.masked_where(datasec<-0.5*1e10,datasec)

   cbaxes = fig.add_axes([0.6775, 0.12, 0.02, 0.55])
   ax = fig.add_subplot(2,1,1)
   ax.set_position([0.075,0.12,0.6,0.8])
   ax.set_facecolor('xkcd:gray')
   P=ax.pcolormesh(x,-intfsec,datasec,cmap=cmap)
   cmax = kwargs.get('cmax',None)
   cmax = datasec.max() if cmax is None else cmax
   cmin = kwargs.get('cmin',None)
   cmin = datasec.min() if cmin is None else cmin
   P.set_clim(cmin,cmax)
#  cbaxes = fig.add_axes([0.6775, 0.12, 0.02, 0.55])
   cb = ax.figure.colorbar(P, cax = cbaxes)

   dlim = kwargs.get('dlim',None)


   if dlim is not None:
      plt.ylim(-dlim,0)



   ax.set_ylabel("Depth (m)")
   ax.set_xlabel(xlab)

   title = kwargs.get('title',None)
   title = '' if title is None else title
   plt.title(title)

 

   axm=fig.add_subplot(212)
   axm.set_position([0.76,0.12,0.225,0.85])
   data2d = field[0,:,:]
   data2d = np.ma.masked_where(depthm.mask,data2d)
   axm.set_facecolor('xkcd:gray')
   if data2d.shape[1] == 216 :
    pmesh = plt.pcolormesh(np.rot90(data2d),cmap=cmap)
    section = plt.plot(J,215-I,'r',lw=2)
   else :
    pmesh = plt.pcolormesh(data2d,cmap=cmap)
    section = plt.plot(I,J,'r',lw=2)
   plt.clim(datasec.min(),datasec.max())

   plt.show()

def allplot(f) :
  dm=max(f.fieldlevels)

  no3 = getvarib(f,'no3',1)
  pho = getvarib(f,'pho',1)
  sil = getvarib(f,'sil',1)
  no3b = getvarib(f,'no3',dm)
  phob = getvarib(f,'pho',dm)
  silb = getvarib(f,'sil',dm)
  chl = getvarib(f,'chl',1)
  dia = getvarib(f,'ECO_dia',1)
  fla = getvarib(f,'ECO_fla',1)
  ccl = getvarib(f,'ECO_ccl',1)
  mic = getvarib(f,'ECO_micr',1)
  mes = getvarib(f,'ECO_meso',1)
  det = getvarib(f,'ECO_det',dm)
  opa = getvarib(f,'ECO_opa',dm)
  oxy = getvarib(f,'ECO_oxy',1)
  if str(f)[21:28] != 'Restart':
    pp  = getvarib(f,'pp',1)
  else:
    pass
  sed1= getvarib(f,'ECO_sed1',0)
  sed2= getvarib(f,'ECO_sed2',0)
  sed3= getvarib(f,'ECO_sed3',0)
  sed4= getvarib(f,'ECO_sed4',0)

  cmap = plt.cm.get_cmap('Spectral_r')

  fig=plt.figure(figsize=(18,10),facecolor='w')

  ax1 = fig.add_subplot(4,6,1)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax1.set_facecolor('xkcd:gray')
  pmesh1 = plt.pcolormesh(no3,cmap=cmap)
  ax1.set_title(r"Surface NO$_{3}$ ($\mu$molN L$^{-1}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh1)

  ax2 = fig.add_subplot(4,6,7)
#  levels = MaxNLocator(nbins=40).tick_values(0,5)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax2.set_facecolor('xkcd:gray')
  pmesh2 = plt.pcolormesh(pho,cmap=cmap)
  ax2.set_title(r"Surface PO$_{4}$ ($\mu$molP L$^{-1}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh2)

  ax3 = fig.add_subplot(4,6,13)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax3.set_facecolor('xkcd:gray')
  pmesh3 = plt.pcolormesh(sil,cmap=cmap)
  ax3.set_title(r"Surface Sil ($\mu$molSi L$^{-1}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh3)

  ax4 = fig.add_subplot(4,6,2)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax4.set_facecolor('xkcd:gray')
  pmesh4 = plt.pcolormesh(no3b,cmap=cmap)
  ax4.set_title(r"Bottom NO$_{3}$ ($\mu$molN L$^{-1}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh4)

  ax5 = fig.add_subplot(4,6,8)
#  levels = MaxNLocator(nbins=40).tick_values(0,5)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax5.set_facecolor('xkcd:gray')
  pmesh5 = plt.pcolormesh(phob,cmap=cmap)
  ax5.set_title(r"Bottom PO$_{4}$ ($\mu$molP L$^{-1}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh5)

  ax6 = fig.add_subplot(4,6,14)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax6.set_facecolor('xkcd:gray')
  pmesh6 = plt.pcolormesh(silb,cmap=cmap)
  ax6.set_title(r"Bottom Sil ($\mu$molSi L$^{-1}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh6)

  ax7 = fig.add_subplot(4,6,3)
#  levels = MaxNLocator(nbins=40).tick_values(0,5)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax7.set_facecolor('xkcd:gray')
  pmesh7 = plt.pcolormesh(chl,cmap=cmap)
  ax7.set_title(r"Surface Chl $a$ ($\mu$g L$^{-1}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh7)

  if str(f)[21:28] != 'Restart':
   ax8 = fig.add_subplot(4,6,9)
#   levels = MaxNLocator(nbins=40).tick_values(0,750)
#   norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
   ax8.set_facecolor('xkcd:gray')
   pmesh8 = plt.pcolormesh(pp*60*60*24,cmap=cmap)
   ax8.set_title(r"Integr. PP (mgC m$^{-2}$ d$^{-1}$)",fontsize=10)
   plt.xticks([])
   plt.yticks([])
   plt.colorbar(pmesh8)
  else:
   pass

  ax9 = fig.add_subplot(4,6,15)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax9.set_facecolor('xkcd:gray')
  pmesh9 = plt.pcolormesh(mic,cmap=cmap)
  ax9.set_title(r"Surface MicroZ (mgC m$^{-3}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh9)

  ax10 = fig.add_subplot(4,6,4)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax10.set_facecolor('xkcd:gray')
  pmesh10 = plt.pcolormesh(dia,cmap=cmap)
  ax10.set_title(r"Surface Diatoms (mgC m$^{-3}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh10)

  ax11 = fig.add_subplot(4,6,10)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax11.set_facecolor('xkcd:gray')
  pmesh11 = plt.pcolormesh(fla,cmap=cmap)
  ax11.set_title(r"Surface Flag. (mgC m$^{-3}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh11)

  ax12 = fig.add_subplot(4,6,16)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax12.set_facecolor('xkcd:gray')
  pmesh12 = plt.pcolormesh(mes,cmap=cmap)
  ax12.set_title(r"Surface MesoZ (mgC m$^{-3}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh12)

  ax13 = fig.add_subplot(4,6,5)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax13.set_facecolor('xkcd:gray')
  pmesh13 = plt.pcolormesh(oxy,cmap=cmap)
  ax13.set_title(r"Surface Oxy. (mgC m$^{-3} ??$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh13)

  ax14 = fig.add_subplot(4,6,11)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax14.set_facecolor('xkcd:gray')
  pmesh14 = plt.pcolormesh(det,cmap=cmap)
  ax14.set_title(r"Bottom Det. (mgC m$^{-3}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh14)

  ax15 = fig.add_subplot(4,6,17)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax15.set_facecolor('xkcd:gray')
  pmesh15 = plt.pcolormesh(opa,cmap=cmap)
  ax15.set_title(r"Bottom Opal (mgC m$^{-3}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh15)

  ax16 = fig.add_subplot(4,6,6)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax16.set_facecolor('xkcd:gray')
  pmesh16 = plt.pcolormesh(sed1,cmap=cmap)
  ax16.set_title(r"Sediment N (gC m$^{-2}$)",fontsize=10)#gC/m2
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh16)

  ax17 = fig.add_subplot(4,6,12)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax17.set_facecolor('xkcd:gray')
  pmesh17 = plt.pcolormesh(sed2,cmap=cmap)
  ax17.set_title(r"Sediment Si (gC m$^{-2}$)",fontsize=10)#gC/m2
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh17)

  ax18 = fig.add_subplot(4,6,18)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax18.set_facecolor('xkcd:gray')
  pmesh18 = plt.pcolormesh(sed3,cmap=cmap)
  ax18.set_title(r"Sediment P (gC m$^{-2}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh18)

  ax19 = fig.add_subplot(4,6,19)
#  levels = MaxNLocator(nbins=40).tick_values(0,500)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax19.set_facecolor('xkcd:gray')
  pmesh19 = plt.pcolormesh(ccl,cmap=cmap)
  ax19.set_title(r"Surface Coccoliths (mgC m$^{-3}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh19)

  ax20 = fig.add_subplot(4,6,20)
#  levels = MaxNLocator(nbins=40).tick_values(0,50)
#  norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
  ax20.set_facecolor('xkcd:gray')
  pmesh20 = plt.pcolormesh(sed4,cmap=cmap)
  ax20.set_title(r"Sediment CaCO3 (gC m$^{-2}$)",fontsize=10)
  plt.xticks([])
  plt.yticks([])
  plt.colorbar(pmesh20)

'''
SOME EXAMPLE GENERAL PURPOSE SCRIPTS 

#-----------------------------------------------------------------
#plot TP2 Ob river currents and salinity - multiple days for animation
#-----------------------------------------------------------------
import os,fnmatch
import abplot
import datetime
import cmocean
from netCDF4 import Dataset as NetCDFFile

directory = '/cluster/work/users/cagyum/TP2a0.10/expt_01.0/data/'
#ndays = numpy.size(sorted(fnmatch.filter(os.listdir(directory), 'archm.%s*.b' %('20[01]*') ) ))

for f1 in sorted(fnmatch.filter(os.listdir(directory), 'archm.%s*.b' %('20[01]*') ) ) :
 print(f1[6:14])
 year = f1[6:10]
 day  = f1[11:14]
 d = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(day) - 1)
 icefile = directory+'cice/iceh.'+str(d.year)+'-'+str(d.month).zfill(2)+'-'+str(d.day).zfill(2)+'.nc'
 nc = NetCDFFile(icefile)
 icea = nc.variables['aice_d'][0,:,:]
 nc.close()


 ab = abplot.openfile(directory+f1)
 sal = abplot.getvarib(ab,"salin",1)
 no3 = abplot.getvarib(ab,"no3",1)
 pho = abplot.getvarib(ab,"pho",1)
 sil = abplot.getvarib(ab,"sil",1)
 chl = abplot.getvarib(ab,"chl",1)
 u = abplot.getvarib(ab,"u-vel.",1)
 v = abplot.getvarib(ab,"v-vel.",1) 
 oxy = abplot.getvarib(ab,"ECO_oxy",1)

 uplot = u[225:275,300:372]
 vplot = v[225:275,300:372]
 fig=plt.figure(figsize=(12,6),facecolor='w')
 ax = fig.add_subplot(2,4,1)
 ax.set_facecolor('xkcd:gray')
 cmap = plt.get_cmap('Spectral_r')
 pmesh1 = plt.pcolormesh((vplot**2+uplot**2)**0.5,cmap=cmap)
 pmesh1.set_clim(0.,0.25)
 cb=plt.colorbar(pmesh1)
 plt.title('Speed '+f1[6:14])

 xind = np.arange(uplot.shape[0])+0.5
 yind = np.arange(uplot.shape[1])+0.5
 xind2d = np.zeros(uplot.shape)
 for y in range(yind.shape[0]):
     xind2d[:,y] = xind
 yind2d = np.zeros(uplot.shape)
 for x in range(xind.shape[0]):
     yind2d[x,:] = yind
 plt.quiver(yind2d,xind2d,uplot,vplot)

 ax2 = fig.add_subplot(2,4,2)
 ax2.set_facecolor('xkcd:gray')
 pmesh2 = plt.pcolormesh(sal[225:275,300:372],cmap=cmap)
 pmesh2.set_clim(0.,10.)
 cb2=plt.colorbar(pmesh2)
 plt.title('Salinity '+f1[6:14])
 
 ax3 = fig.add_subplot(2,4,3)
 ax3.set_facecolor('xkcd:gray')
 pmesh3 = plt.pcolormesh(chl[225:275,300:372],cmap=cmap)
 pmesh3.set_clim(0.,20.)
 cb3=plt.colorbar(pmesh3)
 plt.title('Chlorophyll '+f1[6:14])

 ax4 = fig.add_subplot(2,4,4)
 ax4.set_facecolor('xkcd:gray')
 pmesh4 = plt.pcolormesh(no3[225:275,300:372],cmap=cmap)
 pmesh4.set_clim(0.,30.)
 cb4=plt.colorbar(pmesh4)
 plt.title('Nitrate '+f1[6:14])

 ax5 = fig.add_subplot(2,4,7)
 ax5.set_facecolor('xkcd:gray')
 pmesh5 = plt.pcolormesh(sil[225:275,300:372],cmap=cmap)
 pmesh5.set_clim(0.,100.)
 cb5=plt.colorbar(pmesh5)
 plt.title('Silicate '+f1[6:14])

 ax6 = fig.add_subplot(2,4,8)
 ax6.set_facecolor('xkcd:gray')
 pmesh6 = plt.pcolormesh(pho[225:275,300:372],cmap=cmap)
 pmesh6.set_clim(0.,3.)
 cb6=plt.colorbar(pmesh6)
 plt.title('Phosphate '+f1[6:14])

 ax7 = fig.add_subplot(2,4,5)
 ax7.set_facecolor('xkcd:gray')
 pmesh7 = plt.pcolormesh(oxy[225:275,300:372],cmap=cmap)
 pmesh7.set_clim(300.,475.)
 cb7=plt.colorbar(pmesh7)
 plt.title('Oxygen '+f1[6:14])

 ax8 = fig.add_subplot(2,4,6)
 ax8.set_facecolor('xkcd:gray')
 pmesh8 = plt.pcolormesh(icea[225:275,300:372],cmap=cmap)
 pmesh8.set_clim(0.,1.)
 cb8=plt.colorbar(pmesh8)
 plt.title('Ice area '+f1[6:14])

 plt.tight_layout()
 ab.close()
 
 fig.canvas.print_figure(directory+'vector_plot/vector'+f1[6:14]+'.png',dpi=180)
 plt.close(fig)
'''

'''
#-----------------------------------------------------------------
#plot TP2 Ob river flux at the outer bay
#-----------------------------------------------------------------
import os,fnmatch
import abplot
import numpy as np
import abfile
import datetime

directory = '/cluster/work/users/cagyum/TP2a0.10/expt_01.0/data/'
total_all = []
date = []
abgrid = abfile.ABFileGrid(directory + "../../topo/regional.grid","r")
y=abgrid.read_field("scpy")

for f1 in sorted(fnmatch.filter(os.listdir(directory), 'archm.%s*.b' %('20[01]*') ) ) :

 year = f1[6:10]
 day  = f1[11:14]
 d = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(day) - 1)

 ab = abplot.openfile(directory+f1)
 u  = abplot.get3d(ab,"u-vel.")
 th = abplot.get3d(ab,"thknss")/9806.
 ab.close()

 total = np.sum( u[:,250:255,334] * th[:,250:255,334] * y[250:255,334] )/1e6
 total_all.append(total)
 date.append(d)

 print(f1[6:14],total,sum(total_all))

f = open('Ob_Flux.pckl','wb')
pickle.dump([date,total_all],f)
close(f)
'''
'''
#-----------------------------------------------------------------
#plot TP2 Ob river bay nutrient mass
#-----------------------------------------------------------------
import os,fnmatch
import abplot
import numpy as np
import abfile
import datetime

directory = '/cluster/work/users/cagyum/TP2a0.10/expt_01.0/SCRATCH/'
total_no3 = []; total_sil = []
date = []
abgrid = abfile.ABFileGrid(directory + "../../topo/regional.grid","r")
y=abgrid.read_field("scpy")
x=abgrid.read_field("scpx")
plat = abgrid.read_field("plat")
plon = abgrid.read_field("plon")
jdm,idm=plon.shape
abdepth = abfile.ABFileBathy(directory+"../../topo/regional.depth.b","r",idm=idm,jdm=jdm)
depthm=abdepth.read_field("depth")

masked = np.copy(depthm)
masked = np.ma.masked_where(depthm.mask,masked)
masked = np.ma.masked_where(plat>=74.,masked)
masked = np.ma.masked_where(plon<=70.,masked)
masked = np.ma.masked_where(plon>=90.,masked)

for f1 in sorted(fnmatch.filter(os.listdir(directory), 'archm.%s*.b' %('20[01]*') ) ) :

 year = f1[6:10]
 day  = f1[11:14]
 d = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(day) - 1)

 ab = abplot.openfile(directory+f1)
 sil  = abplot.get3d(ab,"sil")
 no3 = abplot.get3d(ab,"no3")
 sed1 = abplot.getvarib(ab,"ECO_sed1",0)
 sed2 = abplot.getvarib(ab,"ECO_sed2",0)
 th = abplot.get3d(ab,"thknss")/9806.
 ab.close()

 no3_sum = np.sum(no3 * th * x * y,axis=0) + sed1 * x * y
 no3_sum = np.ma.masked_where(masked.mask,no3_sum)
 total_no3.append(np.sum(no3_sum))

 sil_sum = np.sum(sil * th * x * y,axis=0) + sed2 * x * y
 sil_sum = np.ma.masked_where(masked.mask,sil_sum)
 total_sil.append(np.sum(sil_sum))

 date.append(d)

 print(f1[6:14],np.sum(sil_sum)/10e9,np.sum(no3_sum)/10e9)

f = open('Ob_Mass.pckl','wb')
pickle.dump([date,total_no3,total_sil],f)
close(f)
'''
