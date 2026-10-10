import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from netCDF4 import Dataset as NetCDFFile
import datetime
import warnings
import math
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from dateutil.relativedelta import relativedelta
from scipy.interpolate import griddata
import pickle
import os.path
#import seawater as sw
#from seawater.library import T90conv
warnings.filterwarnings('ignore')

# SOME STATISTICS #
def pc_bias(s,o):
        return 100.0*sum(s-o)/sum(o)
def bias(s,o):
        return np.mean(s-o)
def correlation(s,o):
        return np.corrcoef(o, s)[0,1]
def rmse(s,o):
        return np.sqrt(np.mean((s-o)**2))
def nstd(s,o):
        return np.std(s)/np.std(o)
def covar(s,o):
        return np.cov(s,o,bias=True)[0,1]
def varian(s):
        return np.var(s)

def distance_on_unit_sphere(lat1, long1, lat2, long2):

# Convert latitude and longitude to
# spherical coordinates in radians. Assumes the Earth is a perfect sphere.
  degrees_to_radians = math.pi/180.0
# phi = 90 - latitude
  phi1 = (90.0 - lat1)*degrees_to_radians
  phi2 = (90.0 - lat2)*degrees_to_radians
# theta = longitude
  theta1 = long1*degrees_to_radians
  theta2 = long2*degrees_to_radians

  cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
  math.cos(phi1)*math.cos(phi2))
  arc = math.acos( cos )
  arc = arc * 6378.137 # convert to km

  return arc


def readARGO(number, variable,*args, **kwargs):
    # file and variable names are compatible with Copernicus NetCDF files
    directory = "./"
    filename = directory + "GL_PR_PF_" + str(number) + ".nc"
    nc = NetCDFFile(filename)

    pres    = nc.variables['PRES'][:,:]
    presqc  = nc.variables['PRES_QC'][:,:]
    lon     = nc.variables['LONGITUDE'][:]
    lat     = nc.variables['LATITUDE'][:]
 
    if 'JULD' in nc.variables:
       IFREMER = True
       if variable == 'CPHL_ADJUSTED' : variable = 'CHLA_ADJUSTED'
       if variable == 'NTAW_ADJUSTED' : variable = 'NITRATE_ADJUSTED'
       if variable == 'DOX2_ADJUSTED' : variable = 'DOXY_ADJUSTED'       
    else:
       IFREMER = False

    var = []; varqc = []
    if variable in nc.variables:

       var   = nc.variables[variable][:,:]
       varqc = nc.variables[variable+'_QC'][:,:]
#       error = kwargs.get('error',None)
#       if error is not None:
#          varer = nc.variables[variable+'_ERROR'][:,:]

       if variable != 'TEMP_ADJUSTED':
          var[var<=0.] = 0.

       if IFREMER:
          var = np.ma.masked_where(presqc != b'1', var)
          pres = np.ma.masked_where(presqc != b'1', pres)          
       else:
          var = np.ma.masked_where(presqc != 1, var)
          pres = np.ma.masked_where(presqc != 1, pres)

#    else:
#       print(variable + ' is not listed in the NetCDF File')
    elif variable == 'poc':
       for variable in ["CPHL_ADJUSTED"]:
         try:
            chl, chlqc = readARGO(number, variable)
         except Exception:
            pass
       chl  = construct(chl,chlqc,[1,5,8])
       for variable in ["BBP700_ADJUSTED"]:
         try:
            bbp, bbpqc = readARGO(number, variable)
         except Exception:
            pass
       bbp  = construct(bbp,bbpqc,[1,5,8])

       chlm = np.ma.masked_where( bbp.mask, chl)
       bbpm = np.ma.masked_where( chlm.mask, bbp) 

       k1 = 89.423; k2 = 0.1881; k3 = 0.7591; k4 = 0.1934; e1 = 1.626; e2 = 21.2; pocm = 33.4
       poc = k1 * bbpm**k2 * (chlm/bbpm)**k3 * (chlm/bbpm)**(k4*np.log10(bbpm))
       poc[poc<pocm] = e1 * poc[poc<pocm] - e2
       poc = np.ma.masked_where(poc>10000.,poc)
       var = np.ma.masked_where(poc<0.,poc)
       varqc = np.zeros((var.shape))
       varqc[~var.mask] = 1

    elif variable == 'pocK':
       for variable in ["CPHL_ADJUSTED"]:
         try:
            chl, chlqc = readARGO(number, variable)
         except Exception:
            pass
       chl  = construct(chl,chlqc,[1,5,8])
       for variable in ["BBP700_ADJUSTED"]:
         try:
            bbp, bbpqc = readARGO(number, variable)
         except Exception:
            pass
       bbp  = construct(bbp,bbpqc,[1,5,8])

       chlm = np.ma.masked_where( bbp.mask, chl)
       bbpm = np.ma.masked_where( chlm.mask, bbp)


       var, dummy1, dummy2 = Koest23_modelB_700p((bbpm),chlm)
       varqc = np.zeros((var.shape))
       varqc[~var.mask] = 1

    else:
       print(variable + ' is not listed in the NetCDF File')

    if IFREMER:
       time    = nc.variables['JULD'][:]
    else:
       time    = nc.variables['TIME'][:]
    time2d = np.zeros((pres.shape))
    for z in range(pres.shape[1]):
        time2d[:,z] = time
    if IFREMER:
       since = int(nc.variables['JULD'].units[11:15])
    else:
       since = int(nc.variables['TIME'].units[11:15])
    dates = list()
    for t in range(time.shape[0]):
        dates.append( datetime.datetime(since,1,1) + datetime.timedelta(days=time[t]) )
    datestart = dates[0]
    dateend   = dates[-1]

    diff = datestart-datetime.datetime(since,1,1,0,0,0)
    hour1 = int(diff.days * 24. + diff.seconds/60./60.)
    diff = dateend-datetime.datetime(since,1,1,0,0,0)
    hour2 = int(diff.days * 24. + diff.seconds/60./60.)

    timerange = np.concatenate((float(hour1),time*24.),axis=None)
    timerange = np.concatenate((timerange,float(hour2)),axis=None)
    latrange = np.concatenate((lat[0],lat),axis=None)
    latrange = np.concatenate((latrange,lat[-1]),axis=None)
    lonrange = np.concatenate((lon[0],lon),axis=None)
    lonrange = np.concatenate((lonrange,lon[-1]),axis=None)
    latint = np.interp(np.arange(hour1,hour2+1,1),timerange,latrange)
    lonint = np.interp(np.arange(hour1,hour2+1,1),timerange,lonrange)
    timeint = np.arange(hour1,hour2+1,1)


    if kwargs.get('all'):
       return pres, since, time2d, dates, lon, lonint, lat, latint, timeint, var, varqc
    elif kwargs.get('simple'):
       return pres, time2d, dates, var, varqc
    else:
#       if error is not None:
#          return var, varqc, varer
#       else:
       return var, varqc
      


def select_flags(var,qc,flag, *args, **kwargs):
    IFREMER=False
    try :
        type(qc.max())
    except  Exception:
        IFREMER=True
    if IFREMER:
       flag = str(flag).encode()
    new = np.zeros(var.shape) - 1E10
    new[qc==flag] = var[qc==flag]
    new = np.ma.masked_where(new < -1E5, new)
    return new


def construct(var,qc,flags, *args, **kwargs):
    IFREMER=False
    try :
        type(qc.max())
    except  Exception:
        IFREMER=True
    new = np.zeros(var.shape) - 1E10
    n = np.size(flags)
    for i in range(n):
        if IFREMER:
           flag = str(flags[i]).encode()
           new[qc==flag] = var[qc==flag]
        else:
           new[qc==flags[i]] = var[qc==flags[i]]
    new = np.ma.masked_where(new < -1E5, new)
    return new


def interpol(time2d,depth,since,var):
    time = time2d[:,0]
    x = list({int(d):d for d in time}.values()) 

    d1d = depth.flatten()
    v1d = var.flatten()
    t1d = time2d.flatten()

    t1d = t1d[~d1d.mask]
    v1d = v1d[~d1d.mask]
    d1d = d1d[~d1d.mask]

    t1d = t1d[~v1d.mask]
    d1d = d1d[~v1d.mask]
    v1d = v1d[~v1d.mask]

    y = np.arange(0.,float(int(d1d.max())),0.5)     
    xm,ym = np.meshgrid(x,y)

    intp = griddata( (d1d,t1d), v1d, (ym,xm), method='nearest')

    dates = list()
    for t in range(len(x)):
        dates.append( datetime.datetime(since,1,1) + datetime.timedelta(days=x[t]) )

    return dates,-y,intp

def find_next(argodates,enddate,index):
    index_first = index
    while argodates[index] <= enddate:
          if index >= len(argodates):
             index = index_first
             exit()
          index = index + 1
    print('Forcing end date is set to: ',argodates[index])
    print('This is done to enabling the next profile (if available) for interpolation')
    return index

def find_previous(argodates,startdate,index):
    index_first = index
    while argodates[index] >= startdate:
          if index == 0:
             index = index_first
             exit()
          index = index - 1
    print('Forcing start date is set to: ',argodates[index])
    print('This is done to enabling the next profile (if available) for interpolation')
    return index

def get_title(variable):

   if variable == "chl":
      title = r"Chlorophyll $a$ (mgChl m$^{-3}$)"
   elif variable == "nit":
      title = r"Nitrate (mmolN m$^{-3}$)"
   elif variable == "oxy":
      title = r"Oxygen (mmol m$^{-3}$)" 
   elif variable == "poc":
      title = r"POCs (mgC m$^{-3}$)"
   else:
      title = variable
   
   return title

def get_colormap(variable):
    import cmocean
    import nclcmaps
    if variable == "chl" or variable == "total_chlorophyll_calculator_result":
       #cmap = cmocean.cm.algae_r
       cmap = plt.cm.get_cmap('Spectral_r')
       #cmap = nclcmaps.cmap('WhiteBlueGreenYellowRed')
       #cmap = cmocean.cm.curl
       #cmap = cmocean.cm.speed_r
    elif variable == "nit":
       #cmap = cmocean.cm.balance
       cmap = plt.cm.get_cmap('Spectral_r')
    elif variable == "temp":
       cmap = cmocean.cm.thermal
    elif variable == "sal":
       cmap = cmocean.cm.haline
    elif variable == "oxy":
       #cmap = cmocean.cm.curl 
       cmap = plt.cm.get_cmap('Spectral_r')
    elif variable == "poc":
       #cmap = cmocean.cm.curl 
       cmap = plt.cm.get_cmap('Spectral_r')
    else: 
       cmap = plt.cm.get_cmap('Spectral_r')
   
    return cmap

def pmeshargo2(fld,depth,dates, *args, **kwargs):
    plt.style.use('seaborn-v0_8-notebook')

    first = kwargs.get('first',None)
    if kwargs.get('first'):
       y = int(first[0:4]); m = int(first[5:7]); d = int(first[8:10])
       plotmin = datetime.datetime(y,m,d)
    else: 
       plotmin = dates[0]

    last = kwargs.get('last',None)
    if kwargs.get('last'):
       y = int(last[0:4]); m = int(last[5:7]); d = int(last[8:10])
       plotmax = datetime.datetime(y,m,d)
    else: 
       plotmax = dates[-1]

    dtop = kwargs.get('dtop',None)
    dtop = max(depth) if dtop is None else -dtop
    dbot = kwargs.get('dbot',None)
    dbot = min(depth) if dbot is None else -dbot 

    variable = kwargs.get('variable',None) 
    cmap = get_colormap(variable)

    ind1 = np.abs(dtop - depth[:]).argmin() ; ind1 = min(ind1, depth.shape[0]-1)
    ind2 = np.abs(dbot - depth[:]).argmin() ; ind2 = max(ind2, depth.shape[0])
    depth = depth[ind1:ind2+1]
    fld = fld[ind1:ind2+1,:]

    cmax = kwargs.get('cmax',None)
    cmax = fld.max() if cmax is None else cmax
    cmin = kwargs.get('cmin',None)
    cmin = fld.min() if cmin is None else cmin
     
    fig = plt.figure(figsize=(10,3.5),facecolor='w')
    ax  = fig.add_subplot(1,1,1)
    ax.set_position([0.125,0.175,0.75,0.725])
    pmesh = plt.pcolormesh(dates,depth,fld,cmap=cmap,shading='auto')

    ax.xaxis.axis_date()
    ax.yaxis.set_tick_params(labelsize='large')
    ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)

    plt.clim(cmin,cmax)
    plt.clim(cmin,cmax)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dbot,dtop])
    ax.grid(color='k', linestyle=':', linewidth=1)

    if variable is not None:
      plt.title(get_title(variable), loc='left', fontsize=15)

    ax.set_ylabel("Depth (m)",fontsize=15)

    cbaxes = fig.add_axes([0.9, 0.25, 0.025, 0.5])
    cb = plt.colorbar(pmesh, cax = cbaxes)
    cbaxes.yaxis.set_tick_params(labelsize='large')  

def pmeshargo(fld,depth,dates, *args, **kwargs):
    plt.style.use('seaborn-v0_8')
    import cmocean

    if kwargs.get('first'):
       plotmin = datetime.datetime(first,1,1)
    else: 
       plotmin = dates[0]

    if kwargs.get('last'):
       plotmax = datetime.datetime(last+1,1,1)
    else: 
       plotmax = dates[-1]

    nyears = (plotmax - plotmin).days / 365.
    if nyears ==1 :
        interval=1
    elif nyears == 2 or nyears == 3:
        interval = 3
    else :
        interval = 6

    cmax = kwargs.get('cmax',None)
    cmax = fld.max() if cmax is None else cmax
    cmin = kwargs.get('cmin',None)
    cmin = fld.min() if cmin is None else cmin
    dtop = kwargs.get('dtop',None)
    dtop = (depth).max() if dtop is None else -dtop
    dbot = kwargs.get('dbot',None)
    dbot = (depth).min() if dbot is None else -dbot 

    cmapv = kwargs.get('cmapv',None)
    if cmapv != None:
       if cmapv == 'oxygen':
          cmap = cmocean.cm.oxy
       else:
          cmap = plt.cm.get_cmap('Spectral_r')
    else:
       cmap = plt.cm.get_cmap('Spectral_r')


    fig = plt.figure(figsize=(10,3.5),facecolor='w')
    ax  = fig.add_subplot(1,1,1)
    ax.set_position([0.125,0.125,0.75,0.8])
    pmesh = plt.pcolormesh(dates,depth,fld,cmap=cmap,shading='auto')
    plt.clim(cmin,cmax)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dbot,dtop])

    fmt_year = mdates.MonthLocator(interval=interval)
    ax.xaxis.set_major_locator(fmt_year)
    # Minor ticks every month.
    fmt_month = mdates.MonthLocator()
    ax.xaxis.set_minor_locator(fmt_month)
    # Text in the x axis will be displayed in 'YYYY-mm' format.
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

    ax.set_ylabel("Depth (m)",fontsize=13)
    #ax.grid(which='major', axis='both',color='w')
    cbaxes = fig.add_axes([0.9, 0.25, 0.025, 0.5])
    cb = plt.colorbar(pmesh, cax = cbaxes)


def trajectory(number, *args, **kwargs):
    import cartopy.crs as ccrs
    import cartopy.feature

    pres, since, time2d, dates, lon, lonint, lat, latint, timeint,temp, tempqc = readARGO(number, 'TEMP_ADJUSTED',all=True)

    time = time2d[:,0] - (datetime.datetime(1970,1,1)-datetime.datetime(since,1,1)).days

    warnings.filterwarnings('ignore')
    fig = plt.figure(figsize=[8.75, 5])
    ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude=lon.mean(), central_latitude=lat.mean()))
    #ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=lon.mean()))
#    ax.set_extent([-60, 30, 50, 80], ccrs.PlateCarree())
#    ax.set_extent([np.min([-60,lon.min()-1.]), np.min([lat.min()-1.,30]), \
#                  np.max([50,lon.max()+1.]), np.max([80.,lat.max()+1.])], ccrs.PlateCarree())
    ax.set_extent([  np.max( [ lon.min()-5., -180. ] ),\
                     np.min( [ lon.max()+5., 180.  ] ),\
                     np.max( [ lat.min()-5., -90.  ] ),\
                     np.min( [ lat.max()+5., 90.   ] )],\
                     ccrs.PlateCarree())


    ax.set_position([0.125,0.125,0.7,0.8])
    fig.subplots_adjust(bottom=0.075, top=0.925,
                        left=0.08, right=0.92, wspace=0.02)

    gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)


    url = kwargs.get('url',None) 
    if url is not None:
       import cmocean
       cmap = cmocean.cm.gray
       nc=NetCDFFile(url)
       nclat = nc.variables['lat'][:]
       nclon = nc.variables['lon'][:]

       IImin = np.abs(lon.min()-nclon).argmin() - 250
       IImax = np.abs(lon.max()-nclon).argmin() + 250
       JJmin = np.abs(lat.min()-nclat).argmin() - 250
       JJmax = np.abs(lat.max()-nclat).argmin() + 250       

       elevation = nc.variables['elevation'][JJmin:JJmax+250,IImin:IImax+250]
       elevation[elevation>0.]=0.
       print(elevation.min(),elevation.max())
       el = plt.pcolormesh(nclon[IImin:IImax+250],nclat[JJmin:JJmax+250],elevation,transform=ccrs.PlateCarree(),cmap=cmap)

    ax.add_feature(cartopy.feature.LAND)   # land is specified to plot above ...
    ax.coastlines(resolution='50m')

    cmap = plt.get_cmap('Spectral_r')
    p = plt.scatter(lon,lat,s=10,c=time,cmap=cmap,marker='.',transform=ccrs.PlateCarree())
    cbaxes = fig.add_axes([0.88, 0.2, 0.025, 0.6])
    cb = plt.colorbar(p, cax = cbaxes)
    if url is not None:
       cbaxel = fig.add_axes([0.05, 0.2, 0.025, 0.6])
       cbel = plt.colorbar(el, cax = cbaxel)


    fmt_qu_year = mdates.MonthLocator(interval=3)
    cbaxes.yaxis.set_major_locator(fmt_qu_year)
    fmt_month = mdates.MonthLocator()
    cbaxes.yaxis.set_minor_locator(fmt_month)
    cbaxes.yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

    #plt.tight_layout()
    if kwargs.get('print'):
       fig.canvas.print_figure("trajectory.png",dpi=210)
       plt.close(fig)

def get_woa(lat,lon,month,varib,*args, **kwargs):
    if varib == 'temperature': initial = 't'
    if varib == 'salinity': initial = 's'
    if varib == 'nitrate': initial = 'n'
    if varib == 'silicate': initial = 'i'
    if varib == 'phosphate': initial = 'p'
    if varib == 'oxygen': initial = 'o'

    tries = 1000
    if initial == 't' or initial == 's':
       for i in range(tries):
           try:
              ncy = NetCDFFile('https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/'+ varib +\
                '/A5B7/1.00/woa18_A5B7_'+ initial +'00_01.nc')
              nclat = ncy.variables['lat'][:]
              nclon = ncy.variables['lon'][:]
              II = np.abs(lon-nclon).argmin()
              JJ = np.abs(lat-nclat).argmin()
              vary = ncy.variables[initial + '_an'][0,:,JJ,II]
              break
           except Exception:
              if i < tries - 1:
                 print("ERROR, try number : ",i+1)
                 continue
              else:
                 raise RuntimeError("Tried a 1000 times, wait and retry")
              break

       for i in range(tries):
           try:          
              ncm = NetCDFFile('https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/'+ varib +\
                '/A5B7/1.00/woa18_A5B7_'+ initial + str(month).zfill(2) + '_01.nc')
              varm = ncm.variables[initial + '_an'][0,:,JJ,II]
              break
           except Exception:
              if i < tries - 1:
                 print("ERROR, try number : ",i+1)
                 continue
              else:
                 raise RuntimeError("Tried a 1000 times, wait and retry")
              break

    else:
       for i in range(tries):
           try:
              ncy = NetCDFFile('https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/'+ varib +\
                '/all/1.00/woa18_all_'+ initial +'00_01.nc')
              nclat = ncy.variables['lat'][:]
              nclon = ncy.variables['lon'][:]
              II = np.abs(lon-nclon).argmin()
              JJ = np.abs(lat-nclat).argmin()
              vary = ncy.variables[initial + '_an'][0,:,JJ,II]
              break
           except Exception:
              if i < tries - 1:
                 print("ERROR, try number : ",i+1)
                 continue
              else:
                 raise RuntimeError("Tried a 1000 times, wait and retry")
              break

       for i in range(tries):
           try:
              ncm = NetCDFFile('https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/'+ varib +\
                '/all/1.00/woa18_all_'+ initial + str(month).zfill(2) + '_01.nc')
              varm = ncm.variables[initial + '_an'][0,:,JJ,II]
              break
           except Exception:
              if i < tries - 1:
                 print("ERROR, try number : ",i+1)
                 continue
              else:
                 raise RuntimeError("Tried a 1000 times, wait and retry")
              break


    depy = ncy.variables['depth'][:]
    depm = ncm.variables['depth'][:]

    finald = np.copy(depm)
    finalv = np.copy(varm)

    finald = np.concatenate((finald,depy[np.size(depm):]),axis=None)
    finalv = np.concatenate((finalv,vary[np.size(depm):]),axis=None)

    if np.ma.is_masked(vary):
       finald = np.copy(finald[~vary.mask])
       finalv = np.copy(finalv[~vary.mask])

    finald = np.concatenate((finald,-12000.),axis=None)
    finalv = np.concatenate((finalv,finalv[-1]),axis=None)

    #make sure that depth is negative
    finald = -np.abs(finald)
    return finalv, finald

def write_woa(lat,lon,since,timeargo,start,datemid,end,varib):

    # start and end expected format : '2000-01-01'
    starty = int(start[0:4]); startm = int(start[5:7]); startd = int(start[8:10])
    endy   = int(end[0:4]); endm = int(end[5:7]); endd = int(end[8:10])
   # midy   = int(mid[0:4]); midm = int(mid[5:7]); midd = int(mid[9:11])

    datestart = datetime.datetime(starty,startm,startd) - relativedelta(months=1) 
   # datemid = datetime.datetime(midy,midm,midd)
    dateend = datetime.datetime(endy,endm,endd) + relativedelta(months=1)

    data = {}
    for t in range(12):
        print('reading '+varib+' spinup point lat: '+str(round(lat[0],2))+' lon: '+str(round(lon[0],2))+' @month: ',t+1)
        v,d = get_woa(lat[0],lon[0],t+1,varib)
        data[t] = {'v': v, 'd': d }

    if varib == 'temperature':
       filename = 'tprof'
    elif varib == 'salinity':
       filename = 'sprof'
    else:
       filename = varib

    with open(varib+'.dat','w') as f:    

         while datestart <= datemid :
               datemonth = datetime.datetime(datestart.year,datestart.month,15)
               printmonth = datestart.month - 1
               f.write('{} {} {}'.format(str(datemonth),str(data[printmonth]['v'].shape[0]),'2\n'))
               for l in range(data[printmonth]['v'].shape[0]):
                   f.write('{} {}'.format(str(data[printmonth]['d'][l]),str(data[printmonth]['v'][l])+'\n'))
               datestart = datestart + relativedelta(months=1)
         copied=datemonth # store previous loops datemonth 
         while datemid <= dateend + relativedelta(months=1) :               
               datemonth = datetime.datetime(datemid.year,datemid.month,15)
               if datemonth <= copied: 
                  datemonth = datemonth + relativedelta(months=1)
                  datemid = datemid + relativedelta(months=1)
               ttt = np.abs(( datemonth-datetime.datetime(since,1,1) ).days - timeargo ).argmin()
               latmid = np.interp(( datemonth-datetime.datetime(since,1,1) ).days,timeargo,lat)
               lonmid = np.interp(( datemonth-datetime.datetime(since,1,1) ).days,timeargo,lon)
               print('reading '+varib+' @lat: '+str(latmid)+' lon: '+str(lonmid)+' @month: '+str(datemonth.month))
               v,d = get_woa(latmid,lonmid,datemonth.month,varib)
               f.write('{} {} {}'.format(str(datemonth),str(v.shape[0]),'2\n'))
               for l in range(v.shape[0]):
                   f.write('{} {}'.format(str(d[l]),str(v[l])+'\n'))
               datemid = datemid + relativedelta(months=1)

def prepare_woa_input(lat,lon,since,timeargo,start,mid,end):

    for variable in ['nitrate','silicate','phosphate','oxygen','temperature','salinity']:
        if os.path.exists("./"+str(variable)+".dat"):
           print(str(variable)+".dat exists. Make sure that it is correct")
        else:
           write_woa(lat,lon,since,timeargo,start,mid,end,variable)


def get_era5(JJ,II,timecal, *args, **kwargs):
    print(str(timecal),end='\r')
    url = './ERA5.GOTM.'+str(int(timecal.year))+'.nc'
    nc=NetCDFFile(url)

    nctime = nc.variables['time'][:]
    since = int(nc.variables['time'].units[12:16])
    diff = timecal-datetime.datetime(since,1,1,0,0,0)
    hour = int(diff.days * 24. + diff.seconds/60./60.)
    tim1 = np.abs(hour - nctime).argmin() 

#    tp = nc.variables['mtpr'][tim1,JJ,II]; tp = tp / 1000   
#    msl = nc.variables['msl'][tim1,JJ,II]; msl = msl * 0.01
#    t2m = nc.variables['t2m'][tim1,JJ,II]; t2m = t2m - 273.15
#    d2m = nc.variables['d2m'][tim1,JJ,II]; d2m = d2m - 273.15
#    v10 = nc.variables['v10'][tim1,JJ,II]; v10 = max(v10,0.001)
#    u10 = nc.variables['u10'][tim1,JJ,II]; u10 = max(u10,0.001)
#    tcc = nc.variables['tcc'][tim1,JJ,II]; tcc = min(tcc,0.999) 
#    swr = nc.variables['msnswrf'][tim1,JJ,II]; swr = max(swr,0.0)

    tp  = interpolate_forcing(hour,nctime,tim1,nc,'mtpr',II,JJ); tp = tp / 1000.
    msl = interpolate_forcing(hour,nctime,tim1,nc,'msl',II,JJ); msl = msl * 0.01
    t2m = interpolate_forcing(hour,nctime,tim1,nc,'t2m',II,JJ); t2m = t2m - 273.15
    d2m = interpolate_forcing(hour,nctime,tim1,nc,'d2m',II,JJ); d2m = d2m - 273.15
    v10 = interpolate_forcing(hour,nctime,tim1,nc,'v10',II,JJ); v10 = max(v10,0.001)
    u10 = interpolate_forcing(hour,nctime,tim1,nc,'u10',II,JJ); u10 = max(u10,0.001)
    tcc = interpolate_forcing(hour,nctime,tim1,nc,'tcc',II,JJ); tcc = min(tcc,0.999)
    swr = interpolate_forcing(hour,nctime,tim1,nc,'msnswrf',II,JJ); swr = max(swr,0.0)

    nc.close()

    return tp,msl,t2m,d2m,v10,u10,tcc,swr

def get_ecnc(JJ,II,timecal, *args, **kwargs):
    print(str(timecal),end='\r')
    url = './ECNC.GOTM.'+str(int(timecal.year))+'.nc'
    nc=NetCDFFile(url)

    nctime = nc.variables['time'][:]
    since = int(nc.variables['time'].units[12:16])
    diff = timecal-datetime.datetime(since,1,1,0,0,0)
    hour = int(diff.days * 24. + diff.seconds/60./60.)
    tim1 = np.abs(hour - nctime).argmin() 

    tp  = interpolate_forcing(hour,nctime,tim1,nc,'TP',II,JJ); tp = tp /3600. / 6.
    msl = interpolate_forcing(hour,nctime,tim1,nc,'MSL',II,JJ); msl = msl * 0.01
    t2m = interpolate_forcing(hour,nctime,tim1,nc,'T2M',II,JJ); t2m = t2m - 273.15
    d2m = interpolate_forcing(hour,nctime,tim1,nc,'D2M',II,JJ); d2m = d2m - 273.15
    v10 = interpolate_forcing(hour,nctime,tim1,nc,'V10M',II,JJ); v10 = max(v10,0.001)
    u10 = interpolate_forcing(hour,nctime,tim1,nc,'U10M',II,JJ); u10 = max(u10,0.001)
    tcc = interpolate_forcing(hour,nctime,tim1,nc,'TCC',II,JJ); tcc = min(tcc,0.999)
    swr = interpolate_forcing(hour,nctime,tim1,nc,'SSR',II,JJ); swr = swr /3600. / 6.; swr = max(swr,0.0)
    
    nc.close()

    return tp,msl,t2m,d2m,v10,u10,tcc,swr

def interpolate_forcing(desired_hour,ncdf_time,time_index,ncdf_file,variable,Iindex,Jindex):

   pt1 = time_index - 1
   pt3 = time_index + 1
   if time_index == 0: pt1 = time_index
   if time_index == ncdf_time.shape[0]-1: pt3 = ncdf_time.shape[0]-1

   forcing = ncdf_file.variables[variable][pt1:pt3+1,Jindex,Iindex]
   interpolated = np.interp(desired_hour,ncdf_time[pt1:pt3+1],forcing)

   return interpolated
   

def prepare_era5_input(lat,lon,start,end, *args, **kwargs):

    if np.size(lat) == 1: 
        fixed = True
    else:
        fixed = False

    if kwargs.get('spinup'):
       f1name = 'meteo.spinup.dat'
       f2name = 'precip.spinup.dat'
       f3name = 'ssr.spinup.dat'
    else:
       f1name = 'meteo.argo.dat'
       f2name = 'precip.argo.dat'
       f3name = 'ssr.argo.dat'

    if fixed:
       starty = int(start[0:4]); startm = int(start[5:7]); startd = int(start[8:10])
       endy = int(end[0:4]); endm = int(end[5:7]); endd = int(end[8:10])
       datestart = datetime.datetime(starty,startm,startd) - relativedelta(hours=1) 
       dateend = datetime.datetime(endy,endm,endd) + relativedelta(hours=1)
    else:
       datestart = start
       dateend   = end

    t = 0
    with open(f1name,'w') as f, open(f2name,'w') as f2, open(f3name,'w') as f3:
       atm_prefix = 'ERA5'
  #     if datestart >= ( datetime.date.today() - relativedelta(months=6) ):
       if datestart.year >= 2023:
             atm_prefix='ECNC'
             url='./'+atm_prefix+'.GOTM.'+str(datestart.year)+'.nc'
             nc=NetCDFFile(url)
             nclat = nc.variables['lat'][:]
             nclon = nc.variables['lon'][:]
       else:
             url='./'+atm_prefix+'.GOTM.'+str(datestart.year)+'.nc'
             nc=NetCDFFile(url)
             nclat = nc.variables['latitude'][:]
             nclon = nc.variables['longitude'][:]
       target_year = datestart.year 

       while datestart <= dateend :

          if datestart.year != target_year:
             nc.close()
             #if datestart >= ( datetime.date.today() - relativedelta(months=6) ):
             if datestart.year >= 2023:
                atm_prefix='ECNC'
                url='./'+atm_prefix+'.GOTM.'+str(datestart.year)+'.nc'
                nc=NetCDFFile(url)
                nclat = nc.variables['lat'][:]
                nclon = nc.variables['lon'][:]
             else:
                url='./'+atm_prefix+'.GOTM.'+str(datestart.year)+'.nc'
                nc=NetCDFFile(url)
                nclat = nc.variables['latitude'][:]
                nclon = nc.variables['longitude'][:]
             target_year = datestart.year

          if fixed:    
             II = np.abs(lon-nclon).argmin()
             JJ = np.abs(lat-nclat).argmin()
             dist = distance_on_unit_sphere(lat,lon,nclat[JJ],nclon[II])
             if dist > 50. : print("Argo and forcing is: "+str(dist)+" km away "+str(datestart)+\
                                    "Argo lat/lon: "+str(lat[t])+"/"+str(lon[t])+\
                                    " forcing lat/lon: "+str(nclat[JJ])+"/"+str(nclon[II]))
          else:
             if t<lon.shape[0]:
              II = np.abs(lon[t]-nclon).argmin()
              JJ = np.abs(lat[t]-nclat).argmin()
              dist = distance_on_unit_sphere(lat[t],lon[t],nclat[JJ],nclon[II])
              if dist > 50. : print("Argo and forcing is: "+str(dist)+" km away "+str(datestart)+\
                                    "Argo lat/lon: "+str(lat[t])+"/"+str(lon[t])+\
                                    " forcing lat/lon: "+str(nclat[JJ])+"/"+str(nclon[II]))
             else:
              print('end date is within 2 profiles, will use the last coordinates')
              print('should not be a concern if this is for a couple of days (~48 warnings)')
          if atm_prefix == 'ECNC':
             tp,msl,t2m,d2m,v10,u10,tcc,ssr = get_ecnc(JJ,II,datestart)
          else:
             tp,msl,t2m,d2m,v10,u10,tcc,ssr = get_era5(JJ,II,datestart)
          
          f.write('{}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\t{:.8f}\n'\
                  .format(str(datestart),u10,v10,msl,t2m,d2m,tcc))
          f2.write('{}\t{:.8f}\n'.format(str(datestart),tp))
          f3.write('{}\t{:.8f}\n'.format(str(datestart),ssr))
          datestart = datestart + relativedelta(hours=1)
          t = t + 1

def get_bathy(lat,lon,*args, **kwargs):
    url = kwargs.get('url',None)
    if kwargs.get(url) is not None:
       nc=NetCDFFile(url)
       nclat = nc.variables['lat'][:]
       nclon = nc.variables['lon'][:]

       elevation = np.zeros((lat.shape[0]))
       for l in range(lat.shape[0]):
           II = np.abs(lon[l]-nclon).argmin()
           JJ = np.abs(lat[l]-nclat).argmin()
           elevation[l] = nc.variables['elevation'][JJ,II]
           #print(str(l+1)+' of '+str(lat.shape[0])+' depth = '+str(elevation[l])) 
           
       return elevation   
    else:
       print('You need to provide the netcdf file path. e.g. url=PATH')

def prepare_argo_only_relaxation(fld,prs,tim,lat,lon,n):
    array = np.zeros((n+1))
    tt = 0.
    la = 0.
    lo = 0.
    for t in range(tim.shape[0]):
        #print(t)
        #c = np.ma.MaskedArray.count(fld[t,:])
        #if c > n :
        #if c > 0 :
        
            prs1 = prs[t,~fld[t,:].mask]
            fld1 = fld[t,~fld[t,:].mask]
            #print(prs1.min())
            if prs1.max()>np.float(n) and prs1.min()<=15.:
              prs1 = np.ma.concatenate((0.,prs1),axis=None)
              fld1 = np.ma.concatenate((fld1[0],fld1[:]),axis=None)
              arr = np.interp(np.arange(n+1),prs1,fld1)
           #if arr.min() > -5:
              tt = np.ma.concatenate((tt,tim[t,0]),axis=None)
              la = np.ma.concatenate((la,lat[t]),axis=None)
              lo = np.ma.concatenate((lo,lon[t]),axis=None)
              array = np.vstack((array,arr))
    #print(array.shape)            
    array = array[1:,:]
    tt = tt[1:]
    la = la[1:]
    lo = lo[1:]

    return tt,-np.arange(n+1),array,la,lo

def merge_argo_woa(fld,name,time,since,lat,lon,depth):
    array=[]; arrayd=[]; arrayt=[]   
    with open(name+'.argo.dat','w') as f:

      for t in range(time.shape[0]):
        year  = int((datetime.datetime(since,1,1) + datetime.timedelta(days=time[t])).strftime('%Y'))
        month = int((datetime.datetime(since,1,1) + datetime.timedelta(days=time[t])).strftime('%m'))
        day   = int((datetime.datetime(since,1,1) + datetime.timedelta(days=time[t])).strftime('%d'))
        hours = int((datetime.datetime(since,1,1) + datetime.timedelta(days=time[t])).strftime('%H'))
        mins  = int((datetime.datetime(since,1,1) + datetime.timedelta(days=time[t])).strftime('%M'))
        print(str(datetime.datetime(year,month,day,hours,mins))+' step '+str(t)+" of "\
               +str(time.shape[0]-1)+" lat: "+str(lat[t])+" lon: "+str(lon[t]))
        if day >= 15.:
           prof1,d = get_woa(lat[t],lon[t],month,name)
           refday1=((datetime.datetime(year,month,15)-datetime.datetime(since,1,1)).days)
           if month == 12 :
              prof2,d = get_woa(lat[t],lon[t],1,name)
              refday2=((datetime.datetime(year+1,1,15)-datetime.datetime(since,1,1)).days)
           else:
              prof2,d = get_woa(lat[t],lon[t],month+1,name)
              refday2=((datetime.datetime(year,month+1,15)-datetime.datetime(since,1,1)).days)
        else:
           if month==1:
              prof1,d = get_woa(lat[t],lon[t],12,name)
              refday1=((datetime.datetime(year-1,12,15)-datetime.datetime(since,1,1)).days)
           else:
              prof1,d = get_woa(lat[t],lon[t],month-1,name)
              refday1=((datetime.datetime(year,month-1,15)-datetime.datetime(since,1,1)).days)
           prof2,d = get_woa(lat[t],lon[t],month,name)
           refday2=((datetime.datetime(year,month,15)-datetime.datetime(since,1,1)).days)
        
        woa = ((prof2-prof1)/(refday2-refday1+1.))*(time[t]-refday1) + prof1

        if not np.ma.MaskedArray.count(fld[t,:]) == 0:
           prs = -depth[t,~fld[t,:].mask]
           fld1 = fld[t,~fld[t,:].mask]
           if prs.max()>=-20.:
              if np.size(prs) <= 10 : print('sample number too low: '+str(np.size(prs))+' '+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2))
              if prs.min() >= -300. : 
                 print('sample depth too low: '+str(prs.min())+' '+str(year)+'-'+str(month).zfill(2)+'-'+str(day).zfill(2))
              else:
              #fld = np.interp(np.arange(-prs[0],-prs[-1],0.5),-prs,fld)
              #prs = -np.arange(-prs[0],-prs[-1],0.5)
              #prs[0]=0.

                 dfinal = np.concatenate((prs,d[d<prs.min()]),axis=None)
                 vfinal = np.concatenate((fld1,woa[d<prs.min()]),axis=None)

                 # print to file here
                 f.write('{} {} {}'.format(str(datetime.datetime(year,month,day,hours,mins)),str(dfinal.shape[0]),'2\n'))
                 for l in range(dfinal.shape[0]):
                    f.write('{} {}'.format(str(dfinal[l]),str(vfinal[l])+'\n'))

                 if len(array)==0:
                    array  = vfinal
                    arrayd = dfinal
                    arrayt = np.repeat(time[t],arrayd.shape[0])
                 else:
                    array  = np.ma.concatenate((array,vfinal),axis=None)
                    arrayd = np.ma.concatenate((arrayd,dfinal),axis=None)
                    arrayt = np.ma.concatenate((arrayt,np.repeat(time[t],arrayd.shape[0])),axis=None)

    return array,arrayd,arrayt

def vs_sat(time,lat,lon,number,*args, **kwargs):
    '''
    #how to:
    import argoinput
    nc = argoinput.read('FILENAME')
    time = argoinput.getvar(nc,'time',2010,2020)

    sincey  = int(nc['time'].units.split(' ')[2].split('-')[0])
    sincem  = int(nc['time'].units.split(' ')[2].split('-')[1])
    sinced  = int(nc['time'].units.split(' ')[2].split('-')[2])
    sinceh  = int(nc['time'].units.split(' ')[3].split(':')[0])
    sincemi = int(nc['time'].units.split(' ')[3].split(':')[1])
    sinces  = int(nc['time'].units.split(' ')[3].split(':')[2])
    sincedd = datetime.datetime(int(sincey),int(sincem),int(sinced),int(sinceh),int(sincemi),int(sinces))

    argoinput.vs_sat(time,15.0,70.0,'OUTPUTNAME',sinced=sincedd)
    '''
    sst_prefix = 'https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/'
    sst_suffix = '090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'
    bio_prefix = 'https://www.oceancolour.org/thredds/dodsC/cci/v5.0-release/geographic/daily/'
    chl_suffix = '/ESACCI-OC-L3S-CHLOR_A-MERGED-1D_DAILY_4km_GEO_PML_OCx-'
    kd_suffix  = '/ESACCI-OC-L3S-K_490-MERGED-1D_DAILY_4km_GEO_PML_KD490_Lee-'

    # get sst lat,lon
    sample_sst = sst_prefix + str(2014) + '/' + str(1).zfill(3) + '/' \
                      + str(2014) + str(1).zfill(2) + str(1).zfill(2) + sst_suffix
    nc = NetCDFFile(sample_sst)
    sstlat = nc.variables['lat'][:]
    sstlon = nc.variables['lon'][:]
    nc.close()
    # get chl lat,lon
    sample_chl = bio_prefix + 'chlor_a/' + str(2014) + chl_suffix + \
                           str(2014) + str(1).zfill(2) + str(1).zfill(2) + '-fv5.0.nc'
    nc = NetCDFFile(sample_chl)
    chllat = nc.variables['lat'][:]
    chllon = nc.variables['lon'][:]
    nc.close()



    xi = time[:,0]
    if np.size(lon) == 1:
       lon = np.repeat(lon,xi.shape[0])
       lat = np.repeat(lat,xi.shape[0])
    sinced = kwargs.get('sinced',None)
    if_model = kwargs.get('if_model',None)
    if sinced is not None:
       dates = list()

       count = 0
       DATA = {'dates':list(),'lat':[],'lon':[],'sst':[],'chl':[],'kd':[]}
       year_check = int(0); month_check = int(0); day_check = int(0); hour_check = int(0); min_check = int(0)       
       sincedd = sinced
       for t in range(xi.shape[0]):
           if if_model:
              year = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%Y'))
              month = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%m'))
              day   = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%d'))
              hours = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%H'))
              mins  = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%M'))
              nday  = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%j'))
           else:
              year = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%Y'))
              month = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%m'))
              day   = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%d'))
              hours = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%H'))
              mins  = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%M'))
              nday  = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%j'))
           #dates.append(datetime.datetime(year,month,day,hours,mins))

           if int(year) == year_check and int(month) == month_check and int(day) == day_check \
              and int(hours) == hour_check and int(mins) == min_check :

              pass
           else:
             year_check = int(year); month_check = int(month); day_check = int(day); hour_check = int(hours); min_check = int(mins)

             sst_name = sst_prefix + str(year) + '/' + str(nday).zfill(3) + '/' \
                      + str(year) + str(month).zfill(2) + str(day).zfill(2) + sst_suffix
             chl_name = bio_prefix + 'chlor_a/' + str(year) + chl_suffix + \
                           str(year) + str(month).zfill(2) + str(day).zfill(2) + '-fv5.0.nc'
             kd_name  = bio_prefix + 'kd/' + str(year) + kd_suffix + \
                           str(year) + str(month).zfill(2) + str(day).zfill(2) + '-fv5.0.nc'
             print(kd_name)
             # NOW SST
             II = np.abs(lon[t]-sstlon).argmin()
             JJ = np.abs(lat[t]-sstlat).argmin()
             pick = 15 # number of points
             radius = 1 # km

             try:
                nc = NetCDFFile(sst_name)
                satsst = nc.variables['analysed_sst'][0,JJ-pick:JJ+pick,II-pick:II+pick]
                satlatsub = nc.variables['lat'][JJ-pick:JJ+pick]
                satlonsub = nc.variables['lon'][II-pick:II+pick]
                nc.close()

                dist = np.zeros((pick*2,pick*2))
                for i in range(pick*2):
                  for j in range(pick*2):
                     dist[j,i] = distance_on_unit_sphere(lat[t],lon[t],satlatsub[j],satlonsub[i])

                sstsat = np.mean(satsst[dist<=radius]) - 273.15
                DATA['dates'].append(datetime.datetime(year,month,day,hours,mins))
                DATA['lat'] = np.ma.concatenate( (DATA['lat'],lat[t]), axis=None )
                DATA['lon'] = np.ma.concatenate( (DATA['lon'],lon[t]), axis=None )
                DATA['sst'] = np.ma.concatenate( (DATA['sst'],sstsat), axis=None )
             except Exception:
                pass

             # NOW CHL
             is_it_there=True
             try:
                nc = NetCDFFile(chl_name)
             except Exception as e :
                print('File does not exist : chlor_a : ' + \
                  str(year) + str(month).zfill(2) + str(day).zfill(2))
                is_it_there=False
             if is_it_there:
                II = np.abs(lon[t]-chllon).argmin()
                JJ = np.abs(lat[t]-chllat).argmin()
                pick = 5 # number of points
                radius = 5 # km
 
                nc = NetCDFFile(chl_name)
                satchl = nc.variables['chlor_a'][0,JJ-pick:JJ+pick,II-pick:II+pick]
                satlatsub = nc.variables['lat'][JJ-pick:JJ+pick]
                satlonsub = nc.variables['lon'][II-pick:II+pick]
                nc.close()

             is_it_there_kd=True
             try:
                nckd = NetCDFFile(kd_name)
             except Exception as e :
                is_it_there_kd=False

             if is_it_there_kd:
                kd490 = nckd.variables['kd_490'][0,JJ-pick:JJ+pick,II-pick:II+pick]
                nckd.close()

             dist = np.zeros((pick*2,pick*2))
             for i in range(pick*2):
                for j in range(pick*2):
                   dist[j,i] = distance_on_unit_sphere(lat[t],lon[t],satlatsub[j],satlonsub[i])

             chls = np.mean(satchl[dist<=radius]) # satellite chl is the mean within radius
             if is_it_there_kd:
                kd   = 1./np.mean(kd490[dist<=radius]) # mean satellite sample depth 
             if np.ma.is_masked(kd):
                kd = 10.
             else :
                kd = 10.

             DATA['chl'] = np.ma.concatenate( (DATA['chl'],chls), axis=None )
             DATA['kd']  = np.ma.concatenate( (DATA['kd'],kd), axis=None )

             print(DATA['dates'][count],DATA['lat'][count],DATA['lon'][count],DATA['sst'][count],DATA['chl'][count],DATA['kd'][count])
             count += 1
    f=open('satellite_sst_chl_kd_'+str(number)+'.pckl','wb')
    pickle.dump(DATA,f)
    f.close()
    
    return DATA
def getvar(nc,varib, *args, **kwargs):
    sincey = int(nc['time'].units.split(' ')[2].split('-')[0])
    sincem = int(nc['time'].units.split(' ')[2].split('-')[1])
    sinced = int(nc['time'].units.split(' ')[2].split('-')[2])


    time = nc.variables['time'][:]
    depth = nc.variables['z'][:,:,0,0]
    if varib == 'nuh':
       depth = nc.variables['zi'][:,:,0,0]

    timeall = np.zeros(( depth.shape[0],depth.shape[1] ))
    yearall = np.zeros(( depth.shape[0],depth.shape[1] ))
    year = np.zeros(( time.shape[0] ))

    dates = list()
    for t in range(time.shape[0]) :
       yearr = ( datetime.datetime(sincey, sincem, sinced) + datetime.timedelta(int(time[t]/60./60./24.)) )
       year[t] = int( yearr.strftime('%Y') )
       dates.append(yearr)

    first = kwargs.get('first',None)
    first = year[0] if first is None else first

    last = kwargs.get('last',None)
    last = year[-1] if last is None else last


    yindex = np.where((year >= first) & (year <= last)) 
    for d in range(depth.shape[1]) :
       timeall[:,d] = time

    if varib == 'time':
       var = np.copy(timeall)
    elif varib == 'date':
       var = np.copy(dates)
    elif varib == 'depth':
       var = np.copy(depth)
    elif varib == 'depthi':
       var = nc.variables['zi'][:,:,0,0]    
    elif varib == 'mld':
       var = nc.variables['mld_surf'][:,0,0]
       var = var1 + var2
    elif varib == 'chl' :
       var1 = nc.variables['ECO_diachl'][:,:,0,0]
       var2 = nc.variables['ECO_flachl'][:,:,0,0]
       var = var1 + var2
       if 'ECO_coccochl' in nc.variables.keys():
         var3 = nc.variables['ECO_coccochl'][:,:,0,0]
         var = var + var3
       if 'ECO_bg' in nc.variables.keys():
         var4 = nc.variables['ECO_bgchl'][:,:,0,0]
         var = var + var4
    elif varib == 'chl_norwecom' :
       var = nc.variables['ECO_chla'][:,:,0,0]
    elif varib == "nit_norwecom":
       var = nc.variables['ECO_nit'][:,:,0,0]
       var = var / 14.007
    elif varib == "pocs":
       var1 = nc.variables['ECO_microzoo'][:,:,0,0]
       var2 = nc.variables['ECO_det'][:,:,0,0]
       var3 = nc.variables['ECO_dia'][:,:,0,0]
       var4 = nc.variables['ECO_fla'][:,:,0,0]
       var = var1 + var2 + var3 + var4
       if 'ECO_cocco' in nc.variables.keys():
         var5 = nc.variables['ECO_cocco'][:,:,0,0]
         var = var + var5
    elif varib == "pocs_norwecom":
       var1 = nc.variables['ECO_mic'][:,:,0,0]
       var2 = nc.variables['ECO_det'][:,:,0,0]
       var3 = nc.variables['ECO_dia'][:,:,0,0]
       var4 = nc.variables['ECO_fla'][:,:,0,0]
       var = var1 + var2 + var3 + var4
       var = var/14.007*6.625*12.01
    elif varib == "poc":
       var1 = nc.variables['ECO_microzoo'][:,:,0,0]
       var2 = nc.variables['ECO_det'][:,:,0,0]
       var3 = nc.variables['ECO_dia'][:,:,0,0]
       var4 = nc.variables['ECO_fla'][:,:,0,0]
       var = var1 + var2 + var3 + var4
       if 'ECO_cocco' in nc.variables.keys():
         var5 = nc.variables['ECO_cocco'][:,:,0,0]
         var = var + var5
       var6 = nc.variables['ECO_mesozoo'][:,:,0,0]
       var = var + var6
    elif varib == 'oxy_norwecom' :
       var = nc.variables['ECO_oxy'][:,:,0,0]
       var = var/32.0
    else :
         var = nc.variables[varib][:,:,0,0]


    if varib == 'date':
       var = var[yindex]
    elif len(var.shape) == 2:
       var = var[yindex,:]
       var = var[0,:,:]
    elif len(var.shape) == 1:
       var = var[yindex]

    if varib == 'ECO_no3' :
       var = var / 12.01 / 6.625
    if varib == 'ECO_sil' :
       var = var / 12.01 / 6.625
    if varib == 'ECO_pho' : 
       var = var / 12.01 / 106.0
    if varib == 'ECO_primprod' :
       var = var * 86400.
    if varib == 'ECO_secprod' :
       var = var * 86400.

    return var

def get1D(argo_number,variable,d):
   number = int(argo_number)
   pres, since, time2d, dates, lon, lonint, lat, latint, timeint, temp, tempqc = readARGO(number, 'TEMP_ADJUSTED',all=True)
   if variable == "chl":
      for variable in ["CPHL_ADJUSTED"]:
         try:
            var, varqc = readARGO(number, variable)
         except Exception:
            pass
      var  = construct(var,varqc,[1,5,8])

   if variable == "oxy":
      for variable in ["DOXY_ADJUSTED","DOX2_ADJUSTED"]:
         try:
            var, varqc = readARGO(number, variable)
         except Exception:
            pass 
      var  = construct(var,varqc,[1,5,8])

   if variable == "nit":
      for variable in ["NTAW_ADJUSTED","NITRATE_ADJUSTED"]:
         try:
            var, varqc = readARGO(number, variable)
            var  = construct(var,varqc,[1,5,8])
         except Exception:
            pass 
      var  = construct(var,varqc,[1,2,5,8])

   if variable == "poc":
      for variable in ["CPHL_ADJUSTED"]:
         try:
            chl, chlqc = readARGO(number, variable)
         except Exception:
            pass
      chl  = construct(chl,chlqc,[1,5,8])
      for variable in ["BBP700_ADJUSTED"]:
         try:
            bbp, bbpqc = readARGO(number, variable)
         except Exception:
            pass
      bbp  = construct(bbp,bbpqc,[1,5,8])

      chlm = np.ma.masked_where( bbp.mask, chl)
      bbpm = np.ma.masked_where( chlm.mask, bbp) 

      k1 = 89.423; k2 = 0.1881; k3 = 0.7591; k4 = 0.1934; e1 = 1.626; e2 = 21.2; pocm = 33.4
      poc = k1 * bbpm**k2 * (chlm/bbpm)**k3 * (chlm/bbpm)**(k4*np.log10(bbpm))
      poc[poc<pocm] = e1 * poc[poc<pocm] - e2
      poc = np.ma.masked_where(poc>10000.,poc)
      var = np.ma.masked_where(poc<0.,poc)
 
   datea,deptha,vara=interpol(time2d,pres,since,var)
   var1d = []
#   print(-d,abs((abs(deptha)-d)).argmin() )
#   print(abs(deptha)-d)
   for t in range(len(datea)):
      #var1d.append( np.interp( -d, deptha, np.flipud(vara[:,t]) ))
      var1d.append( np.interp( d, abs(deptha), vara[:,t] ))
#   print(deptha.min(),deptha.max(),deptha.shape,vara.shape,np.array(var1d).max(),   np.flipud(vara)[(abs(deptha)-d).argmin(),t] ) 
   var1d = running_mean(var1d, 3)
   return datea,var1d

def get2D(argo_number,variable):
   number = int(argo_number)
   pres, since, time2d, dates, lon, lonint, lat, latint, timeint, temp, tempqc = readARGO(number, 'TEMP_ADJUSTED',all=True)
   if variable == "chl":
      for variable in ["CPHL_ADJUSTED"]:
         try:
            var, varqc = readARGO(number, variable)
         except Exception:
            pass
      var  = construct(var,varqc,[1,5,8])
   elif variable == "oxy":
      for variable in ["DOXY_ADJUSTED","DOX2_ADJUSTED"]:
         try:
            var, varqc = readARGO(number, variable)
         except Exception:
            pass 
      var  = construct(var,varqc,[1,5,8])
   elif variable == "nit":
      for variable in ["NTAW_ADJUSTED","NITRATE_ADJUSTED"]:
         try:
            var, varqc = readARGO(number, variable)
         except Exception:
            pass 
      var  = construct(var,varqc,[1,2,5,8])
   elif variable == "temp":
      var, varqc = readARGO(number, "TEMP_ADJUSTED")
      var  = construct(var,varqc,[1,5,8])
   elif variable == "sal":
      var, varqc = readARGO(number, "PSAL_ADJUSTED")
      var  = construct(var,varqc,[1,5,8])
   elif variable == "poc":
      for variable in ["CPHL_ADJUSTED"]:
         try:
            chl, chlqc = readARGO(number, variable)
         except Exception:
            pass
      chl  = construct(chl,chlqc,[1,5,8])
      for variable in ["BBP700_ADJUSTED"]:
         try:
            bbp, bbpqc = readARGO(number, variable)
         except Exception:
            pass
      bbp  = construct(bbp,bbpqc,[1,5,8])

      chlm = np.ma.masked_where( bbp.mask, chl)
      bbpm = np.ma.masked_where( chlm.mask, bbp) 

      k1 = 89.423; k2 = 0.1881; k3 = 0.7591; k4 = 0.1934; e1 = 1.626; e2 = 21.2; pocm = 33.4
      poc = k1 * bbpm**k2 * (chlm/bbpm)**k3 * (chlm/bbpm)**(k4*np.log10(bbpm))
      poc[poc<pocm] = e1 * poc[poc<pocm] - e2
      poc = np.ma.masked_where(poc>10000.,poc)
      var = np.ma.masked_where(poc<0.,poc)

   else:
      var, varqc = readARGO(number, variable)
      var  = construct(var,varqc,[1,5,8])

   datea,deptha,vara=interpol(time2d,pres,since,var)

   return datea,deptha,vara

def contourdata(x,y,z, *args, **kwargs):
    '''
    typical use:
    
    for ARGO:
    contourdata(timeargo,pres,flag_temp,nlevels=11,dbot=1000,sinced=[since,1,1,0,0,0],type='inferno') 

    for MODEL:
    contourdata(0,0,0,nlevels=11,dbot=1000,type='inferno',model=nc,variable='ECO_no3',sinced=0,minyear=2014,maxyear=2020)
    # default colormap = inferno, no need to use type='inferno' in that case
    '''
    import cmocean
    #plt.style.use('seaborn')

    model = kwargs.get('model',None)
    if model is not None:
       variable = kwargs.get('variable')
       sincey  = int(model['time'].units.split(' ')[2].split('-')[0])
       sincem  = int(model['time'].units.split(' ')[2].split('-')[1])
       sinced  = int(model['time'].units.split(' ')[2].split('-')[2])
       sinceh  = int(model['time'].units.split(' ')[3].split(':')[0])
       sincemi = int(model['time'].units.split(' ')[3].split(':')[1])
       sinces  = int(model['time'].units.split(' ')[3].split(':')[2])
       sincedd = datetime.datetime(int(sincey),int(sincem),int(sinced),int(sinceh),int(sincemi),int(sinces))
       ttt     = model['time'][:]
       first = int((datetime.datetime(sincey,sincem,sinced,sinceh,sincemi,sinces) + \
                  datetime.timedelta(seconds=ttt[0])).strftime('%Y'))
       last  = int((datetime.datetime(sincey,sincem,sinced,sinceh,sincemi,sinces) + \
                  datetime.timedelta(seconds=ttt[-1])).strftime('%Y'))
       z = getvar(model,variable)
       if variable == 'nuh':
          y = -(getvar(model,'depthi'))
       else:
          y = -(getvar(model,'depth'))
       x = getvar(model,'time')

       scale_factor = kwargs.get('scale_factor',None)
       if scale_factor is not None:
          z = z * float(scale_factor)

       #xnum = int(z.shape[0]/5.)
       #ynum = 300 

       #x = x.flatten()
       #y = y.flatten()
       #z = z.flatten()

       xi = x[:,0]
       yi = y[0,:]
       xm = x
       ym = y
       intp = z.T;

    else:
       xnum = int(z.shape[0]/1.)
       ynum = int(z.shape[1]/1.)

       x = x.flatten()
       y = y.flatten()
       z = z.flatten()

       x = x[~y.mask]
       z = z[~y.mask]
       y = y[~y.mask]

       x = x[~z.mask]
       y = y[~z.mask]
       z = z[~z.mask]

       xi = np.linspace(x.min(),x.max(),num=xnum)
       yi = np.linspace(y.min(),y.max(),num=ynum)
       xm,ym = np.meshgrid(xi,yi)

       intp = griddata( (x,y), z, (xm,ym), method='nearest')

    dtop = kwargs.get('dtop',None)
    dtop = (-yi).max() if dtop is None else -dtop
    dbot = kwargs.get('dbot',None)
    dbot = (-yi).min() if dbot is None else -dbot

    botloc = np.abs(yi-(-dbot)).argmin()
    toploc = np.abs(yi-(-dtop)).argmin()

    cmax = kwargs.get('cmax',None)
    cmax = intp[toploc:botloc,:].max() if cmax is None else cmax
    cmin = kwargs.get('cmin',None)
    cmin = intp[toploc:botloc,:].min() if cmin is None else cmin

    cmap = plt.get_cmap('inferno')
    cmaptype = kwargs.get('type',None)
    if cmaptype is not None:
       if cmaptype == 'algae':
          cmap = cmocean.cm.algae
       elif cmaptype == 'saline':
          cmap = cmocean.cm.saline
       elif cmaptype == 'solar':
          cmap = cmocean.cm.solar
       elif cmaptype == 'matter':
          cmap = cmocean.cm.matter
       elif cmaptype == 'amp':
          cmap = cmocean.cm.amp
       elif cmaptype == 'thermal':
          cmap = cmocean.cm.thermal
       elif cmaptype == 'inferno':
          cmap = plt.get_cmap('inferno')
       elif cmaptype == 'magma':
          cmap = plt.get_cmap('magma')
       elif cmaptype == 'spectral':
          cmap = plt.cm.get_cmap('Spectral_r')
       else:
          cmap = cmocean.cm.thermal
    #cmap = cmocean.cm.thermal if cmap is None else cmap #plt.get_cmap('Spectral_r')
    nlevels = kwargs.get('nlevels',None)
    nlevels = 10 if nlevels is None else nlevels
    levels = MaxNLocator(nbins=nlevels).tick_values( cmin,cmax )
#    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    sinced = kwargs.get('sinced',None)
    if sinced is not None:
       dates = list()
       if model is None:
          sincedd = datetime.datetime(int(sinced[0]),int(sinced[1]),int(sinced[2]),int(sinced[3]),int(sinced[4]),int(sinced[5]))
       for t in range(xi.shape[0]):
           if model is None:
              year = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%Y'))
              month = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%m'))
              day   = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%d'))
              hours = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%H'))
              mins  = int((sincedd + datetime.timedelta(days=xi[t])).strftime('%M'))
           else:
              year = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%Y'))
              month = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%m'))
              day   = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%d'))
              hours = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%H'))
              mins  = int((sincedd + datetime.timedelta(seconds=xi[t])).strftime('%M'))
           dates.append(datetime.datetime(year,month,day,hours,mins))

           nyears = dates[-1].year - dates[0].year + 1
           if nyears ==1 : 
              interval=1
           elif nyears == 2 or nyears == 3:
              interval = 3
           else :
              interval = 6
           minyear = kwargs.get('minyear',None)
           maxyear = kwargs.get('maxyear',None)
           if maxyear is None: 
              plotmin=datetime.datetime(dates[0].year,1,1)
              plotmax=datetime.datetime(dates[-1].year+1,1,1,0,0,0)
           else:
              plotmin=datetime.datetime(int(minyear),1,1)
              plotmax=datetime.datetime(int(maxyear)+1,1,1,0,0,0)   

    fig = plt.figure(figsize=(8.5,3.5),facecolor='w')
    ax  = fig.add_subplot(1,1,1)

    if sinced is None:
       cs = plt.contourf(xi,-yi,intp,levels=levels,cmap=cmap,extend='both')
    else:
       cs = plt.contourf(dates,-yi,intp,levels=levels,cmap=cmap,extend='both')
       ax.set_xticklabels(ax.get_xticks(), rotation = 45);
       fmt_half_year = mdates.MonthLocator(interval=interval)
       ax.xaxis.set_major_locator(fmt_half_year)
       fmt_month = mdates.MonthLocator()
       ax.xaxis.set_minor_locator(fmt_month)
       ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
       plt.xlim([plotmin,plotmax])
    cs.cmap.set_under(cmap(0))
    cs.cmap.set_over(cmap(255))
    plt.ylim([dbot,dtop])
    plt.colorbar()
    plt.tight_layout()

#    sprb = kwargs.get('spring_bloom')
#    if sprb:
#       fig = plt.figure(figsize=(7.5,3))
#       ax=fig.add_subplot(111)
#       if sinced is None:
#          for t in range(xi.shape[0]):
#              plt.plot(xi,intp[:,t].mean(),'.',color='k')

    return intp,dates,-yi

def scatterplot(model,obs,obsd,obst,variable, depth, *args, **kwargs):
# HOW TO
# scatterplot(model,flag_temp,pres,timeargo,'temp', 750.,colored='depth')

    minyear = kwargs.get('minyear',None)
    maxyear = kwargs.get('maxyear',None)
    if maxyear is None:
       sincey  = int(model['time'].units.split(' ')[2].split('-')[0])
       sincem  = int(model['time'].units.split(' ')[2].split('-')[1])
       sinced  = int(model['time'].units.split(' ')[2].split('-')[2])
       sinceh  = int(model['time'].units.split(' ')[3].split(':')[0])
       sincemi = int(model['time'].units.split(' ')[3].split(':')[1])
       sinces  = int(model['time'].units.split(' ')[3].split(':')[2])
       sincedd = datetime.datetime(int(sincey),int(sincem),int(sinced),int(sinceh),int(sincemi),int(sinces))
       ttt     = model['time'][:]
       first = int((datetime.datetime(sincey,sincem,sinced,sinceh,sincemi,sinces) + \
                  datetime.timedelta(seconds=ttt[0])).strftime('%Y'))
       last  = int((datetime.datetime(sincey,sincem,sinced,sinceh,sincemi,sinces) + \
                  datetime.timedelta(seconds=ttt[-1])).strftime('%Y'))
    else:
       first = maxyear
       last  = minyear
    z = getvar(model,variable)
    y = -(getvar(model,'depth'))
    x = getvar(model,'time')

    scale_factor = kwargs.get('scale_factor',None)
    if scale_factor is not None:
       z = z*scale_factor

    for p in range(x.shape[0]):
        timediff = ( datetime.datetime(sincey,sincem,sinced,sinceh,sincemi,sinces) + \
                  datetime.timedelta(seconds=x[p,0]) ) - \
                  datetime.datetime(1950,1,1)

        x[p,:] = timediff.days + timediff.seconds/86400.

    obst = obst.flatten()
    obsd = obsd.flatten()
    obs = obs.flatten()
    obst = obst[~obsd.mask]
    obs = obs[~obsd.mask]
    obsd = obsd[~obsd.mask]
    obst = obst[~obs.mask]
    obsd = obsd[~obs.mask]
    obs = obs[~obs.mask]

    intp = griddata( (obst,obsd), obs, (x,y), method='linear')

    monthly = kwargs.get('monthly',None)
    depthmin = kwargs.get('depthmin',None)       

    if monthly is not None:
          dummy = datetime.datetime(int(kwargs.get('styear')),int(kwargs.get('stmonth')),1) - datetime.datetime(1950,1,1)
          stattimemin = dummy.days + dummy.seconds/86400.
          dummy = datetime.datetime(int(kwargs.get('styear')),int(kwargs.get('stmonth')),1) + relativedelta(months=+1) - datetime.datetime(1950,1,1)  
          stattimemax = dummy.days + dummy.seconds/86400.
          index = np.logical_and( \
                     np.logical_and(y <= depth , y >= depthmin),\
                     np.logical_and(x < stattimemax , x >= stattimemin) \
                     )
    else:
          index = np.logical_and(y <= depth , y >= depthmin)

    plotx = z[index].flatten()
    ploty = intp[index].flatten()
    index2 = ~np.isnan(ploty)
    plotx = plotx[index2]
    ploty = ploty[index2]
  
    mask_below = kwargs.get('mask_below',None)
    if mask_below is not None:
       index5 = ploty >= mask_below
       plotx = plotx[index5]
       ploty = ploty[index5]
 
    xtitle = "Model"
    ytitle = "ARGO"
    if variable == 'chl':
       take_log = kwargs.get('take_log',None)
       if take_log is not None:
         plotx = np.log10(plotx)
         ploty = np.log10(ploty)
         xtitle = r"log$_{10}$(Model)"
         ytitle = r"log$_{10}$(ARGO)"
    index3 = np.isfinite(ploty)
    plotx = plotx[index3]
    ploty = ploty[index3]
    index4 = np.isfinite(plotx)
    plotx = plotx[index4]
    ploty = ploty[index4]


    stbias = bias(plotx,ploty)
    stcorr = correlation(plotx,ploty)
    strmse = rmse(plotx,ploty)
    stnstd = nstd(plotx,ploty)
    stcov  = covar(plotx,ploty)
    stvar  = np.var(plotx)
    p = plotx.shape[0]

    printplot = kwargs.get('printplot',None)
    if printplot is not None:
       limmin = np.min([plotx.min(),ploty.min()])
       limmin = np.max([limmin,-4])
       limmax = np.max([plotx.max(),ploty.max()])

       plt.style.use('seaborn-v0_8')
       colored = kwargs.get('colored',None)

       if colored is None:
          fig = plt.figure(figsize=(4,4))
       else:
          cmap = plt.get_cmap('Spectral_r')
          fig = plt.figure(figsize=(4.5,4))
       ax=fig.add_subplot(111)

       if colored == 'depth':
          plt.scatter(plotx,ploty,c=y[index][index2][index3][index4],s=1,cmap=cmap)
          plt.colorbar()
       elif colored == 'time':
          plt.scatter(plotx,ploty,c=x[index][index2][index3][index4],s=1,cmap=cmap)
          plt.colorbar()
       else:
          plt.scatter(plotx,ploty,s=1)

       ax.set_ylabel(ytitle,fontsize=12)
       ax.set_xlabel(xtitle,fontsize=12)
       plt.xlim([limmin,limmax])
       plt.ylim([limmin,limmax])

       plt.plot([limmin,limmax],[limmin,limmax])
       stats = kwargs.get('stats',None)
       if stats:
         ax.text(0.625, 0.175, '%bias: '+str(round(pc_bias(plotx,ploty),3)), fontsize=8,ha='left', va='top',\
          transform=fig.transFigure)
         ax.text(0.625, 0.2, 'bias: '+str(round(bias(plotx,ploty),3)), fontsize=8,ha='left', va='top',\
          transform=fig.transFigure)
         ax.text(0.625, 0.225, 'corr: '+str(round(correlation(plotx,ploty),3)), fontsize=8,ha='left', va='top',\
          transform=fig.transFigure)
         ax.text(0.625, 0.25, 'nstd: '+str(round(nstd(plotx,ploty),3)), fontsize=8,ha='left', va='top',\
          transform=fig.transFigure)
         ax.text(0.625, 0.275, 'rmse: '+str(round(rmse(plotx,ploty),3)), fontsize=8,ha='left', va='top',\
          transform=fig.transFigure)
    
       plt.tight_layout()

    if monthly is not None:
       return stbias, stcorr, strmse, stnstd, p
    if kwargs.get('get_stats'):
       return stbias, stcorr, strmse, stnstd, stvar, p, plotx, ploty
    if kwargs.get('get_interpolated'):
       return x,-y,z,intp

def statmap(model,obs,obsd,obst,variable, minyear, maxyear, depth, *args, **kwargs):

    increment = 10
    numbers = depth/increment

    nyear = maxyear - minyear + 1
    biasmap = np.zeros((int(numbers),nyear*12))
    corrmap = np.zeros((int(numbers),nyear*12))
    rmsemap = np.zeros((int(numbers),nyear*12))
    nstdmap = np.zeros((int(numbers),nyear*12))
    pmap    = np.zeros((int(numbers),nyear*12))
    dates   = list()



    t = 0
    for y in range(minyear,maxyear+1):
        for m in range(12):
            dates.append(datetime.datetime(y,m+1,15))
            for d in range(int(numbers)):

                print('stats for',y,m+1,d*increment,d*increment+increment)

                biasmap[d,t], corrmap[d,t], rmsemap[d,t], nstdmap[d,t], pmap[d,t] = \
                    scatterplot(model,obs,obsd,obst,variable, d*increment+increment,depthmin=d*increment,monthly=True,styear=y,stmonth=m+1)

            t = t + 1
    
    name = kwargs.get('name',None)
    varib = kwargs.get('varib',None)
    if variable is not None:
       if variable == 'chl':
         title_suffix = r" (mg m$^{-3}$)"
         title_prefix = 'Chlorophyll a '
         biasmin=-2.5; biasmax=2.5
         rmsemin=0.;rmsemax=2.5
       if variable == "ECO_oxy":
         title_suffix = r" (mmol m$^{-3}$)"
         title_prefix = 'Oxygen '
         biasmin=-25.; biasmax=25.
         rmsemin=0.;rmsemax=25. 
    else:
       title_suffix = ''
       title_prefix = '' 
    if name is not None:
       fsave = open('StatMAP_'+name+'.pckl','wb')
    else: 
       fsave = open('stat.pckl','wb')
    pickle.dump([biasmap,corrmap,rmsemap,nstdmap,pmap,dates],fsave)  
    fsave.close()

    if_plot = kwargs.get('if_plot',None)
    if if_plot:
       plt.style.use('seaborn-v0_8-notebook')
       import cmocean
       #cmap = plt.get_cmap('Spectral_r')
       cmap = cmocean.cm.balance
       depthh = -(np.arange(numbers)*increment+increment/2)
       fig = plt.figure(figsize=(6,9))

       ax1=fig.add_subplot(411)
       p1=plt.pcolormesh(dates,depthh,biasmap[0:int(numbers),:],cmap=cmap,shading='auto') 
       plt.xlim([datetime.datetime(minyear,1,1),datetime.datetime(maxyear+1,1,1)])
       plt.ylim([-depth,0])
       ax1.grid(which='major', axis='both',color='k')
       plt.clim([biasmin,biasmax])
       plt.title(title_prefix+'Bias'+title_suffix)
       plt.xticks(rotation = 25)
       plt.colorbar(p1)

       ax2=fig.add_subplot(412)
       p2=plt.pcolormesh(dates,depthh,corrmap[0:int(numbers),:],cmap=cmap,shading='auto')
       plt.xlim([datetime.datetime(minyear,1,1),datetime.datetime(maxyear+1,1,1)])
       plt.ylim([-depth,0])
       ax1.grid(which='major', axis='both',color='k')
       plt.clim([-1.,1.])
       plt.title(title_prefix+'Correlation') 
       plt.xticks(rotation = 25)
       plt.colorbar(p2)

       cmap = cmocean.cm.matter
       ax3=fig.add_subplot(413)
       p3=plt.pcolormesh(dates,depthh,rmsemap[0:int(numbers),:],cmap=cmap,shading='auto') 
       plt.xlim([datetime.datetime(minyear,1,1),datetime.datetime(maxyear+1,1,1)])
       plt.ylim([-depth,0])
       ax1.grid(which='major', axis='both',color='k')
       plt.clim([rmsemin,rmsemax])
       plt.title(title_prefix+'RMSE'+title_suffix)
       plt.xticks(rotation = 25)
       plt.colorbar(p3)

       ax4=fig.add_subplot(414)
       p4=plt.pcolormesh(dates,depthh,nstdmap[0:int(numbers),:],cmap=cmap,shading='auto') 
       plt.xlim([datetime.datetime(minyear,1,1),datetime.datetime(maxyear+1,1,1)])
       plt.ylim([-depth,0])
       ax1.grid(which='major', axis='both',color='k')
       plt.clim([0,2])
       plt.title(title_prefix+'Normalized Standard Dev.')
       plt.xticks(rotation = 25)
       plt.colorbar(p4)     

       plt.tight_layout()

    if name is not None:
       fig.canvas.print_figure('StatMAP_'+name+'.png',dpi=180)
    else: 
       fig.canvas.print_figure('StatMAP.png',dpi=180)

    return biasmap, corrmap, rmsemap, nstdmap, pmap, dates

def plot_argo_model_difference(x,y,m,a, *args, **kwargs):
    
    import cmocean
    cmap = cmocean.cm.balance

    plt.style.use('seaborn-v0_8-notebook')

    depth = y[0,:]
    dtop = kwargs.get('dtop',None)
    dtop = max(depth) if dtop is None else -dtop
    dbot = kwargs.get('dbot',None)
    dbot = min(depth) if dbot is None else -dbot

    first = kwargs.get('first',None)
    if kwargs.get('first'):
       year = int(first[0:4]); month = int(first[5:7]); day= int(first[8:10])
       plotmin = datetime.datetime(year,month,day)
    else:
       plotmin = x[0]

    last = kwargs.get('last',None)
    if kwargs.get('last'):
       year = int(last[0:4]); month = int(last[5:7]); day = int(last[8:10])
       plotmax = datetime.datetime(year,month,day)
    else:
       plotmax = x[-1]

    ind1 = np.logical_and(depth>=dbot,depth<=dtop)
    depth = depth[ind1]
    tx = list()
    for t in range(len(x)):
        tx.append(datetime.datetime(1950,1,1)+datetime.timedelta(x[t,0]))
    ind2 = np.logical_and(np.array(tx)>= plotmin, np.array(tx)<= plotmax)
    x = mdates.date2num(tx)
    x = x[ind2]
    fld = m - a
    fld = fld[:,ind1][ind2,:]


    cmax = kwargs.get('cmax',None)
    cmax = fld.max() if cmax is None else cmax
    cmin = kwargs.get('cmin',None)
    cmin = fld.min() if cmin is None else cmin

    fig = plt.figure(figsize=(10,3.5),facecolor='w')
    ax  = fig.add_subplot(1,1,1)
    ax.set_position([0.125,0.175,0.75,0.725])
    pmesh = plt.pcolormesh(x,depth,fld.T,cmap=cmap,shading='auto')

    ax.xaxis.axis_date()
    ax.yaxis.set_tick_params(labelsize='large')
    ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)

    plt.clim(cmin,cmax)
    plt.clim(cmin,cmax)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dbot,dtop])
    ax.grid(color='k', linestyle=':', linewidth=1)

#    if variable is not None:
#      plt.title(get_title(variable), loc='left', fontsize=15)

    ax.set_ylabel("Depth (m)",fontsize=15)
    title = kwargs.get('title',None)
    if title is not None: 
      plt.title("filename: "+title,loc="right", fontsize=10)

    cbaxes = fig.add_axes([0.9, 0.25, 0.025, 0.5])
    cb = plt.colorbar(pmesh, cax = cbaxes)
    cbaxes.yaxis.set_tick_params(labelsize='large')

def TSplot(temp,salt,pres, *args, **kwargs):
   smin = salt.min() - (0.01 * salt.min())
   smax = salt.max() + (0.01 * salt.max())
   tmin = temp.min() - (0.1 * temp.max())
   tmax = temp.max() + (0.1 * temp.max())

   xdim = int(round((smax-smin)/0.1+1,0))
   ydim = int(round((tmax-tmin)+1,0))

   mask_below = kwargs.get('mask_below',None)
   if mask_below is not None:
      temp = np.ma.masked_where(pres<mask_below, temp)
      salt = np.ma.masked_where(pres<mask_below, salt)
      pres = np.ma.masked_where(pres<mask_below, pres)
   # Create empty grid of zeros
   dens = np.zeros((ydim,xdim))
 
   # Create temp and salt vectors of appropiate dimensions
   ti = np.linspace(1,ydim-1,ydim)+tmin
   si = np.linspace(1,xdim-1,xdim)*0.1+smin
   
   for j in range(0,int(ydim)):
    for i in range(0, int(xdim)):
        dens[j,i]=sw.pden(si[i], ti[j], 0.0)

   dens = dens - 1000
   fig1 = plt.figure(facecolor='w')
   ax1 = fig1.add_subplot(111)
   CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
   plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level
 
   plt.scatter(salt,temp,c=pres,s=1)
   plt.colorbar()
 
   ax1.set_xlabel('Salinity')
   ax1.set_ylabel('Temperature (C)')


def plotModelvsArgo(ncname,chla,oxya,nita,poca,depth,dates, *args, **kwargs):
    import modelinput
    plt.style.use('seaborn-v0_8-notebook')
    #plt.style.use('Solarize_Light2')

    first = kwargs.get('first',None)
    if kwargs.get('first'):
       y = int(first[0:4]); m = int(first[5:7]); d = int(first[8:10])
       plotmin = datetime.datetime(y,m,d)
    else: 
       plotmin = dates[0]

    last = kwargs.get('last',None)
    if kwargs.get('last'):
       y = int(last[0:4]); m = int(last[5:7]); d = int(last[8:10])
       plotmax = datetime.datetime(y,m,d)
    else: 
       plotmax = dates[-1]

    nc  = modelinput.read(ncname)
    depthm  = modelinput.getvar(nc,"depth")
    time  = modelinput.getvar(nc,"time")
    norwecom = kwargs.get('norwecom',None)
    if kwargs.get('norwecom'):
      chlm = modelinput.getvar(nc,"ECO_chla")
      oxym = modelinput.getvar(nc,"ECO_oxy")/31.7
      nitm = modelinput.getvar(nc,"ECO_nit")/14.01
      detm = modelinput.getvar(nc,"ECO_det")/14.01 * 6.625 * 12.01
      flam = modelinput.getvar(nc,"ECO_fla")/14.01 * 6.625 * 12.01
      diam = modelinput.getvar(nc,"ECO_dia")/14.01 * 6.625 * 12.01
      micm = modelinput.getvar(nc,"ECO_mic")/14.01 * 6.625 * 12.01
      pocm = flam + diam + detm + micm
    else:       
      chlm = modelinput.getvar(nc,"chl")
      oxym = modelinput.getvar(nc,"oxy")
      nitm = modelinput.getvar(nc,"nit")
      pocm = modelinput.getvar(nc,"pocs")

    dtopc = kwargs.get('dtopc',None)
    dtopc = max(depth) if dtopc is None else -dtopc
    dbotc = kwargs.get('dbotc',None)
    dbotc = min(depth) if dbotc is None else -dbotc
    if chla is not None:
        ind1c = np.abs(dtopc - depth[:]).argmin() ; ind1c = min(ind1c, depth.shape[0]-1)
        ind2c = np.abs(dbotc - depth[:]).argmin() ; ind2c = max(ind2c, depth.shape[0])
        depthc = depth[ind1c:ind2c+1]
        fldchl = chla[ind1c:ind2c+1,:]
    ind1cm = np.abs(dtopc - depthm[0,:]).argmin() ; ind1cm = min(ind1cm + 1, depthm.shape[1]-1)
    ind2cm = np.abs(dbotc - depthm[0,:]).argmin() ; ind2cm = max(ind2cm - 1, 0)
    timecm = time[:,ind2cm:ind1cm+1]
    depthcm = depthm[:,ind2cm:ind1cm+1]
    fldchlm = chlm[:,ind2cm:ind1cm+1]
    cmaxc = kwargs.get('cmaxc',None)
    cmaxc = fldchlm.max() if cmaxc is None else cmaxc
    cminc = kwargs.get('cminc',None)
    cminc = fldchlm.min() if cminc is None else cminc


    dtopo = kwargs.get('dtopo',None)
    dtopo = max(depth) if dtopo is None else -dtopo
    dboto = kwargs.get('dboto',None)
    dboto = min(depth) if dboto is None else -dboto
    if oxya is not None:
        ind1o = np.abs(dtopo - depth[:]).argmin() #; ind1o = min(ind1o, depth.shape[0]-1)
        ind2o = np.abs(dboto - depth[:]).argmin() #; ind2o = max(ind2o, depth.shape[0])
        deptho = depth[ind1o:ind2o+1]
        fldoxy = oxya[ind1o:ind2o+1,:]
    ind1om = np.abs(dtopo - depthm[0,:]).argmin() ; ind1om = min(ind1om + 1, depthm.shape[1]-1)
    ind2om = np.abs(dboto - depthm[0,:]).argmin() ; ind2om = max(ind2om - 1, 0)
    timeom = time[:,ind2om:ind1om+1]
    depthom = depthm[:,ind2om:ind1om+1]
    fldoxym = oxym[:,ind2om:ind1om+1]
    cmaxo = kwargs.get('cmaxo',None)
    maxo = fldoxym.max() if cmaxo is None else cmaxo
    cmino = kwargs.get('cmino',None)
    cmino = fldoxym.min() if cmino is None else cmino


    dtopn = kwargs.get('dtopn',None)
    dtopn = max(depth) if dtopn is None else -dtopn
    dbotn = kwargs.get('dbotn',None)
    dbotn = min(depth) if dbotn is None else -dbotn
    if nita is not None:
        ind1n = np.abs(dtopn - depth[:]).argmin() #; ind1n = min(ind1n, depth.shape[0]-1)
        ind2n = np.abs(dbotn - depth[:]).argmin() #; ind2n = max(ind2n, depth.shape[0])
        depthn = depth[ind1n:ind2n+1]
        fldnit = nita[ind1n:ind2n+1,:]
    ind1nm = np.abs(dtopn - depthm[0,:]).argmin() ; ind1nm = min(ind1nm + 1, depthm.shape[1]-1)
    ind2nm = np.abs(dbotn - depthm[0,:]).argmin() ; ind2nm = max(ind2nm - 1, 0)
    timenm = time[:,ind2nm:ind1nm+1]
    depthnm = depthm[:,ind2nm:ind1nm+1]
    fldnitm = nitm[:,ind2nm:ind1nm+1]
    cmaxn = kwargs.get('cmaxn',None)
    cmaxn = fldnitm.max() if cmaxn is None else cmaxn
    cminn = kwargs.get('cminn',None)
    cminn = fldnitm.min() if cminn is None else cminn

    dtopp = kwargs.get('dtopp',None)
    dtopp = max(depth) if dtopp is None else -dtopp
    dbotp = kwargs.get('dbotp',None)
    dbotp = min(depth) if dbotp is None else -dbotp
    if poca is not None:
        ind1p = np.abs(dtopp - depth[:]).argmin() #; ind1n = min(ind1n, depth.shape[0]-1)
        ind2p = np.abs(dbotp - depth[:]).argmin() #; ind2n = max(ind2n, depth.shape[0])
        depthp = depth[ind1p:ind2p+1]
        fldpoc = poca[ind1p:ind2p+1,:]
    ind1pm = np.abs(dtopp - depthm[0,:]).argmin() ; ind1pm = min(ind1pm + 1, depthm.shape[1]-1)
    ind2pm = np.abs(dbotp - depthm[0,:]).argmin() ; ind2pm = max(ind2pm - 1, 0)
    timepm = time[:,ind2pm:ind1pm+1]
    depthpm = depthm[:,ind2pm:ind1pm+1]
    fldpocm = pocm[:,ind2pm:ind1pm+1]
    cmaxp = kwargs.get('cmaxp',None)
    cmaxp = fldpocm.max() if cmaxp is None else cmaxp
    cminp = kwargs.get('cminp',None)
    cminp = fldpocm.min() if cminp is None else cminp




    fig = plt.figure(figsize=(17,6),facecolor='w')

    if chla is not None:
       ax  = fig.add_subplot(2,4,1)

       cmap = get_colormap("chl")
       pmesh = plt.pcolormesh(dates,depthc,fldchl,cmap=cmap,shading='auto')
       ax.xaxis.axis_date()
       #ax.yaxis.set_tick_params(labelsize='large')
       #ax.xaxis.set_tick_params(labelsize='large')
       plt.xticks(rotation = 25)

       plt.clim(cminc,cmaxc)
       ax.set_xlim([plotmin,plotmax])
       ax.set_ylim([dbotc,dtopc])
       ax.grid(color='k', linestyle=':', linewidth=1)
       plt.title(get_title("chl"), loc='left', fontsize=15)
       #plt.xticks([])
       ax.set_ylabel("Depth (m)",fontsize=15)
       plt.colorbar()

    if oxya is not None:
       ax  = fig.add_subplot(2,4,2)

       cmap = get_colormap("oxy")
       pmesh = plt.pcolormesh(dates,deptho,fldoxy,cmap=cmap,shading='auto')
       ax.xaxis.axis_date()
       #ax.yaxis.set_tick_params(labelsize='large')
       #ax.xaxis.set_tick_params(labelsize='large')
       plt.xticks(rotation = 25)

       plt.clim(cmino,cmaxo)
       ax.set_xlim([plotmin,plotmax])
       ax.set_ylim([dboto,dtopo])
       ax.grid(color='k', linestyle=':', linewidth=1)
       plt.title(get_title("oxy"), loc='left', fontsize=15)
       #ax.set_ylabel("Depth (m)",fontsize=15)
       #plt.xticks([])
       plt.colorbar()

    if nita is not None:
       ax  = fig.add_subplot(2,4,3)

       cmap = get_colormap("nit")
       pmesh = plt.pcolormesh(dates,depthn,fldnit,cmap=cmap,shading='auto')
       ax.xaxis.axis_date()
       #ax.yaxis.set_tick_params(labelsize='large')
       #ax.xaxis.set_tick_params(labelsize='large')
       plt.xticks(rotation = 25)

       plt.clim(cminn,cmaxn)
       ax.set_xlim([plotmin,plotmax])
       ax.set_ylim([dbotn,dtopn])
       ax.grid(color='k', linestyle=':', linewidth=1)
       plt.title(get_title("nit"), loc='left', fontsize=15)
       #ax.set_ylabel("Depth (m)",fontsize=15)
       #plt.xticks([])
       plt.colorbar()

    if poca is not None:
       ax  = fig.add_subplot(2,4,4)

       cmap = get_colormap("poc")
       pmesh = plt.pcolormesh(dates,depthp,fldpoc,cmap=cmap,shading='auto')
       ax.xaxis.axis_date()
       #ax.yaxis.set_tick_params(labelsize='large')
       #ax.xaxis.set_tick_params(labelsize='large')
       plt.xticks(rotation = 25)

       plt.clim(cminp,cmaxp)
       ax.set_xlim([plotmin,plotmax])
       ax.set_ylim([dbotp,dtopp])
       ax.grid(color='k', linestyle=':', linewidth=1)
       plt.title(get_title("poc"), loc='left', fontsize=15)
       #ax.set_ylabel("Depth (m)",fontsize=15)
       #plt.xticks([])
       plt.colorbar()

    ax  = fig.add_subplot(2,4,5)
    cmap = get_colormap("chl")
    pmesh = plt.pcolormesh(timecm,depthcm,fldchlm,cmap=cmap,shading='auto')
    ax.xaxis.axis_date()
    #ax.yaxis.set_tick_params(labelsize='large')
    #ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)
    plt.clim(cminc,cmaxc)
    plt.clim(cminc,cmaxc)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dbotc,dtopc])
    ax.grid(color='k', linestyle=':', linewidth=1)
    ax.set_ylabel("Depth (m)",fontsize=15)
    plt.colorbar()

    ax  = fig.add_subplot(2,4,6)
    cmap = get_colormap("oxy")
    pmesh = plt.pcolormesh(timeom,depthom,fldoxym,cmap=cmap,shading='auto')
    ax.xaxis.axis_date()
    #ax.yaxis.set_tick_params(labelsize='large')
    #ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)
    plt.clim(cmino,cmaxo)
    plt.clim(cmino,cmaxo)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dboto,dtopo])
    ax.grid(color='k', linestyle=':', linewidth=1)
    #ax.set_ylabel("Depth (m)",fontsize=15)
    plt.colorbar()

    ax  = fig.add_subplot(2,4,7)
    cmap = get_colormap("nit")
    pmesh = plt.pcolormesh(timenm,depthnm,fldnitm,cmap=cmap,shading='auto')
    ax.xaxis.axis_date()
    #ax.yaxis.set_tick_params(labelsize='large')
    #ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)
    plt.clim(cminn,cmaxn)
    plt.clim(cminn,cmaxn)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dbotn,dtopn])
    ax.grid(color='k', linestyle=':', linewidth=1)
    #ax.set_ylabel("Depth (m)",fontsize=15)
    plt.colorbar()

    ax  = fig.add_subplot(2,4,8)
    cmap = get_colormap("poc")
    pmesh = plt.pcolormesh(timepm,depthpm,fldpocm,cmap=cmap,shading='auto')
    ax.xaxis.axis_date()
    #ax.yaxis.set_tick_params(labelsize='large')
    #ax.xaxis.set_tick_params(labelsize='large')
    plt.xticks(rotation = 25)
    plt.clim(cminp,cmaxp)
    plt.clim(cminp,cmaxp)
    ax.set_xlim([plotmin,plotmax])
    ax.set_ylim([dbotp,dtopp])
    ax.grid(color='k', linestyle=':', linewidth=1)
    #ax.set_ylabel("Depth (m)",fontsize=15)
    plt.colorbar()

    plt.tight_layout()

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    dummy = (cumsum[N:] - cumsum[:-N]) / float(N)
    dummy2 = list()
    dummy2.append(x[0])
    for t in range(len(dummy)):
       dummy2.append(dummy[t])
    dummy2.append(x[-1])
    return dummy2

def bindata(y,x,v,*args, **kwargs):
    # for now the support is only for upper 500 meters
    # default binning is 1 meter
    match = kwargs.get('matchArgo',None)
    if match is not None:
        xold=x.copy()
        x=list()
        vold = np.copy(v)
#        v = np.zeros((len(xold),vold.shape[1]))
        v=list()
        for t in range(len(match)):
            index = np.abs(match[t] - np.array(xold)).argmin()
            if np.abs(match[t] - xold[index]).days <= 1.:
                x.append(xold[index])
                v.append(vold[index, :])
                #v[t,:] = vold[index,:]
                #print(match[t],"  ",xold[index])
                #print(np.abs(match[t] - x[index]).days))
            # else:
            #     x.append(xold[index])
            #     v[t,:] = -9999.
        v = np.array(np.vstack(v))

    length = 1.
    binnedt = list()
    binned = np.zeros((len(x),int(2500.0/length))) - 9999.
    binnedd = np.zeros((int(2500.0/length))) 

    for t in range(len(x)):
        binnedt.append(datetime.datetime(x[t].year,x[t].month,x[t].day))
        #length = 1.
        dmin = 0.
        dmax = dmin + length
        for d in range(binned.shape[1]):
            v1d = v[t,:]
            y1d = y[t,:]
            index = np.logical_and(y1d >= dmin, y1d < dmax)
            try:
                binned[t,d] = np.nanmean(v1d[index])
                if any(index) != True:
                   binned[t,d] = -9999.0 
            except:
                binned[t,d] = -9999.0
            binnedd[d] = (dmax + dmin)/2.
        #    if dmax > 75.:
        #        length = 3.
        #    if dmax > 200.:
        #        length = 5.
            dmin = dmin + length
            dmax = dmax + length
    binned = np.ma.masked_where(binned < -9990.,binned)

    first = kwargs.get('first',None)
    if first is not None:
       ye = int(first[0:4]); mo = int(first[5:7]); da = int(first[8:10])
       index = np.array(binnedt) >= datetime.datetime(ye,mo,da)
       binned = binned[index,:]
       binnedt = np.array(binnedt)[index]
    last = kwargs.get('last',None)
    if last is not None:
       ye = int(last[0:4]); mo = int(last[5:7]); da = int(last[8:10])
       index = np.array(binnedt) <= datetime.datetime(ye,mo,da)
       binned = binned[index,:]
       binnedt = np.array(binnedt)[index]
    dtop = kwargs.get('dtop',None)
    if dtop is not None:
       index = binnedd >= dtop
       binned = binned[:,index]
       binnedd = binnedd[index]
    dbot = kwargs.get('dbot',None)
    if dbot is not None:
       index = binnedd <= dbot
       binned = binned[:,index]
       binnedd = binnedd[index]
#    print(y.shape, v.shape, len(x),binned.shape)
    return binnedt,binnedd,binned

def Koest23_modelB_700p(bbp,chla,alignment=1,bias='power',alpha=0.125):
    """
    Function to derive POC [mg m^-3] (conversion) using vertical profile of particulate
    backscattering at 700 nm [m^-1] and chlorophyll-a concentration [mg m^-3]
    based on global dataset described in Koestner et al. 2023.

    Koestner D, Stramski D and Reynolds RA (2024) Improved multivariable
    algorithms for estimating oceanic particulate organic carbon
    concentration from optical backscattering and chlorophyll-a measurements.
    Front. Mar. Sci. 10:1197953. doi: 10.3389/fmars.2023.1197953

    Inputs: column vector of particulate backscattering at 700 nm (bbp
    [m^-1]) and corresponding chlorophyll-a concentration (chla [mg m^-3])

    Optional Inputs: alignment factor for difference between sensors,
    HyrdoScats were used for model coefficients (0.9 is probably appropriate
    if using ECO sensors but may be as low as 0.75 from Erickson et al.
    2022); low POC bias correction function (linear or power); alpha for
    prediction intervals

    Outputs: POC and prediction interval of POC in [mg m^-3]

    """

    import numpy as np
    from scipy.stats import t



#    if len(np.shape(bbp)) != 1 or len(np.shape(chla)) != 1 or len(bbp) != len(chla):
#        raise ValueError("bbp and chla must be 1-dimensional arrays of the same length.")

    shape1 = int(bbp.shape[0])
    shape2 = int(bbp.shape[1])

    if len(bbp.shape) != 1:
        bbp = bbp.flatten()

    if len(chla.shape) != 1:
        chla = chla.flatten()

    if len(bbp) != len(chla):
        raise ValueError("bbp and chla must be of the same length.")

    if np.any(bbp < 0):
        bbp[bbp < 0] = 0

    if np.any(chla < 0):
        chla[chla < 0] = 0

    # Alignment adjustment for use with different bbp sensor
    bbp = bbp * alignment
    # "comp" computation
    comp = chla / bbp

    # added DK 20240129, if few comp values over 0 exist, fix minimum value to
    # be 10. This could be updated in future; 14 seems to be a more common
    # minimum. Not sure it makes a huge difference.


    if np.sum(comp > 0) > 10:
        comp[comp == 0] = np.min(comp[comp > 0])
    else:
        comp[comp == 0] = 10

    # fix large comp values  DK 20230902
    comp[comp > 2000] = 2000

    # model B coefficients from Koestner et al. 2022
    # k_modelB=[89.4229308852505; 0.188074098663738; 0.759083352712280; 0.193397796405782]; %2022 version
    # bias_modelB_coef=[1.63582494446266; -21.2269735484462]; %2022 version
    # New model B coefficients from Koestner et al. 2023

    k_modelB = np.array([52.8187501942431, 0.135288289126603, 0.884851394851513, 0.226810214797258])

    # calculate POC from multivariable model
    pocB = k_modelB[0] * (bbp ** k_modelB[1]) * (comp ** k_modelB[2]) * (comp ** (k_modelB[3] * np.log10(bbp)))

    # continious bias correction for modeled POC < bias_lim
    if bias[:3] == 'pow':
        bias_modelB_coef = np.array([1.46918996207386, -0.734453171035830])
        bias_lim = 36.8
        pocB[pocB < bias_lim] = (pocB[pocB < bias_lim] ** bias_modelB_coef[0]) * 10 ** bias_modelB_coef[1]
    elif bias[:3] == 'lin':
        bias_modelB_coef = np.array([1.42263442605374, -14.7264359822160])
        bias_lim = 34.8
        pocB[pocB < bias_lim] = (pocB[pocB < bias_lim] * bias_modelB_coef[0]) + bias_modelB_coef[1]

    # prediction interval process
    S = np.array([[0.0248767643354040, 0.00967335656242312, -0.0105363562351115, -0.00404078515530673],
                  [0.00967335656242312, 0.00405822189121109, -0.00402605221279587, -0.00166729112989270],
                  [-0.0105363562351115, -0.00402605221279587, 0.00483637706215940, 0.00180593541611042],
                  [-0.00404078515530673, -0.00166729112989270, 0.00180593541611042, 0.000727608368784931]])
    mse = 0.0311478519483073
    df = len(bbp) - 4 - 2
    t_value = t.ppf(1 - alpha, df)


    x = np.column_stack([np.ones_like(bbp), np.log10(bbp), np.log10(comp), np.log10(comp) * np.log10(bbp)])
    e = np.zeros_like(bbp)  # Initialize prediction bounds
    for i in range(len(bbp)):
        e[i] = t_value * np.sqrt(mse + np.dot(x[i], np.dot(S, x[i])))

    # calculate prediction bounds log10(y +/- e), where e is prediction bound defined
    # as e = t * sqrt(MSE + xSx); x is the row vector of the design
    # matrix for new predictors ( i.e., [1 log10(bbp) log10(comp)
    # [log10(comp)*log10(bbp)]; % note that e will be in logspace
    # x=[ones(size(bbp)) log10(bbp) log10(comp) log10(comp).*log10(bbp)]


    # Rename final data and convert to mg m-3
    POC = pocB
    POC[POC < 0] = 0

    pi = np.column_stack((np.log10(pocB) - e, np.log10(pocB) + e))
    PI = 10 ** pi

    POC = POC.reshape((shape1, shape2))

    return POC,PI,comp
