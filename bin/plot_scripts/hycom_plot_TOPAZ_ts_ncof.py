#/busr/bin/env python
from netCDF4 import Dataset
import modeltools.hycom
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile.abfile as abf
import numpy as np
import logging
import datetime
import re
import scipy.interpolate
import os.path
import matplotlib.pyplot as plt
import mod_hyc2plot
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

"""
#usage
      python plot_scripts/hycom_plot_TOPAZ_ts_ncof.py --clim=-3,3  temp  1 ../data/monthly/archm.2009_{01..12}.a --filename2 /OSTIA_Regrid/TOPAZ_ncof_sst_2009_{01..12}.nc
"""
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
    fig = matplotlib.pyplot.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.
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


def open_file(myfile0,filetype,fieldname,fieldlevel,datetime1=None,datetime2=None,vector="",
      idm=None,
      jdm=None) :

   logger.info("Now processing  %s"%myfile0)
   m=re.match("(.*)\.[ab]",myfile0)
   if m :
      myfile=m.group(1)
   else :
      myfile=myfile0

   ab2=None
   rdtimes=[]
   if filetype == "archive" :
      ab = abf.ABFileArchv(myfile,"r")
      n_intloop=1

   elif filetype == "regional.depth" :
      ab = abf.ABFileBathy(myfile,"r",idm=idm,jdm=jdm)
      n_intloop=1
   elif filetype == "forcing" :
      ab = abf.ABFileForcing(myfile,"r",idm=idm,jdm=jdm)
      if vector :
         file2=myfile.replace(fieldname,vector)
         logger.info("Opening file %s for vector component nr 2"%file2)
         ab2=abf.ABFileForcing(file2,"r",idm=idm,jdm=jdm)
      if datetime1 is None or datetime2 is None :
         raise NameError("datetime1 and datetime2 must be specified when plotting forcing files")
      else :
         iday1,ihour1,isec1 = modeltools.hycom.datetime_to_ordinal(datetime1,3)
         rdtime1 = modeltools.hycom.dayfor(datetime1.year,iday1,ihour1,3)
         #
         iday2,ihour2,isec2 = modeltools.hycom.datetime_to_ordinal(datetime2,3)
         rdtime2 = modeltools.hycom.dayfor(datetime2.year,iday2,ihour2,3)
         rdtimes=sorted([elem for elem in ab.field_times if elem >rdtime1 and elem < rdtime2])
         n_intloop=len(rdtimes)
   else :
      raise NotImplementedError("Filetype %s not implemented"%filetype)
   # Check that fieldname is actually in file
   if fieldname not in  ab.fieldnames :
      logger.error("Unknown field %s at level %d"%(fieldname,fieldlevel))
      logger.error("Available fields : %s"%(" ".join(ab.fieldnames)))
      raise ValueError("Unknown field %s at level %d"%(fieldname,fieldlevel))

   return n_intloop,ab,ab2,rdtimes

def main(myfiles,fieldname,fieldlevel,
      idm=None,
      jdm=None,
      clim=None,
      filetype="archive",
      window=None,
      cmap="jet",
      datetime1=None,
      datetime2=None,
      vector="",
      tokml=False,
      masklim=None,
      filename2='',
      filename5='',
      dpi=180) :


   LinDic=mod_hyc2plot.cmap_dict('sawtooth_fc100.txt')
   if 'temp' in fieldname:
      cmap= matplotlib.colors.LinearSegmentedColormap('my_colormap',LinDic)
   else:
      cmap=matplotlib.pyplot.get_cmap("jet")
   if tokml :
      ab = abf.ABFileGrid("regional.grid","r")
      plon=ab.read_field("plon")
      plat=ab.read_field("plat")
      ab.close()

   
   TP5grid='/cluster/work/users/achoth/TP5a0.06/topo'
   ab = abf.ABFileGrid("regional.grid","r")
   plon=ab.read_field("plon")
   plat=ab.read_field("plat")
   scpx=ab.read_field("scpx")
   scpy=ab.read_field("scpy")
   target_lonlats=[plon,plat]
   abdpth = abf.ABFileBathy('regional.depth',"r",idm=ab.idm,jdm=ab.jdm)
   mdpth=abdpth.read_field('depth')
   maskd=mdpth.data
   maskd[maskd>1e29]=np.nan
   #Region_mask=True
   Region_mask=False
   if Region_mask:
      maskd[plat>70]=np.nan
      #maskd[plat<50]=np.nan
      maskd[plon>20]=np.nan
      maskd[plon<-30]=np.nan

   Nordic_mask=maskd   

   proj=ccrs.Stereographic(central_latitude=90.0,central_longitude=-40.0)
   pxy = proj.transform_points(ccrs.PlateCarree(), plon, plat)
   px=pxy[:,:,0]
   py=pxy[:,:,1]
   x,y=np.meshgrid(np.arange(plon.shape[1]),np.arange(plon.shape[0]))


   if vector :
      logger.info("Vector component 1:%s"%fieldname)
      logger.info("Vector component 2:%s"%vector) 

   #---------------first read and compute clim  
   Err_map=1
   freezp=-1.8
   freezp=-2.0
   sum_fld1=maskd
   sum_fld1[~np.isnan(sum_fld1)]=0.0
   Clim_arr=np.zeros((plon.shape[0],plon.shape[1],12))
   if 'tem' in fieldname:
      counter=0
      rlxfile0="/cluster/work/users/achoth/TP5a0.06/relax/050/relax_tem.a"
      rlx_afile = abf.AFile(ab.idm,ab.jdm,rlxfile0,"r")
      lyr=fieldlevel
      record_num=1
      record_var=record_num-1 
      fld = rlx_afile.read_record(record_var)
      print('mn,mx  data='),print(fld.min()),print(fld.max())
      kdm=50
      dt_clim=np.zeros(12)
      for mnth in range(12) :
          fld1=rlx_afile.read_record(mnth*kdm+lyr-1)
          print('record, mn,mx  data='), print(kdm*mnth, fld1.min(),fld1.max())
          print('record, mn,mx  data='), print(kdm*mnth, fld1.min(),fld1.max())
          # Intloop used to read more fields in one file. Only for forcing for now
          dt_clim[mnth]=mod_hyc2plot.spatiomean(fld1,maskd)
          sum_fld1=sum_fld1+fld1
          Clim_arr[:,:,mnth]=fld1[:,:]
          counter=counter+1
          print('counter=',counter)
          del  fld1
      Clim_Avg=sum_fld1/counter
      del  sum_fld1

   #--------------- 
   # compute for TP5 files
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   counter=0
   file_count=0
   sum_fld1=maskd
   sum_fld1[~np.isnan(sum_fld1)]=0.0
   if filename5:
       dt_5=np.zeros(len(filename5))
       diff_dt_5=np.zeros(len(filename5))
       rmse_dt_5=np.zeros(len(filename5))
       SSTost_arr=np.zeros((plon.shape[0],plon.shape[1],len(filename5)))
       Labl5="filename5"
       yyyy1=filename5[0][-10:-6]
       print("filename5[0]=",filename5[0][-10:-6])
       print("yyy1=", yyyy1)
       tid_5=np.array([datetime.datetime(int(yyyy1), 1, 15) + \
           relativedelta(months=i) for i in range(len(filename5))])
       Labl5="OSTIA SST"
       print("Labl5=", Labl5)
       for ncfile0 in filename5 :
             logger.info("Now processing  %s"%ncfile0)
             fh = Dataset(ncfile0, mode='r')
             fld = fh.variables['analysed_sst'][:]
             print("fld.shpe", fld.shape)
             fh.close()
             fld=fld[0,:,:]
             SSTost_arr[:,:,file_count]=fld[:,:]

             # Create scalar field for vectors
             J,I=np.meshgrid(np.arange(fld.shape[0]),np.arange(fld.shape[1]))
             print('mn,mx  data=',fld.min(),fld.max())
             print('mn,mx  data=',fld.min(),fld.max())
             sum_fld1=sum_fld1+fld
             dt_5[counter]=mod_hyc2plot.spatiomean(fld,Nordic_mask)
             diff_dt_5[counter]=mod_hyc2plot.spatiomean(fld[:,:]-Clim_arr[:,:,counter%12],Nordic_mask)
             print("diff_dt_5[counter]=",diff_dt_5[counter])
             rmse_dt_5[counter]=np.sqrt(mod_hyc2plot.spatiomean((fld[:,:]-Clim_arr[:,:,counter%12])**2,Nordic_mask))
             # Apply mask if requested
             counter=counter+1
             file_count=file_count+1
             del  fld

       print('Computing the avearge of file_counter= ', file_count, 'counter=',counter,'(counter-1)%12=', (counter-1)%12)
       if file_count> 0:
          Ncof_Avg=sum_fld1/file_count
          Clim_Avg=Ncof_Avg
       if Err_map:
          cmap=cmocean.cm.balance
          fld_diff=Ncof_Avg - Clim_Avg
          print("fld_diff.shape=", fld_diff.shape)

       P=ax.pcolormesh(x[J,I],y[J,I],fld_diff[J,I],cmap=cmap)
       if 'temp' in fieldname:
        P1=plt.contour(x[J,I],y[J,I],fld_diff[J,I],levels=[-1.,1,4.0,8],colors=('w',),linestyles=('-',),linewidths=(1.5,))
        matplotlib.pyplot.clabel(P1, fmt = '%2.1d', colors = 'w', fontsize=10)


       # Print figure.
       aspect = 40
       pad_fraction = 0.25
       divider = make_axes_locatable(ax)
       width = axes_size.AxesY(ax, aspect=1./aspect)
       pad = axes_size.Fraction(pad_fraction, width)
       cax = divider.append_axes("right", size=width, pad=pad)
       cb=ax.figure.colorbar(P,cax=cax,extend='both')
       if clim is not None : P.set_clim(clim)
       ax.set_title('Diff: OSTIA SST - CLIM SST '+"%s:%s(%d)"%(ncfile0,'analysed_sst',0))
       # Print figure.
       fnamepng_template=filename5[0][1:8].replace("/",'')+yyyy1+"_NCOF_SST_%s_%d_%03d_iii%03d_Avg.png"
       if Region_mask:
          fnamepng_template='Region_'+yyyy1+'_'+filename5[0][1:8].replace("/",'') \
                +"_NCOF_SST_%s_%d_%03d_iii%03d_Avg.png"
       fnamepng=fnamepng_template%('analysed_sst',0,counter,file_count)
       logger.info("output in  %s"%fnamepng)
       figure.canvas.print_figure(fnamepng,bbox_inches='tight',dpi=dpi)
       ax.clear()
       cb.remove()
       datmen=np.nanmean(fld_diff)
       print('-----------mean error data=',  datmen)
       del sum_fld1


   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   ax.set_facecolor('xkcd:gray')
   onemm=9.806
   counter=0
   file_count=0
   sum_fld1=maskd
   dt_cnl=np.zeros(len(myfiles))
   diff_dt_cnl=np.zeros(len(myfiles))
   rmse_dt_cnl=np.zeros(len(myfiles))
   Labl1="Model SST"
   if "SPRBAS_0" in myfiles[0]:
      Labl1="Model: prsbas=0 "
   if filename2:
      dt_2=np.zeros(len(filename2))
      diff_dt_2=np.zeros(len(filename2))
      rmse_dt_2=np.zeros(len(filename2))
      yyyy1=filename2[0][-10:-6]
      print("filename2[0]=",filename2[0][-10:-6])
      print("filename2[0]=",filename2[0])
      print("yyy1=", yyyy1)
      tid_2=np.array([datetime.datetime(int(yyyy1), 1, 15) \
         + relativedelta(months=i) for i in range(len(filename2))])
      Labl2="filename2"
      Labl2="Corrected"
      if "erai" in filename2[0]:
         Labl2="Model: prsbas=1e5 "

   yyyy1cnt=myfiles[0][-9:-5]
   print("myfiles[0]=",myfiles[0][-9:-5])
   print("myfiles[0]=",myfiles[0])
   print("yyy1cnt=", yyyy1cnt)
   base = datetime.datetime(int(yyyy1cnt), 1, 15)
   tid=np.array([base + relativedelta(months=i) for i in range(len(myfiles))])
   if len(myfiles)==36:
      base = datetime.datetime(int(yyyy1cnt), 1, 15)
      tid=np.array([base + relativedelta(months=i) for i in range(len(myfiles))])
   #tid=mdates.date2num(tidstr)
   nmexp=1
   if filename2:
      nmexp=nmexp+1
   print('processing data from No runs ==##############>>>>>>>>>>>>>>>>>>>>>>>', nmexp)
   for iii in range(nmexp):
      counter=0
      file_count=0
      sum_fld1=maskd
      sum_fld1[~np.isnan(sum_fld1)]=0.0
      if iii==1 and filename2:
         myfiles=filename2
      else:
         logger.info(">>>>>--------------------------Processing the first files=  %d<<<<"%iii)

      logger.info(">>>>>--------------------------Processing the first files=  %d<<<<"%iii)
      for myfile0 in myfiles :
         # Open files, and return some useful stuff.
         # ab2 i used in case of vector
         # rdtimes is used for plotting forcing fields
         n_intloop,ab,ab2,rdtimes = open_file(myfile0,filetype,fieldname,fieldlevel,\
               datetime1=datetime1,datetime2=datetime2,vector=vector,idm=idm,jdm=jdm)
         # Intloop used to read more fields in one file. Only for forcing for now
         for i_intloop in range(n_intloop) :
            # Read ab file of different types
            if filetype == "archive" :
               fld1 = ab.read_field(fieldname,fieldlevel)
            elif filetype == "forcing" :
               fld1 = ab.read_field(fieldname,rdtimes[i_intloop])
               if vector :fld2=ab2.read_field(vector,rdtimes[i_intloop])
               logger.info("Processing time %.2f"%rdtimes[i_intloop])
            else :
               raise NotImplementedError("Filetype %s not implemented"%filetype)
            if not window :
               J,I=np.meshgrid(np.arange(fld1.shape[0]),np.arange(fld1.shape[1]))
            else :
               J,I=np.meshgrid( np.arange(window[1],window[3]),np.arange(window[0],window[2]))
            # Create scalar field for vectors
            if vector : 
               fld = np.sqrt(fld1**2+fld2**2)
            else :
               fld=fld1

            print('---------mn,mx  data=',fld.min(),fld.max())

            print('---------mn,mx  data=',fld.min(),fld.max())
            sum_fld1=sum_fld1+fld

            if iii==0:
               dt_cnl[counter]=mod_hyc2plot.spatiomean(fld,Nordic_mask)
               diff_dt_cnl[counter]=mod_hyc2plot.spatiomean(fld[:,:]-SSTost_arr[:,:,counter],Nordic_mask)
               rmse_dt_cnl[counter]=np.sqrt(mod_hyc2plot.spatiomean((fld[:,:]-SSTost_arr[:,:,counter])**2,Nordic_mask))
               Labl=Labl1
            if iii==1 and filename2:
               dt_2[counter]=mod_hyc2plot.spatiomean(fld,Nordic_mask)
               diff_dt_2[counter]=mod_hyc2plot.spatiomean(fld[:,:]-SSTost_arr[:,:,counter],Nordic_mask)
               rmse_dt_2[counter]=np.sqrt(mod_hyc2plot.spatiomean((fld[:,:]-SSTost_arr[:,:,counter])**2,Nordic_mask))
               Labl=Labl2
            # Apply mask if requested
            counter=counter+1
            file_count=file_count+1
            del fld
         # End i_intloop
      print('Computing the avearge of file_counter= ', file_count, 'counter=',counter)
      if file_count> 0:
         fld_Avg=sum_fld1/file_count
      if Err_map:
         cmap=cmocean.cm.balance
         fld_diff=fld_Avg -Clim_Avg

      if fieldname=='k.e.' :
            P=ax.pcolormesh(x[J,I],y[J,I],np.log10(fld_Avg[J,I]),cmap=cmap)
      elif fieldname=='srfhgt' :
            P=ax.pcolormesh(x[J,I],y[J,I],(fld_Avg[J,I]/onemm),cmap=cmap)
      else :
            P=ax.pcolormesh(x[J,I],y[J,I],fld_diff[J,I],cmap=cmap)
            
      if 'temp' in fieldname:
        P1=ax.contour(x[J,I],y[J,I],fld_diff[J,I],levels=[-1.,1,4.0,8],colors=('w',),linestyles=('-',),linewidths=(1.5,))
        matplotlib.pyplot.clabel(P1, fmt = '%2.1d', colors = 'w', fontsize=10) #contour line labels

       ##
      # Print figure.
      ax.set_facecolor('xkcd:gray')
      aspect = 40
      pad_fraction = 0.25
      divider = make_axes_locatable(ax)
      width = axes_size.AxesY(ax, aspect=1./aspect)
      pad = axes_size.Fraction(pad_fraction, width)
      cax = divider.append_axes("right", size=width, pad=pad)
      cb=ax.figure.colorbar(P,cax=cax,extend='both')
      if clim is not None : P.set_clim(clim)
      #ax.set_title('Diff: Model SST - Ostia SST :'+myfiles[0])
      ax.set_title(Labl + '  :( Model SST - Ostia SST )')
      # Print figure.
      fnamepng_template=myfiles[0][-16:-5].replace("/",'') +"_Avg_TP2_%s_%d_%03d_iii%03d.png"
      if Region_mask:
         fnamepng_template='Region_'+yyyy1cnt+myfiles[0][7:8].replace("/",'') \
       +"_Avg_TP2_%s_%d_%03d_iii%03d.png"
      fnamepng=fnamepng_template%(fieldname,fieldlevel,counter,iii)
      logger.info("output in  %s"%fnamepng)
      figure.canvas.print_figure(fnamepng,bbox_inches='tight',dpi=dpi)
      ax.clear()
      cb.remove()
      datmen=np.nanmean(fld_diff)
      spatiodatmen=mod_hyc2plot.spatiomean(fld_diff,Nordic_mask)
      print('-----------mean diff data, spatio=', datmen,spatiodatmen)
      del sum_fld1
      #---------------------------------------


   #tid_clim=tid[::31]+14
   print('tid len=', tid.shape)
   if filename2:
      print('dt_2=', dt_2.shape)
   tid_clim=np.array([base + relativedelta(months=i) for i in range(12)])
   #figure, ax = matplotlib.pyplot.figure()
   figure, ax = plt.subplots()
   rpt=len(dt_cnl)/12
   dt_clim_cat=dt_clim
   for ii in range(int(rpt-1)):
      print("concatenate ")
      dt_clim_cat=np.concatenate([dt_clim_cat,dt_clim])
      
     
   years = YearLocator()   # every year
   months = MonthLocator()  # every month
   yearsFmt = DateFormatter('%Y')
   #ax=figure.add_subplot(111)
   nplts=1  
   ax.plot_date(tid, dt_cnl, '-o',color='g',ms=3,  label=Labl1)
   if filename2:
      ax.plot_date(tid_2, dt_2, '-v',color='orange',ms=3,  label=Labl2)
      nplts=nplts+1  
   if filename5:
      ax.plot_date(tid_5, dt_5,'--', color='m', label=Labl5)
      nplts=nplts+1  
   if 'tem' in fieldname:
      ax.plot_date(tid[0:len(dt_cnl)], dt_clim_cat[:],':' ,color='black', label='WOA18-Clim.')
      nplts=nplts+1  
   ax.xaxis.set_major_locator(years)
   ax.xaxis.set_major_formatter(yearsFmt)
   ax.xaxis.set_minor_locator(months)
   ax.autoscale_view()
   # format the coords message box
   def price(x):
       return '$%1.2f' % x
   ax.fmt_xdata = DateFormatter('%Y-%m-%d')
   ax.fmt_ydata = price
   ax.grid(True)
   figure.autofmt_xdate()
   legend = plt.legend(loc='upper left',fontsize=8)
   plt.title("Area-averaged: %s(%d)"%(fieldname,fieldlevel))
   #plt.xlabel('dayes')
   plt.ylabel("%s(%d)"%(fieldname,fieldlevel))
   #plt.title('Pakistan India Population till 2007')
   ts_fil="time_series_model_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   if Region_mask:
      ts_fil='Region_'+ts_fil
   figure.canvas.print_figure(ts_fil,bbox_inches='tight',dpi=dpi)
   logger.info("Successfull printing:  %s"%ts_fil)

   fnm='TS'
   # plot short
   figure, ax = plt.subplots()
   nplts=1  
   ll=-1*len(tid)
   if filename2:
      ll=-1*len(tid_2)
   ax.plot_date(tid[ll:], dt_cnl[ll:], '-o',color='g',ms=3,  label=Labl1)
   if filename2:
      ax.plot_date(tid_2, dt_2, '-v',color='orange',ms=3,  label=Labl2)
      nplts=nplts+1 
      fnm=fnm+'filenme2'
   if filename5:
      ax.plot_date(tid_5, dt_5,'-+', color='m', label=Labl5)
      nplts=nplts+1  

   ax.xaxis.set_major_locator(years)
   ax.xaxis.set_major_formatter(yearsFmt)
   ax.xaxis.set_minor_locator(months)
   ax.autoscale_view()
   # format the coords message box
   def price(x):
       return '$%1.2f' % x
   ax.fmt_xdata = DateFormatter('%Y-%m-%d')
   ax.fmt_ydata = price
   ax.grid(True)
   figure.autofmt_xdate()
   legend = plt.legend(loc='upper left',fontsize=8)
   plt.title("Area-averaged: %s(%d)"%(fieldname,fieldlevel))
   #plt.xlabel('dayes')
   plt.ylabel("%s(%d)"%(fieldname,fieldlevel))
   #plt.title('Pakistan India Population till 2007')
   ts_fil=fnm+'_'+yyyy1cnt+"_model_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   if Region_mask:
      ts_fil='Region_'+fnm+'_'+yyyy1cnt+"_model_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   figure.canvas.print_figure(ts_fil,bbox_inches='tight',dpi=dpi)
   logger.info("Successfull printing:  %s"%ts_fil)

   # plot short  mean error
   figure, ax = plt.subplots()
   print("diff_dt_cnl[:]=",diff_dt_cnl[:])
   nplts=1  
   ll=-1*len(tid)
   if filename2:
      ll=-1*len(tid)
   ax.plot_date(tid[ll:], diff_dt_cnl[ll:], '-o',color='g',ms=3,  label=Labl1)
   if filename2:
      ax.plot_date(tid_2, diff_dt_2, '-v',color='orange',ms=3,  label=Labl2)
      nplts=nplts+1  
 
   ax.xaxis.set_major_locator(years)
   ax.xaxis.set_major_formatter(yearsFmt)
   ax.xaxis.set_minor_locator(months)
   ax.autoscale_view()
   # format the coords message box
   def price(x):
       return '$%1.2f' % x
   ax.fmt_xdata = DateFormatter('%Y-%m-%d')
   ax.fmt_ydata = price
   ax.grid(True)
   figure.autofmt_xdate()
   legend = plt.legend(loc='upper left',fontsize=8)
   plt.title("Mean diff: %s(%d)"%(fieldname,fieldlevel))
   plt.ylabel("diff:%s(%d)"%(fieldname,fieldlevel))
   ts_fil='Mdiff_'+fnm+"_model_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   if Region_mask:
      ts_fil='Region_Mdiff_'+fnm+"_model_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   figure.canvas.print_figure(ts_fil,bbox_inches='tight',dpi=dpi)
   logger.info("Successfull printing:  %s"%ts_fil)




   # plot rooot mean square RMSE  error
   figure, ax = plt.subplots()
   nplts=1  
   ll=-1*len(tid)
   if filename2:
      ll=-1*len(tid_2)
   ax.plot_date(tid[ll:], rmse_dt_cnl[ll:], '-o',color='g',ms=3,  label=Labl1)
   if filename2:
      ax.plot_date(tid_2, rmse_dt_2, '-v',color='orange',ms=3,  label=Labl2)
      nplts=nplts+1  
  
   ax.xaxis.set_major_locator(years)
   ax.xaxis.set_major_formatter(yearsFmt)
   ax.xaxis.set_minor_locator(months)
   ax.autoscale_view()
   # format the coords message box
   def price(x):
       return '$%1.2f' % x
   ax.fmt_xdata = DateFormatter('%Y-%m-%d')
   ax.fmt_ydata = price
   ax.grid(True)
   figure.autofmt_xdate()
   legend = plt.legend(loc='upper left',fontsize=8)
   plt.title("RMSE: %s(%d)"%(fieldname,fieldlevel))
   #plt.xlabel('dayes')
   plt.ylabel("RMSE:%s(%d)"%(fieldname,fieldlevel))
   #plt.title('Pakistan India Population till 2007')
   ts_fil='RMSE_'+fnm+"_model_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   if Region_mask:
      ts_fil='Region_RMSE_'+fnm+"_model_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   figure.canvas.print_figure(ts_fil,bbox_inches='tight',dpi=dpi)
   logger.info("Successfull printing:  %s"%ts_fil)

if __name__ == "__main__" :
   class ClimParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       setattr(args, self.dest, tmp)
   class WindowParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [int(elem) for elem in tmp[0:4]]
       setattr(args, self.dest, tmp)
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--dpi',        type=int, default=180)
   parser.add_argument('--clim',       action=ClimParseAction,default=None,help="range of colormap")
   parser.add_argument('--masklim',    action=ClimParseAction,default=None,help="mask limits")
   parser.add_argument('--cmap',       type=str,default="jet",help="matplotlib colormap to use")
   parser.add_argument('--filetype',   type=str, default="archive",help="filetype to plot (archive by default)")
   parser.add_argument('--window',     action=WindowParseAction, help='firsti,firstj,lasti,lastj', default=None)
   parser.add_argument('--idm',        type=int, help='Grid dimension 1st index []', default=None)
   parser.add_argument('--jdm',        type=int, help='Grid dimension 2nd index []', default=None)
   parser.add_argument('--datetime1',  action=DateTimeParseAction)
   parser.add_argument('--datetime2',  action=DateTimeParseAction)
   parser.add_argument('--vector',     type=str,default="",help="Denotes second vector component")
   parser.add_argument('--tokml',      action="store_true", default=False)
   parser.add_argument('--filename2',   help="",nargs="+")
   parser.add_argument('--filename5',   help="",nargs="+")
   parser.add_argument('fieldname',  type=str)
   parser.add_argument('fieldlevel', type=int)
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()

   main(args.filename,args.fieldname,args.fieldlevel,
         idm=args.idm,jdm=args.jdm,clim=args.clim,filetype=args.filetype,
         window=args.window,cmap=args.cmap,
         datetime1=args.datetime1,datetime2=args.datetime2,
         vector=args.vector,
         tokml=args.tokml,
         masklim=args.masklim,
         filename2=args.filename2,
         filename5=args.filename5,
         dpi=args.dpi)


