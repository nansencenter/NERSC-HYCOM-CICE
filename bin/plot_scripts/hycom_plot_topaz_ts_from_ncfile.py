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
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size


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

#Example:
#python  ./hycom_plot_topaz_ts_from_ncfile.py  "srfhgt00"  ./Spinup_1/TP6archm_monthly_201{7..9}{01..12}_mean.nc --filename2 ./Spin_up2/TP6archm_monthly_201{7..9}{01..12}_mean.nc
#
def main(myfiles,fieldname,
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
      dpi=180) :


   cmap=matplotlib.pyplot.get_cmap("jet")
   if tokml :
      ab = abf.ABFileGrid("regional.grid","r")
      plon=ab.read_field("plon")
      plat=ab.read_field("plat")
      ab.close()

   
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
      maskd[plat>80]=np.nan
      maskd[plat<50]=np.nan
      maskd[plon>60]=np.nan
      maskd[plon<-50]=np.nan

   Nordic_mask=maskd   

   proj=ccrs.Stereographic(central_latitude=90.0,central_longitude=-40.0)
   pxy = proj.transform_points(ccrs.PlateCarree(), plon, plat)
   px=pxy[:,:,0]
   py=pxy[:,:,1]
   x,y=np.meshgrid(np.arange(plon.shape[1]),np.arange(plon.shape[0]))


   if vector :
      logger.info("Vector component 1:%s"%fieldname)
      logger.info("Vector component 2:%s"%vector) 

   #---------------
   fieldlevel=0
   Err_map=1
   #freezp=-2.5
   freezp=-1.8
   Point_tid=True
   Point_tid=False
   if Point_tid:
      ix=1394
      jy=267
   sum_fld1=maskd
   sum_fld1[~np.isnan(sum_fld1)]=0.0
   Clim_arr=np.zeros((plon.shape[0],plon.shape[1],12))
   #--------------- 
   # compute for TP6 files
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   counter=0
   file_count=0
   sum_fld1=maskd
   sum_fld1[~np.isnan(sum_fld1)]=0.0
       
   #-----------------------------------------
   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   onemm=9.806
   counter=0
   file_count=0
   sum_fld1=maskd
   sum_fld1[~np.isnan(sum_fld1)]=0.0
   dt_cnl=np.zeros(len(myfiles))
   diff_dt_cnl=np.zeros(len(myfiles))
   rmse_dt_cnl=np.zeros(len(myfiles))
   Labl1=myfiles[0][:28]
   #Labl1="CNTL: New prsbas=0"
   yyyy1=myfiles[0][-14:-10]
   print("myfiles[0]=",myfiles[0])
   print("yyy1=", yyyy1)
   base = datetime.datetime(int(yyyy1), 1, 15)
   tid=np.array([base + relativedelta(months=i) for i in range(len(myfiles))])
   counter=0
   file_count=0
   sum_fld1=maskd
   sum_fld1[~np.isnan(sum_fld1)]=0.0
   logger.info(">>>>>--------------------------Processing the first files=  myfiles")
   if "salin" in fieldname:
      fieldname="salin01"
   for ncfile0 in myfiles :
       logger.info("Now processing  %s"%ncfile0)
       fh = Dataset(ncfile0, mode='r')
       fld_arr = fh.variables[fieldname][:]
       if "srfhgt" in fieldname:
          #convert to "m"
          fld_arr= fld_arr/9.806
       print("fld_arr.shpe", fld_arr.shape)
       tot=fld_arr.shape[0]
       fh.close()
       for ii in range(tot):
          fld=fld_arr[ii,:,:]
          print('mn,mx=',fld.min(),fld.max(), 'count=',counter)
          dt_cnl[counter]=np.nanmean(fld)
          if Point_tid:
             dt_cnl[counter]=fld[jy,ix]
          print("fld.shape", fld.shape)
          print("Nordic_mask.shape", Nordic_mask.shape)
          counter=counter+1
          sum_fld1=sum_fld1+fld
          del fld
         # End i_intloop      
   print('Computing the avearge of file_counter= ', file_count, 'counter=',counter)
   #next experminet
   if filename2:
      dt_2=np.zeros(len(filename2))
      diff_dt_2=np.zeros(len(filename2))
      rmse_dt_2=np.zeros(len(filename2))
      yyyy2=filename2[0][-14:-10]
      print("filename2[0]=",filename2[0])
      print("yyy1=", yyyy2)
      tid_2=np.array([datetime.datetime(int(yyyy2), 1, 15)+relativedelta(months=i) for i in range(len(filename2))])
      Labl2=filename2[0][:28]
      counter=0
      file_count=0
      sum_fld1=maskd
      sum_fld1[~np.isnan(sum_fld1)]=0.0
      logger.info(">>>>>--------------------------Processing the first files=  myfiles")
      for ncfil in filename2:
          logger.info("Now processing  %s"%ncfil)
          fh = Dataset(ncfil, mode='r')
          fld_arr = fh.variables[fieldname][:]
          if "srfhgt" in fieldname:
             fld_arr= fld_arr/9.806
          print("fld_arr.shpe", fld_arr.shape)
          tot=fld_arr.shape[0]
          fh.close()
          for ii in range(tot):
             fld=fld_arr[ii,:,:]
             #fld=np.ma.masked_where(fld<freezp,fld)
             print('mn,mx=',fld.min(),fld.max(), 'count=',counter)
             dt_2[counter]=np.nanmean(fld)
             if Point_tid:
                dt_2[counter]=fld[jy,ix]
             counter=counter+1
             sum_fld1=sum_fld1+fld
             del fld

   #---------------------------------------
   figure, ax = plt.subplots()
   years = YearLocator()   # every year
   months = MonthLocator()  # every month
   yearsFmt = DateFormatter('%Y')
   #ax=figure.add_subplot(111)
   nplts=1  
   ax.plot_date(tid, dt_cnl, '-o',color='g',ms=3,  label=Labl1)
   if filename2:
      ax.plot_date(tid_2, dt_2, '-v',color='blue',ms=3,  label=Labl2)
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
   legend = plt.legend(loc='upper right',fontsize=8)
   if Point_tid:
      plt.title("Point:(lon,lat)=("+str(plon[jy,ix])+','+str(plat[jy,ix])+"): %s(%d)"%(fieldname,fieldlevel))
   else:
      plt.title("Area-averaged: %s(%d)"%(fieldname,fieldlevel))
   #plt.xlabel('dayes')
   if "srfhgt" in fieldname:
      plt.ylabel("%s[m]"%(fieldname))
   else:
      plt.ylabel("%s(%d)"%(fieldname,fieldlevel))
   #plt.title('Pakistan India Population till 2007')
   ts_fil="Time_series_cntl%s_%02d_%02d"%(fieldname,fieldlevel,len(myfiles))
   if Region_mask:
      ts_fil='Region_'+ts_fil
   if Point_tid:
      ts_fil='Point_ix'+str(ix)+'jy'+str(jy)+ts_fil
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
   parser.add_argument('fieldname',  type=str)
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()

   main(args.filename,args.fieldname,
         idm=args.idm,jdm=args.jdm,clim=args.clim,filetype=args.filetype,
         window=args.window,cmap=args.cmap,
         datetime1=args.datetime1,datetime2=args.datetime2,
         vector=args.vector,
         tokml=args.tokml,
         masklim=args.masklim,
         filename2=args.filename2,
         dpi=args.dpi)


