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
python ./hycom_plot_TOPAZ_vs_clim.py --clim=-2,2  temp  level TP6monthly_201{6..7}*.a
#
python ./hycom_plot_TOPAZ_vs_clim.py --clim=-2,2  temp  1 TP6monthly_201{6..7}*.a
python ./hycom_plot_TOPAZ_vs_clim.py --clim=-2,2  salin  level TP6monthly_201{6..7}*.a
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
   if 'temp' or 'sal' in fieldname:
      cmap= matplotlib.colors.LinearSegmentedColormap('my_colormap',LinDic)
   else:
      cmap= matplotlib.colors.LinearSegmentedColormap('my_colormap',LinDic)
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
   Region_mask=True
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
   sum_fld1=maskd
   sum_fld1[~np.isnan(sum_fld1)]=0.0
   Clim_arr=np.zeros((plon.shape[0],plon.shape[1],12))
   if 'tem' or 'sal' in fieldname:
      counter=0
      if 'tem' in fieldname:
         rlxfile0="/cluster/work/users/achoth/TP5a0.06/relax/050/relax_tem.a"
      if 'sal' in fieldname:
         rlxfile0="/cluster/work/users/achoth/TP5a0.06/relax/050/relax_sal.a"
      rlx_afile = abf.AFile(ab.idm,ab.jdm,rlxfile0,"r")
      lyr=fieldlevel
      record_num=1
      record_var=record_num-1 
      fld = rlx_afile.read_record(record_var)
      print('mn,mx  data='),fld.min(),fld.max()
      kdm=50
      dt_clim=np.zeros(12)
      for mnth in range(12) :
          fld1=rlx_afile.read_record(mnth*kdm+lyr-1)
          print('record, mn,mx  data='), kdm*mnth, fld1.min(),fld1.max()
          # Intloop used to read more fields in one file. Only for forcing for now
          dt_clim[mnth]=mod_hyc2plot.spatiomean(fld1,maskd)
          sum_fld1=sum_fld1+fld1
          Clim_arr[:,:,mnth]=fld1[:,:]
          counter=counter+1
          print('counter='),counter
          del  fld1
      Clim_Avg=sum_fld1/counter
      del  sum_fld1

   #---------------filename 
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
   Labl1="Model: "+fieldname
   if "SPRBAS_0" in myfiles[0]:
      Labl1="CNTL: prsbas=0 "
   if filename2:
      dt_2=np.zeros(len(filename2))
      diff_dt_2=np.zeros(len(filename2))
      rmse_dt_2=np.zeros(len(filename2))
      yyyy1=filename2[0][-9:-5]
      print("filename2[0]="),filename2[0][-9:-5]
      print("filename2[0]="),filename2[0]
      print("yyy1="), yyyy1
      tid_2=np.array([datetime.datetime(int(yyyy1), 1, 15) \
         + relativedelta(months=i) for i in range(len(filename2))])
      Labl2="filename2"
      Labl2="Corrected"
      if "erai" in filename2[0]:
         Labl2="CNTL: prsbas=1e5 "
   
   yyyy1cnt=myfiles[0][-8:-5]
   print("myfiles[0]="),myfiles[0][-9:-5]
   print("myfiles[0]="),myfiles[0]
   print("yyy1cnt="), print(yyyy1cnt)
   base = datetime.datetime(int(yyyy1cnt), 1, 15)
   tid=np.array([base + relativedelta(months=i) for i in range(len(myfiles))])
   if len(myfiles)==36:
      base = datetime.datetime(int(yyyy1cnt), 1, 15)
      tid=np.array([base + relativedelta(months=i) for i in range(len(myfiles))])

   nmexp=1
   if filename2:
      nmexp=nmexp+1
   print('processing data from No runs ==##############>>>>>>>>>>>>>>>>>>>>>>>'), nmexp
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
            print('---------mn,mx  data='),fld.min(),fld.max()
            sum_fld1=sum_fld1+fld
            cindx=np.remainder(counter,12)
            print("counter"), counter, print("cindx="), cindx
            if iii==0:
               dt_cnl[counter]=mod_hyc2plot.spatiomean(fld,Nordic_mask)
               diff_dt_cnl[counter]=mod_hyc2plot.spatiomean(fld[:,:]-Clim_arr[:,:,cindx],Nordic_mask)
               rmse_dt_cnl[counter]=np.sqrt(mod_hyc2plot.spatiomean((fld[:,:]-Clim_arr[:,:,cindx])**2,Nordic_mask))
               Labl=Labl1
            if iii==1 and filename2:
               dt_2[counter]=mod_hyc2plot.spatiomean(fld,Nordic_mask)
               diff_dt_2[counter]=mod_hyc2plot.spatiomean(fld[:,:]-Clim_arr[:,:,cindx],Nordic_mask)
               rmse_dt_2[counter]=np.sqrt(mod_hyc2plot.spatiomean((fld[:,:]-Clim_arr[:,:,cindx])**2,Nordic_mask))
               Labl=Labl2
            # Apply mask if requested
            counter=counter+1
            file_count=file_count+1
            del fld
         # End i_intloop
      print('Computing the avearge of file_counter= '), print(file_count), print('counter='),print(counter), print('cindx='),print(cindx)
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

      # Print figure.
      aspect = 40
      pad_fraction = 0.25
      divider = make_axes_locatable(ax)
      width = axes_size.AxesY(ax, aspect=1./aspect)
      pad = axes_size.Fraction(pad_fraction, width)
      cax = divider.append_axes("right", size=width, pad=pad)
      cb=ax.figure.colorbar(P,cax=cax,extend='both')
      if clim is not None : P.set_clim(clim)
      ax.set_title("Diff:%s(%d)"%(fieldname,fieldlevel)+' :( Model - Clim )')
      # Print figure.
      fnamepng_template=myfiles[0][-20:-5].replace("/",'') +"_Avg_TP6_%s_%d_%03d_iii%03d_Avg.png"
      if Region_mask:
         fnamepng_template='Region_'+yyyy1cnt+myfiles[0][1:11].replace("/",'') \
               +"Avg_TP5_%s_%d_%03d_iii%03d_Avg.png"
      fnamepng=fnamepng_template%(fieldname,fieldlevel,counter,iii)
      logger.info("output in  %s"%fnamepng)
      figure.canvas.print_figure(fnamepng,bbox_inches='tight',dpi=dpi)
      ax.clear()
      cb.remove()
      datmen=np.nanmean(fld_diff)
      spatiodatmen=mod_hyc2plot.spatiomean(fld_diff,Nordic_mask)
      print('-----------mean diff data, spatio='), datmen,spatiodatmen
      del sum_fld1
      #---------------------------------------


   print('tid len='), print(tid.shape)
   if filename2:
      print('dt_2='), print(dt_2.shape)
   tid_clim=np.array([base + relativedelta(months=i) for i in range(12)])
   figure, ax = plt.subplots()
   rpt=len(dt_cnl)/12
   dt_clim_cat=dt_clim
   for ii in range(int(rpt-1)):
      print("concatenate ")
      dt_clim_cat=np.concatenate([dt_clim_cat,dt_clim])
      
     
   years = YearLocator()   # every year
   months = MonthLocator()  # every month
   yearsFmt = DateFormatter('%Y')
   nplts=1  
   ax.plot_date(tid, dt_cnl, '-o',color='g',ms=3,  label=Labl1)
   if filename2:
      ax.plot_date(tid_2, dt_2, '-v',color='orange',ms=3,  label=Labl2)
      nplts=nplts+1  
   if filename5:
      ax.plot_date(tid_5, dt_5,'--', color='m', label=Labl5)
      nplts=nplts+1  
   if 'tem' or 'sal' in fieldname:
      ax.plot_date(tid[0:len(dt_cnl)], dt_clim_cat[:],':' ,color='black', label='Clim:'+fieldname)
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
   plt.ylabel("%s(%d)"%(fieldname,fieldlevel))
   ts_fil="time_series_cntl_flx_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   if Region_mask:
      ts_fil='Region_'+ts_fil
   figure.canvas.print_figure(ts_fil,bbox_inches='tight',dpi=dpi)
   logger.info("Successfull printing:  %s"%ts_fil)

   #-----------------
   # plot short  mean error
   figure, ax = plt.subplots()
   print("diff_dt_cnl[:]="),diff_dt_cnl[:]
   nplts=1  
   ll=-1*len(tid)
   if filename2:
      ll=-1*len(tid)
   ax.plot_date(tid[ll:], diff_dt_cnl[ll:], '-o',color='g',ms=3)
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
   plt.title("Mean diff:Model-Clim: %s(%d)"%(fieldname,fieldlevel))
   plt.ylabel("diff:%s(%d)"%(fieldname,fieldlevel))
   ts_fil='Mdiff'+"ST_cntl_flx_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   if Region_mask:
      ts_fil='Region_Mdiff'+"ST_cntl_flx_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   figure.canvas.print_figure(ts_fil,bbox_inches='tight',dpi=dpi)
   logger.info("Successfull printing:  %s"%ts_fil)




   # plot rooot mean square RMSE  error
   figure, ax = plt.subplots()
   nplts=1  
   ll=-1*len(tid)
   if filename2:
      ll=-1*len(tid_2)
   ax.plot_date(tid[ll:], rmse_dt_cnl[ll:], '-o',color='g',ms=3)
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
   plt.title("RMSE: (Model-Clim) %s(%d)"%(fieldname,fieldlevel))
   plt.ylabel("RMSE:%s(%d)"%(fieldname,fieldlevel))
   ts_fil='RMSE'+"ST_cntl_flx_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   if Region_mask:
      ts_fil='Region_RMSE'+"ST2007_cntl_flx_%s_%02d_%02d"%(fieldname,fieldlevel,nplts)
   figure.canvas.print_figure(ts_fil,bbox_inches='tight',dpi=dpi)
   logger.info("Successfull printing:  %s"%ts_fil)
   logger.info("End --------------------------------------------------------------- printing:  %s"%ts_fil)

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


