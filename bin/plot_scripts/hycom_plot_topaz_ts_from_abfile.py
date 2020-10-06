#/busr/bin/env python
from netCDF4 import Dataset
import modeltools.hycom
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
import logging
import datetime
import re
import scipy.interpolate
import os.path
import matplotlib.pyplot as plt
#import mod_HYCOM_utils as mhu
#import mod_reading as mr
#import myMOD
import matplotlib.dates as mdates
from dateutil.relativedelta import relativedelta
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import cmocean

'''
Plot timeseries of the domain-wide average of a given field: examples 
python  ./hycom_plot_topaz_ts_from_abfile.py  temp  1 ./Spin_up1/TP6monthly_2019_{01..12}.a  
python  ./hycom_plot_topaz_ts_from_abfile.py  temp  1 ./Spin_up1/TP6monthly_2019_{01..12}.a  --filename2 ./NEMO_data2/NM8monthly_2019_{01..12}.a
python  ./hycom_plot_topaz_ts_from_abfile.py  temp  1 ./Spin_up1/TP6monthly_2019_{01..12}.a  --filename2 ./NEMO_data2/NM8monthly_2019_{01..12}.a
The script will also the corespondat fields from climatolgy file.
'''

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
      ab = abfile.ABFileArchv(myfile,"r")
      n_intloop=1
   #elif filetype == "restart" :
   #   tmp = abfile.ABFileRestart(myfile,"r",idm=gfile.idm,jdm=gfile.jdm)
   elif filetype == "regional.depth" :
      ab = abfile.ABFileBathy(myfile,"r",idm=idm,jdm=jdm)
      n_intloop=1
   elif filetype == "forcing" :
      ab = abfile.ABFileForcing(myfile,"r",idm=idm,jdm=jdm)
      if vector :
         file2=myfile.replace(fieldname,vector)
         logger.info("Opening file %s for vector component nr 2"%file2)
         ab2=abfile.ABFileForcing(file2,"r",idm=idm,jdm=jdm)
      if datetime1 is None or datetime2 is None :
         raise NameError,"datetime1 and datetime2 must be specified when plotting forcing files"
      else :
         iday1,ihour1,isec1 = modeltools.hycom.datetime_to_ordinal(datetime1,3)
         rdtime1 = modeltools.hycom.dayfor(datetime1.year,iday1,ihour1,3)
         #
         iday2,ihour2,isec2 = modeltools.hycom.datetime_to_ordinal(datetime2,3)
         rdtime2 = modeltools.hycom.dayfor(datetime2.year,iday2,ihour2,3)
         rdtimes=sorted([elem for elem in ab.field_times if elem >rdtime1 and elem < rdtime2])
         n_intloop=len(rdtimes)
   else :
      raise NotImplementedError,"Filetype %s not implemented"%filetype
   # Check that fieldname is actually in file
   if fieldname not in  ab.fieldnames :
      logger.error("Unknown field %s at level %d"%(fieldname,fieldlevel))
      logger.error("Available fields : %s"%(" ".join(ab.fieldnames)))
      raise ValueError,"Unknown field %s at level %d"%(fieldname,fieldlevel)

   return n_intloop,ab,ab2,rdtimes


#python  ./hycom_plot_topaz_ts_from_abfile.py  temp  1 ./DAILY_archMEAN/Spin_up1/TP6monthly_2019_{01..12}.a  
#python  ./hycom_plot_topaz_ts_from_abfile.py  temp  1 ./DAILY_archMEAN/Spin_up1/TP6monthly_2019_{01..12}.a  --filename2 ./NEMO_data2/NM8monthly_2019_{01..12}.a
#python  ./hycom_plot_topaz_ts_from_abfile.py  temp  1 ./DAILY_archMEAN/Spin_up1/TP6monthly_2019_{01..12}.a  --filename2 ./NEMO_data2/NM8monthly_2019_{01..12}.a


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
      dpi=180) :


   #cmap=matplotlib.pyplot.get_cmap("jet")
   #cmap=matplotlib.pyplot.get_cmap(cmap)
   cmap=matplotlib.pyplot.get_cmap("jet")

   # Prelim support for projections. import basemap only if needed since its usually not needed
   # aaaand installation can sometimes be a bit painful .... bmn is None now, define it if needed
   #bm=None
   from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
   from mpl_toolkits.basemap import Basemap
   
   ab = abfile.ABFileGrid("regional.grid","r")
   plon=ab.read_field("plon")
   plat=ab.read_field("plat")
   scpx=ab.read_field("scpx")
   scpy=ab.read_field("scpy")
   target_lonlats=[plon,plat]
   abdpth = abfile.ABFileBathy('regional.depth',"r",idm=ab.idm,jdm=ab.jdm)
   mdpth=abdpth.read_field('depth')
   maskd=mdpth.data
   maskd[maskd>1e29]=numpy.nan
   #Region_mask=True
   Region_mask=False
   if Region_mask:
      maskd[plat>80]=numpy.nan
      maskd[plat<50]=numpy.nan
      maskd[plon>60]=numpy.nan
      maskd[plon<-50]=numpy.nan
      #bm = Basemap(width=9000000,height=9000000,
   Nordic_mask=maskd   

   bm = Basemap(width=7400000,height=7400000, \
         resolution='i',projection='stere',\
         lat_ts=70,lat_0=85,lon_0=-40.)
   x,y=bm(plon,plat)

   if vector :
      logger.info("Vector component 1:%s"%fieldname)
      logger.info("Vector component 2:%s"%vector) 

   #---------------first read and compute clim  
   Err_map=1
   #freezp=-2.5
   freezp=-1.8
   sum_fld1=maskd
   sum_fld1[~numpy.isnan(sum_fld1)]=0.0
   Clim_arr=numpy.zeros((plon.shape[0],plon.shape[1],12))
   #--------------- 
   # compute for TP6 files
   #-----------------------------------------
   #---------------------------------------------------------------------------------
   #---------------------------------------------------------------------------------
   # filename2
   onemm=9.806
   counter=0
   file_count=0
   sum_fld1=maskd
   sum_fld1[~numpy.isnan(sum_fld1)]=0.0
   dt_cnl=numpy.zeros(len(myfiles))
   diff_dt_cnl=numpy.zeros(len(myfiles))
   rmse_dt_cnl=numpy.zeros(len(myfiles))
   #Labl1=myfiles[0][1:12]
   Labl1="CNTL SST"
   #Labl1="CNTL: prsbas=0"
   Labl1=myfiles[0][:28]
   yyyy1=myfiles[0][-9:-5]
   if "archv." in  myfiles[0]:
      yyyy1=myfiles[0][-13:-9]
   print "myfiles[0]=",myfiles[0]
   print "yyy1=", yyyy1
   if filename2:
      dt_2=numpy.zeros(len(filename2))
      diff_dt_2=numpy.zeros(len(filename2))
      rmse_dt_2=numpy.zeros(len(filename2))
      tid_2=numpy.array([datetime.datetime(int(yyyy1), 1, 15) \
         + relativedelta(months=i) for i in xrange(len(filename2))])
      Labl2=filename2[0][:28]
      counter=0
      file_count=0
      sum_fld1=maskd
      sum_fld1[~numpy.isnan(sum_fld1)]=0.0
      if "srfhgt" in fieldname:
         fieldname="srfhgt"
      elif "temp" in fieldname:
         fieldname="temp"
      for fil0 in filename2 :
           logger.info("Now processing  %s"%fil0)
           n_intloop,ab,ab2,rdtimes = open_file(fil0,filetype,fieldname,fieldlevel,\
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
                 raise NotImplementedError,"Filetype %s not implemented"%filetype
              # Create scalar field for vectors
              print '---------mn,mx  data=',fld1.min(),fld1.max()
              #if "srfhgt" in fieldname:
              #   fld1= fld1/9.806
              print "fld1.shpe", fld1.shape
              print 'mn,mx=',fld1.min(),fld1.max(), 'count=',counter
              dt_2[counter]=numpy.nanmean(fld1)
              counter=counter+1
              sum_fld1=sum_fld1+fld1
              del fld1
   #---------------------------------------------------------------------------------
   #---------------------------------------------------------------------------------
   # filename main
   #tid=numpy.zeros(len(myfiles))
   #tid[:]=range(len(myfiles))
   #base = datetime.datetime(2007, 1, 15)
   base = datetime.datetime(int(yyyy1), 1, 15)
   tid=numpy.array([base + relativedelta(months=i) for i in xrange(len(myfiles))])
   if "archv." in  myfiles[0]:
      tid=numpy.array([base + relativedelta(days=i) for i in xrange(len(myfiles))])
   nmexp=1
   if filename2:
      nmexp=nmexp+1
   print 'processing data from No runs ==##############>>>>>>>>>>>>>>>>>>>>>>>', nmexp
   whole_domain=True
   whole_domain=False
   #
   counter=0
   file_count=0
   sum_fld1=maskd
   sum_fld1[~numpy.isnan(sum_fld1)]=0.0
   logger.info(">>>>>--------------------------Processing the first files=  myfiles")
   for myfile0 in myfiles :
       logger.info("Now processing  %s"%myfile0)
       n_intloop,ab,ab2,rdtimes = open_file(myfile0,filetype,fieldname,fieldlevel,\
            datetime1=datetime1,datetime2=datetime2,vector=vector,idm=idm,jdm=jdm)
       # Intloop used to read more fields in one file. Only for forcing for now
       for i_intloop in range(n_intloop) :
          # Read ab file of different types
          if filetype == "archive" :
             fld1 = ab.read_field(fieldname,fieldlevel)
             if ('temp' in fieldname) and whole_domain:
                vert_fld_sum=0
                for lvl in range(50):
                   print 'lvl=',lvl, fieldlevel
                   fld_lvl = ab.read_field(fieldname,lvl+1)
                   vert_fld_sum=vert_fld_sum + numpy.nanmean(fld_lvl)
                vert_fld_avg=vert_fld_sum/50.0

          elif filetype == "forcing" :
             fld1 = ab.read_field(fieldname,rdtimes[i_intloop])
             if vector :fld2=ab2.read_field(vector,rdtimes[i_intloop])
             logger.info("Processing time %.2f"%rdtimes[i_intloop])
          else :
             raise NotImplementedError,"Filetype %s not implemented"%filetype
          # Create scalar field for vectors
          print '---------mn,mx  data=',fld1.min(),fld1.max()
          #if "srfhgt" in fieldname:
          #   fld1= fld1/9.806
          print "fld1.shpe", fld1.shape
          print 'mn,mx=',fld1.min(),fld1.max(), 'count=',counter
          if ('temp' in fieldname) and whole_domain:
             dt_cnl[counter]=vert_fld_avg
          else:
             dt_cnl[counter]=numpy.nanmean(fld1)
          counter=counter+1
          sum_fld1=sum_fld1+fld1
          del fld1
          # End i_intloop
       print 'Computing the avearge of file_counter= ', file_count, 'counter=',counter
   #---------------------------------------
   #---------------------------------------
   #plot_climatology
   Clim_arr=numpy.zeros((plon.shape[0],plon.shape[1],12))
   if 'tem' in fieldname:
      counter=0
      rlxfile0="/home/sm_alfal/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/relax/070/relax_tem.a"
      rlx_afile = abfile.AFile(ab.idm,ab.jdm,rlxfile0,"r")
      lyr=fieldlevel
      record_num=lyr
      record_var=record_num-1 
      fld = rlx_afile.read_record(record_var)
      print 'mn,mx  data=',fld.min(),fld.max()
      kdm=50
      dt_clim=numpy.zeros(12)
      for mnth in range(12) :
          fld1=rlx_afile.read_record(mnth*kdm+lyr-1)
          logger.debug("File %s, record_var/mnth*kdm %03d/%03d"%(rlxfile0,record_var,mnth*kdm))
          print 'record, mn,mx  data=', kdm*mnth, fld1.min(),fld1.max()
          #fld1=numpy.ma.masked_where(fld1<freezp,fld1)
          #fld1=numpy.where(fld1<freezp,freezp,fld1)
          print 'record, mn,mx  data=', kdm*mnth, fld1.min(),fld1.max()
          # Intloop used to read more fields in one file. Only for forcing for now
          dt_clim[mnth]=numpy.nanmean(fld1)
          #Clim_arr[:,:,mnth]=fld1[:,:]
          counter=counter+1
          print 'counter=',counter
          del  fld1
      #
      tid_clim=numpy.array([base + relativedelta(months=i) for i in xrange(12)])
      #figure, ax = matplotlib.pyplot.figure()
      rpt=len(dt_cnl)/12
      dt_clim_cat=dt_clim
      for ii in range(rpt-1):
         print "concatenate "
         dt_clim_cat=numpy.concatenate([dt_clim_cat,dt_clim])
   #
   #---------------------------------------
   #---------------------------------------
   #if "srfhgt" in fieldname:
   #   srfhgt_mean=numpy.nanmean(dt_cnl)
   #   logger.info("Srfhgt mean ----->>= %.2f"%srfhgt_mean)
   #   dt_cnl[:]=dt_cnl[:]-srfhgt_mean

   #tid_clim=tid[::31]+14
   #figure, ax = matplotlib.pyplot.figure()
   figure, ax = plt.subplots()
   years = YearLocator()   # every year
   months = MonthLocator()  # every month
   yearsFmt = DateFormatter('%Y')
   #ax=figure.add_subplot(111)
   nplts=1  
   ax.plot_date(tid, dt_cnl, '-o',color='g',ms=3,  label=Labl1)
   if 'tem' in fieldname:
      ax.plot_date(tid[0:len(dt_cnl)], dt_clim_cat[:],':' ,color='black', label='Phc-Clim.')
   if filename2:
      ax.plot_date(tid_2, dt_2, '-v',color='orange',ms=3,  label=Labl2)
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
   if "k.e" in fieldname:
      fieldname="KE"
   if "u-vel" in fieldname:
      fieldname="u-vel"
   ts_fil="Time_series_cntl%s_%02d_%02d"%(fieldname,fieldlevel,counter)
   if Region_mask:
      ts_fil='Region_'+ts_fil
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
         dpi=args.dpi)


