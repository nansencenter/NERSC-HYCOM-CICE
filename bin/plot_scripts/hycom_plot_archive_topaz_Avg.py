#!/usr/bin/env python
import modeltools.hycom
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import abfile.abfile as abf
import numpy as np
import logging
import datetime
import re
import scipy.interpolate
import os.path
import cmocean
import mod_hyc2plot
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

"""
Script for computing and ploting the time-averaged of 2D field from agiven ab files:
  python ./hycom_plot_archive_topaz_Avg.py  --clim=26,36 salin  1 --exptid="Spinup1_" TP6archm.2020_00?_12.a
  python ./hycom_plot_archive_topaz_Avg.py --clim=-1.8,13 temp  1 --exptid="Spinup1_" TP6archm.2020_00?_12.a
   
  #plot speed and velcity vectors:
  python ./hycom_plot_archive_topaz_Avg.py --clim=0,0.3 --vector=v-vel.  u-vel.  1 --exptid="Spinup1" TP6archm.2020_0??_12.a

"""


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
      exptid="",
      masklim=None,
      dpi=180) :

   cmap=plt.get_cmap("jet")
   LinDic=mod_hyc2plot.cmap_dict('sawtooth_fc100.txt')
   my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',LinDic)
   cmap=my_cmap

   if vector:
      cmap=cmocean.cm.speed

   if tokml :
      ab = abf.ABFileGrid("regional.grid","r")
      plon=ab.read_field("plon")
      plat=ab.read_field("plat")
      ab.close()

   ab = abf.ABFileGrid("regional.grid","r")
   plon=ab.read_field("plon")
   plat=ab.read_field("plat")
   ab.close()

   
   proj=ccrs.Stereographic(central_latitude=90.0,central_longitude=-40.0)
   pxy = proj.transform_points(ccrs.PlateCarree(), plon, plat)
   px=pxy[:,:,0]
   py=pxy[:,:,1]
   x,y=np.meshgrid(np.arange(plon.shape[1]),np.arange(plon.shape[0]))

   

   if vector :
      logger.info("Vector component 1:%s"%fieldname)
      logger.info("Vector component 2:%s"%vector) 

   figure = plt.figure(figsize=(8,8))
   ax=figure.add_subplot(111)
   ax.set_facecolor('xkcd:gray')
   onemm=9.806
   counter=0
   sum_fld=np.zeros(plon.shape)
   sumf1=np.zeros(plon.shape)
   sumf2=np.zeros(plon.shape)
   file_count=0
   for myfile0 in myfiles :

      # Open files, and return some useful stuff.
      # ab2 i used in case of vector
      # rdtimes is used for plotting forcing fields
      n_intloop,ab,ab2,rdtimes = open_file(myfile0,filetype,fieldname,fieldlevel,datetime1=datetime1,datetime2=datetime2,vector=vector,idm=idm,jdm=jdm)
      # Intloop used to read more fields in one file. Only for forcing for now
      for i_intloop in range(n_intloop) :

         # Read ab file of different types
         if filetype == "archive" :
            fld1 = ab.read_field(fieldname,fieldlevel)
            if vector :fld2=ab.read_field(vector,fieldlevel)
         elif filetype == "regional.depth" :
            fld1 = ab.read_field(fieldname)
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

         print('mnfld1,mxfld1=',fld1.min(),fld1.max()) 
         # Create scalar field for vectors
         if vector : 
            fld = np.sqrt(fld1**2+fld2**2)
            print('mnfld2,mxfld2=',fld2.min(),fld2.max()) 
         else :
            fld=fld1

         # Apply mask if requested
         if masklim :
            fld = np.ma.masked_where(fld<=masklim[0],fld)
            fld = np.ma.masked_where(fld>=masklim[1],fld)

         sum_fld=sum_fld+fld
         counter=counter+1
         if vector :
            sumf1=sumf1+fld1
            sumf2=sumf2+fld2
      file_count=file_count+1
      # End i_intloop
   print('Computing the avearge of file_counter= ', file_count, 'counter=',counter)
   if file_count> 0:
      fld_Avg=sum_fld/file_count
      print('mn_Avg_fld,mx_Avg_fld=',fld_Avg.min(),fld_Avg.max()) 
   if vector:
      f1_avg= sumf1/file_count
      f2_avg= sumf2/file_count
      print('mn_avg_fld1,mx_avg_fld1=',f1_avg.min(),f1_avg.max())  
      print('mn_avg_fld2,mx_avg_fld2=',f2_avg.min(),f2_avg.max())  

   if fieldname=='k.e.' :
         P=plt.pcolormesh(x[J,I],y[J,I],np.log10(fld_Avg[J,I]),cmap=cmap, shading='auto')
   elif fieldname=='srfhgt' :
         P=plt.pcolormesh(x[J,I],y[J,I],(fld_Avg[J,I]/onemm),cmap=cmap, shading='auto')
   else :
         P=plt.pcolormesh(x[J,I],y[J,I],fld_Avg[J,I],cmap=cmap, shading='auto')
   
   if 'temp' in fieldname:
            P1=plt.contour(x[J,I],y[J,I],fld_Avg[J,I],levels=[-1.,1,4.0,8], \
                       colors=('w',),linestyles=('-',),linewidths=(1.5,))
            plt.clabel(P1, fmt = '%2.1d', colors = 'w', fontsize=10) #contour line labels
      

   if vector: 
      skip=10
      logger.info("ploting quiver .......>>> %s"%vector)
      I2=I[::skip,::skip]
      J2=J[::skip,::skip]
      plt.quiver(x[J2,I2],y[J2,I2],f1_avg[J2,I2],f2_avg[J2,I2])
                  
   ##
   # Print figure and remove wite space.
   ax.set_facecolor('xkcd:gray')
   aspect = 40
   pad_fraction = 0.25
   divider = make_axes_locatable(ax)
   width = axes_size.AxesY(ax, aspect=1./aspect)
   pad = axes_size.Fraction(pad_fraction, width)
   cax = divider.append_axes("right", size=width, pad=pad)
   if vector: 
      cb=ax.figure.colorbar(P,cax=cax,extend='max')
   else:
      cb=ax.figure.colorbar(P,cax=cax,extend='both')
   if clim is not None : P.set_clim(clim)
   ax.set_title('Avg'+"%s:%s(%d)"%(myfile0,fieldname,fieldlevel))
   # Print figure.
   fnamepng_template=exptid+"Avg_%s_%d_%03d_Avg.png"
   fnamepng=fnamepng_template%(fieldname,fieldlevel,counter)
   logger.info("output in  %s"%fnamepng)
   figure.canvas.print_figure(fnamepng,bbox_inches='tight',dpi=dpi)
   ax.clear()
   cb.remove()



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
   parser.add_argument('--exptid'      , type=str,default="")
   parser.add_argument('fieldname',  type=str)
   parser.add_argument('fieldlevel', type=int)
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()

   main(args.filename,args.fieldname,args.fieldlevel,
         idm=args.idm,jdm=args.jdm,clim=args.clim,filetype=args.filetype,
         window=args.window,cmap=args.cmap,
         datetime1=args.datetime1,datetime2=args.datetime2,
         vector=args.vector,
         tokml=args.tokml,exptid=args.exptid,
         masklim=args.masklim,
         dpi=args.dpi)


