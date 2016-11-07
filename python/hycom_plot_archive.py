#!/usr/bin/env python
import modeltools.hycom
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
from mpl_toolkits.basemap import Basemap
import logging
import datetime
import re

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


def main(myfiles,fieldname,fieldlevel,idm=None,jdm=None,clim=None,filetype="archive",window=None,
      cmap="jet",datetime1=None,datetime2=None,vector="") :


   #cmap=matplotlib.pyplot.get_cmap("jet")
   cmap=matplotlib.pyplot.get_cmap(cmap)

   # Prelim support for projections
   bm=None
   #ab = abfile.ABFileGrid("regional.grid","r")
   #plon=ab.read_field("plon")
   #plat=ab.read_field("plat")
   #bm = Basemap(width=9000000,height=9000000,
   #      resolution='i',projection='stere',\
   #      lat_ts=70,lat_0=90,lon_0=-40.)
   #x,y=bm(plon,plat)

   if vector :
      logger.info("Vector component 1:%s"%fieldname)
      logger.info("Vector component 2:%s"%vector) 

   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)

   counter=0
   for myfile0 in myfiles :


      logger.info("Now processing  %s"%myfile0)
      m=re.match("(.*)\.[ab]",myfile0)
      if m :
         myfile=m.group(1)
      else :
         myfile=myfile0

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
            raise NotImplementedError,"Filetype %s not implemented"%filetype

         if not window :
            J,I=numpy.meshgrid(numpy.arange(fld1.shape[0]),numpy.arange(fld1.shape[1]))
         else :
            J,I=numpy.meshgrid( numpy.arange(window[1],window[3]),numpy.arange(window[0],window[2]))

         # Create simple figure
         #figure = matplotlib.pyplot.figure(figsize=(8,8))
         #ax=figure.add_subplot(111)
         if vector : 
            fld = numpy.sqrt(fld1**2+fld2**2)
         else :
            fld=fld1
         #print fld.min(),fld.max()

         if bm is not None :
            P=bm.pcolormesh(x[J,I],y[J,I],fld[J,I],cmap=cmap)
            bm.drawcoastlines()
            bm.fillcontinents(color='.5',lake_color='aqua')
            # draw parallels and meridians.
            bm.drawparallels(numpy.arange(-80.,81.,20.))
            bm.drawmeridians(numpy.arange(-180.,181.,20.))
            bm.drawmapboundary() #fill_color='aqua')
         else :
            P=ax.pcolormesh(I,J,fld[J,I],cmap=cmap)
         cb=ax.figure.colorbar(P)
         if clim is not None : P.set_clim(clim)
         ax.set_title("%s:%s(%d)"%(myfile0,fieldname,fieldlevel))

         if vector: 
            skip=10
            I2=I[::skip,::skip]
            J2=J[::skip,::skip]
            ax.quiver(I2,J2,fld1[J2,I2],fld2[J2,I2])


         # Print figure.
         figure.canvas.print_figure("tst%03d.png"%counter,dpi=180)

         ax.clear()
         cb.remove()
         counter=counter+1

      # End i_intloop



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
   parser.add_argument('--clim',       action=ClimParseAction,default=None)
   parser.add_argument('--cmap',       type=str,default="jet")
   parser.add_argument('--filetype',   type=str, help='',default="archive")
   parser.add_argument('--window',     action=WindowParseAction, help='firsti,firstj,lasti,lastj', default=None)
   parser.add_argument('--idm',        type=int, help='Grid dimension 1st index []', default=None)
   parser.add_argument('--jdm',        type=int, help='Grid dimension 2nd index []', default=None)
   parser.add_argument('--datetime1',  action=DateTimeParseAction)
   parser.add_argument('--datetime2',  action=DateTimeParseAction)
   parser.add_argument('--vector',     type=str,default="")
   parser.add_argument('fieldname',  type=str)
   parser.add_argument('fieldlevel', type=int)
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()

   main(args.filename,args.fieldname,args.fieldlevel,
         idm=args.idm,jdm=args.jdm,clim=args.clim,filetype=args.filetype,
         window=args.window,cmap=args.cmap,
         datetime1=args.datetime1,datetime2=args.datetime2,
         vector=args.vector)


