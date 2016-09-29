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


def main(myfiles,fieldname,fieldlevel,idm=None,jdm=None,clim=None,filetype="archive",window=None,cmap="jet") :


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

   figure = matplotlib.pyplot.figure(figsize=(8,8))
   ax=figure.add_subplot(111)

   for i,myfile0 in enumerate(myfiles) :


      logger.info("Now processing  %s"%myfile0)
      m=re.match("(.*)\.[ab]",myfile0)
      if m :
         myfile=m.group(1)
      else :
         myfile=myfile0

      if filetype == "archive" :
         ab = abfile.ABFileArchv(myfile,"r")
      #elif filetype == "restart" :
      #   tmp = abfile.ABFileRestart(myfile,"r",idm=gfile.idm,jdm=gfile.jdm)
      elif filetype == "regional.depth" :
         ab = abfile.ABFileBathy(myfile,"r",idm=idm,jdm=jdm)
      else :
         raise NotImplementedError,"Filetype %s not implemented"%filetype

      # Check that fieldname is actually in file
      if fieldname not in  ab.fieldnames :
         logger.error("Unknown field %s at level %d"%(fieldname,fieldlevel))
         logger.error("Available fields : %s"%(" ".join(ab.fieldnames)))
         raise ValueError,"Unknown field %s at level %d"%(fieldname,fieldlevel)

      # Read ab file of different types
      if filetype == "archive" :
         fld = ab.read_field(fieldname,fieldlevel)
      elif filetype == "regional.depth" :
         fld = ab.read_field(fieldname)
      else :
         raise NotImplementedError,"Filetype %s not implemented"%filetype

      if not window :
         J,I=numpy.meshgrid(numpy.arange(fld.shape[0]),numpy.arange(fld.shape[1]))
      else :
         J,I=numpy.meshgrid( numpy.arange(window[1],window[3]),numpy.arange(window[0],window[2]))



      # Create simple figure
      #figure = matplotlib.pyplot.figure(figsize=(8,8))
      #ax=figure.add_subplot(111)
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

      # Print figure.
      figure.canvas.print_figure("tst%03d.png"%i,dpi=180)

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
       print tmp
       tmp = [int(elem) for elem in tmp[0:4]]
       print tmp
       setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--clim',       action=ClimParseAction,default=None)
   parser.add_argument('--cmap',       type=str,default="jet")
   parser.add_argument('--filetype',   type=str, help='',default="archive")
   parser.add_argument('--window',     action=WindowParseAction, help='firsti,firstj,lasti,lastj', default=None)
   parser.add_argument('--idm',        type=int, help='Grid dimension 1st index []', default=None)
   parser.add_argument('--jdm',        type=int, help='Grid dimension 2nd index []', default=None)
   parser.add_argument('fieldname',  type=str)
   parser.add_argument('fieldlevel', type=int)
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()

   main(args.filename,args.fieldname,args.fieldlevel,idm=args.idm,jdm=args.jdm,clim=args.clim,filetype=args.filetype,window=args.window,cmap=args.cmap)


