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


def main(myfiles,fieldname,fieldlevel,idm=None,jdm=None,clim=None,filetype="archive") :




   cmap=matplotlib.pyplot.get_cmap("jet")

   for i,myfile0 in enumerate(myfiles) :

      m=re.match("(.*)\.[ab]",myfile0)
      if m :
         myfile=m.group(1)
      else :
         myfile=myfile0

      if filetype == "archive" :
         ab = abfile.ABFileArchv(myfile,"r")
      #elif filetype == "restart" :
      #   tmp = abfile.ABFileRestart(myfile,"r",idm=gfile.idm,jdm=gfile.jdm)
      else :
         raise NotImplementedError,"Filetype %s not implemented"%filetype

      if fieldname not in  ab.fieldnames :
         logger.error("Unknown field %s at level %d"%(fieldname,fieldlevel))
         logger.error("Available fields : %s"%(" ".join(ab.fieldnames)))
         raise ValueError,"Unknown field %s at level %d"%(fieldname,fieldlevel)

       
      fld = ab.read_field(fieldname,fieldlevel)



      figure = matplotlib.pyplot.figure(figsize=(8,8))
      ax=figure.add_subplot(111)
      #P=ax.pcolormesh(fld)
      #P=ax.pcolormesh(fld[2200:2800,3500:4500],cmap=cmap)
      P=ax.pcolormesh(fld,cmap=cmap)
      ax.figure.colorbar(P)
      if clim is not None : P.set_clim(clim)

      ax.set_title("%s:%s(%d)"%(myfile0,fieldname,fieldlevel))


      #figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=w5.min(), vmax=w5.max()))
      #ax.contour(w5)#,[-10.,-100.,-500.,-1000.])
      #ax.set_title("Slope fac in color, depth contours in black")
      #logger.info("Slope factor in slopefac.png")
      figure.canvas.print_figure("tst%03d.png"%i)



if __name__ == "__main__" :
   class ClimParseAction(argparse.Action) :
     def __call__(self, parser, args, values, option_string=None):
       tmp = values.split(",")
       tmp = [float(elem) for elem in tmp[0:2]]
       setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--clim',     action=ClimParseAction,default=None)
   parser.add_argument('--filetype',   type=str, help='',default="archive")
   parser.add_argument('--idm',        type=int, help='Grid dimension 1st index []', default=None)
   parser.add_argument('--jdm',        type=int, help='Grid dimension 2nd index []', default=None)
   parser.add_argument('fieldname',  type=str)
   parser.add_argument('fieldlevel', type=int)
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()

   main(args.filename,args.fieldname,args.fieldlevel,idm=args.idm,clim=args.clim,filetype=args.filetype)


