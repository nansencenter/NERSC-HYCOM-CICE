#!/usr/bin/env python
import modeltools.hycom
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
#import modeltools.forcing.bathy
#import modeltools.hycom.io
import abfile
#import modeltools.cice.io
import numpy
from mpl_toolkits.basemap import Basemap
import netCDF4
import logging

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

if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('idm',     type=int, help='Grid dimension 1st index []')
   parser.add_argument('jdm',     type=int, help='Grid dimension 2nd index []')
   parser.add_argument('filename', help="")
   parser.add_argument('records',  nargs="+",    type=int)

   args = parser.parse_args()

   #afile = modeltools.hycom.io.AFile(args.idm,args.jdm,args.filename,"r")
   afile = abfile.AFile(args.idm,args.jdm,args.filename,"r")

   cmap=matplotlib.pyplot.get_cmap("jet")
   #cmap=matplotlib.pyplot.get_cmap("YlOrRd")

    
   for record in args.records :

      fld = afile.read_record(record)
      figure = matplotlib.pyplot.figure(figsize=(8,8))
      ax=figure.add_subplot(111)
      #P=ax.pcolormesh(fld)
      #P=ax.pcolormesh(fld[2200:2800,3500:4500],cmap=cmap)
      P=ax.pcolormesh(fld,cmap=cmap)
      ax.figure.colorbar(P)


      #figure.colorbar(P,norm=matplotlib.colors.LogNorm(vmin=w5.min(), vmax=w5.max()))
      #ax.contour(w5)#,[-10.,-100.,-500.,-1000.])
      #ax.set_title("Slope fac in color, depth contours in black")
      #logger.info("Slope factor in slopefac.png")
      figure.canvas.print_figure("tst%03d.png"%record)




