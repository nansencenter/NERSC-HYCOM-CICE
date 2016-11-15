#!/usr/bin/env python
##!/usr/bin/python -E
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import modeltools.hycom
import logging
import argparse
import datetime
import numpy
import os
import netCDF4
import re
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

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


def main(blkdatfiles):

   x      = numpy.linspace(0., 30.*1000,120)
   nx=x.size
   bottom=numpy.zeros((nx))
   bottom[:nx/2] = numpy.linspace(5., 100.,nx/2) + 80.*numpy.sin(x[0:nx/2] * 2*numpy.pi / 20000.)
   bottom[nx/2:] = numpy.linspace(bottom[nx/2-1], 3000+500.,nx/2)
   bottom[nx/2:] = bottom[nx/2:] + 250.*numpy.sin((x[nx/2:]-x[nx/2-1]) * 2*numpy.pi / 20000.)

   f,ax = plt.subplots(2,2,sharey=False)
   ax1=ax[0,0]
   ax2=ax[1,0]
   ax3=ax[0,1]
   cols=["r","g","b","m","k","c","y"]

   f2,ax9 = plt.subplots(1,figsize=(10,5))



   for i,blkdatfile in enumerate(blkdatfiles) :

      # Open blkdat.input
      bp=modeltools.hycom.BlkdatParser(blkdatfile)
      dp0k = bp["dp0k"]
      ds0k = bp["ds0k"]
      sigma = bp["sigma"]
      kdm=bp["kdm"]
      nhybrd=bp["nhybrd"]
      nsigma=bp["nsigma"]
      thflag=bp["thflag"]
      mysig   = modeltools.hycom.Sigma(thflag)

      intf,masks = bp.intf_min_profile(numpy.array([3000.]))
      intf=numpy.squeeze(intf)
      intfmid=(intf[1:]+intf[:-1])*.5

      # Plot vertical profile of sig-2 vs depth amd index ?
      ax1.plot(sigma,-intfmid,lw=2,label=os.path.basename(blkdatfile),color=cols[i])
      ax1.set_ylabel("Min interface depth")
      ax1.set_xlabel("sigma")
      #
      ax2.plot(sigma,-numpy.arange(intfmid.size),lw=2,color=cols[i])
      ax2.set_ylabel("Layer index")
      ax2.set_xlabel("sigma")
      #
      ax3.plot(sigma,-intfmid,lw=2,label=os.path.basename(blkdatfile),color=cols[i])
      ax3.set_ylabel("Min interface depth")
      ax3.set_xlabel("sigma")
      ax3.set_ylim(-500,0)

      intf,masks = bp.intf_min_profile(numpy.array(bottom))

      # Plot layers
      first=False
      for k in range(intf.shape[1]) :
         if (k+1)%5 == 0 :
            ax9.plot(x,-intf[:,k],color=cols[i],linestyle="--",lw=.5)
         else:
            if first :
               label=os.path.basename(blkdatfile)
               first=False
            else :
               label=None
            ax9.plot(x,-intf[:,k],color=cols[i],lw=.5,label=label)
         if k>=1 :
            xpos = int(x.size*.7 + float(i)*x.size*.1/len(blkdatfiles))
            textx = x[xpos]
            texty = -0.5*(intf[xpos,k-1] + intf[xpos,k])
            #print k,textx,texty,intf[xpos,k-1],intf[xpos,k]
            ax9.text(textx,texty,str(k),verticalalignment="center",horizontalalignment="center",fontsize=6,color=cols[i])
   ax9.plot(x,-bottom,lw=4,color="k")



   # 
   f.savefig("mindp_sigma.png",dpi=180)
   f2.savefig("mindp_sigma_sec.png",dpi=180)
   ax9.set_ylim([-500.,0.])
   f2.savefig("mindp_sigma_sec_500m.png",dpi=180)
   ax9.set_ylim([-100.,0.])
   f2.savefig("mindp_sigma_sec_100m.png",dpi=180)








if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('blkdatfile' , nargs="+", help="")

   args = parser.parse_args()
   main(args.blkdatfile)



