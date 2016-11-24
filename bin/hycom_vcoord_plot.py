#!/usr/bin/env python
##!/usr/bin/python -E
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches
import modeltools.hycom
import logging
import argparse
import datetime
import numpy
import os

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

def my_fill_between(x, y1, y2=0, ax=None, **kwargs):
    """Plot filled region between `y1` and `y2`.

    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.
    """
    ax = ax if ax is not None else plt.gca()
    label=kwargs["label"]
    del kwargs["label"]
    ax.fill_between(x, y1, y2, **kwargs)
    if "where" in kwargs.keys() : del kwargs["where"]
    del kwargs["interpolate"]
    #print kwargs
    kwargs["label"]=label
    p = plt.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)
    return p


def main(blkdat_file):

   bp=modeltools.hycom.BlkdatParser(blkdat_file)

   kdm=bp["kdm"]
   nhybrd=bp["nhybrd"]
   nsigma=bp["nsigma"]
   logger.info("%d layers, %d hybrid layers, %d sigma layers"%(kdm,nhybrd,nsigma))

   # Get min interface thickness profile 
   #tmp,tmpmasks=bp.intf_min_profile(numpy.array([10000]))
   #maxdepth=numpy.max(tmp)
   maxdepth=3000.

   # x and depth arrays
   x      = numpy.linspace(0., 30.*1000,120)
   nx=x.size
   bottom=numpy.zeros((nx))
   bottom[:nx/2] = numpy.linspace(5., 100.,nx/2) + 80.*numpy.sin(x[0:nx/2] * 2*numpy.pi / 20000.)
   bottom[nx/2:] = numpy.linspace(bottom[nx/2-1], maxdepth+500.,nx/2)
   bottom[nx/2:] = bottom[nx/2:] + 250.*numpy.sin((x[nx/2:]-x[nx/2-1]) * 2*numpy.pi / 20000.)

   # Get a min interface section and masks describing "regime"
   intf,masks=bp.intf_min_profile(bottom)
   Imask=masks["Shallow z"]
   Jmask=masks["Deep z"]
   Kmask=masks["Sigma"]

   ax=plt.gca()
   ax.hold(True)

   # Plot colors for regimes
   my_fill_between(x,intf[:,nsigma]*-1.,0.,ax,color="g",interpolate=False,label="Sigma")
   my_fill_between(x,intf[:,nsigma]*-1.,0.,ax,color="r",interpolate=False,where=Imask,label="Shallow z")
   my_fill_between(x,intf[:,nhybrd]*-1.,0.,ax,color="b",interpolate=False,where=Jmask,label="Deep z")
   #my_fill_between(x,intf[:,nsigma]*-1.,0.,ax,color="g",interpolate=False,where=Kmask,label="Sigma")
   if nhybrd <> kdm :
      my_fill_between(x,-bottom,intf[:,nhybrd]*-1.,ax,color="c",interpolate=False,where=~Jmask,label="Isopycnal")

   # Plot layers
   for k in range(intf.shape[1]) :
      #if (k+1)%5 == 0 :
      #   plt.plot(x,-intf[:,k],color="k",linestyle="--",label=str(k+1))
      #else:
      #   plt.plot(x,-intf[:,k],color=".5")
      plt.plot(x,-intf[:,k],color=".5")

      if k>=1 :
         xpos = int(x.size*.8)
         textx = x[xpos]
         texty = -0.5*(intf[xpos,k-1] + intf[xpos,k])
         ax.text(textx,texty,str(k),verticalalignment="center",horizontalalignment="center",fontsize=6,fontweight="bold")
         #
         xpos = int(x.size*.2)
         textx = x[xpos]
         texty = -0.5*(intf[xpos,k-1] + intf[xpos,k])
         ax.text(textx,texty,str(k),verticalalignment="center",horizontalalignment="center",fontsize=6,fontweight="bold")

      if (k+1)%2 == 0 :
         pc=ax.fill_between(x,-intf[:,k-1],-intf[:,k],color="none")
         for path in pc.get_paths()  :
            patch = matplotlib.patches.PathPatch(path, hatch='//', facecolor='none',linewidth=.1)
            ax.add_patch(patch)


   # Save to full-depth file
   ax.plot(x,-bottom,lw=4,color="k")
   ax.legend(fontsize=6)
   fname="vcoord.png"
   logger.info("Min thickness section in %s"%fname)
   plt.gcf().savefig(fname,dpi=180)
      
   # Save to top 200 m depth file
   ax.set_ylim(-200,0)
   fname="vcoord200.png"
   logger.info("Min thickness section [top 200 m] in %s"%fname)
   plt.gcf().savefig(fname,dpi=180)
      
   # Save to top 200 m depth file
   ax.set_ylim(-50,0)
   fname="vcoord050.png"
   logger.info("Min thickness section [top 50 m] in %s"%fname)
   plt.gcf().savefig(fname,dpi=180)


if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='Creates an illustration of the setting of the hybrid coordinates, as specified in blkdat.input')
   parser.add_argument('blkdat_file' , help="blkdat file")
   args = parser.parse_args()
   main(args.blkdat_file)



