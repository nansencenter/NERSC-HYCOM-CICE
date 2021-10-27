#!/usr/bin/env python
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
##AO from mpl_toolkits.basemap import Basemap
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


# Northern line 
northern_limit_longitudes = [-35 ,-16, -6.6, -4.2, -0.3,   4., 31.6, 97., 172.,-172,
   -150., -129, -101.,-88, -83.84, -82.7, -75., -46.]
northern_limit_latitudes  = [67.5, 65,  62.,  57.7, 51.0, 49.4, 51.5, 57., 66.7, 66.1,
   65.,    64.,  64., 66.6, 67.81,  73.7, 81.9, 80.5]


def main(blkdat1,blkdat2,nsmooth=0) :

   # Open blkdat files 
   bp1=modeltools.hycom.BlkdatParser(blkdat1)
   bp2=modeltools.hycom.BlkdatParser(blkdat2)

   sig1 = bp1["sigma"]
   sig2 = bp2["sigma"]

   # Read plon,plat
   gfile=abfile.ABFileGrid("regional.grid","r")
   plon=gfile.read_field("plon")
   plat=gfile.read_field("plat")
   gfile.close()

   # Define area with sig2


   # smooth transition (sets weights)

   # Go through layers, and apply weight for different sigma regions




   # Read input bathymetry
   bfile=abfile.ABFileBathy("regional.depth","r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   in_depth_m=bfile.read_field("depth")
   bfile.close()

   # Starting point  (State 1 for Atlantic)
   kapref=numpy.ones(plat.shape)*1.0 
   print kapref.min(),kapref.max()

   # Find regions north of northern limit. Assumes 
   for segment in range(len(northern_limit_longitudes)) :

      ind1=segment
      ind2=(segment+1)%len(northern_limit_longitudes)

      lo1=northern_limit_longitudes[ind1]
      la1=northern_limit_latitudes [ind1]
      lo2=northern_limit_longitudes[ind2]
      la2=northern_limit_latitudes [ind2]

      
      tmp1=numpy.mod(plon+360-lo1,360.)
      tmp2=numpy.mod(lo2+360-lo1,360.)
      J=tmp1<=tmp2
      #print numpy.count_nonzero(J)

      # Linear weights and latitude in selected points
      w2 = tmp1 / tmp2
      w1 = 1.-w2
      la = la2 * w2 + la1 * w1
      
      kapref[J] = numpy.where(plat[J] > la[J],2.0,kapref[J])


   import scipy.ndimage
   kapref = scipy.ndimage.gaussian_filter(kapref, sigma=20)

   #print in_depth_m.min(),type(in_depth_m)
   kaprefplot=numpy.ma.masked_where(in_depth_m.mask,kapref)
   figure = matplotlib.pyplot.figure()
   ax=figure.add_subplot(111)
   P=ax.pcolormesh(kaprefplot)
   figure.colorbar(P,ax=ax)
   figure.canvas.print_figure("kapref.png")

   af = abfile.AFile(plon.shape[1],plon.shape[0],"tbaric.a","w")
   hmin,hmax = af.writerecord(kapref,None,record=0)
   af.close()
   bf = open("tbaric.b","w")
   bf.write("tbaric.b\n")
   bf.write("\n")
   bf.write("\n")
   bf.write("\n")
   bf.write("i/jdm =  %5d %5d\n"%(plon.shape[1],plon.shape[0]))
   bf.write("tbaric: range = %14.6e%14.6e\n"%(hmin,hmax))
   bf.close()





if __name__ == "__main__" :

   raise NameError,"in testing"


   parser = argparse.ArgumentParser(description='Prepare spatially varying sigma coordinate setup')
   parser.add_argument('--nsmooth',     action=ClimParseAction,default=None)
   parser.add_argument('blkdat1',     type=int, help='')
   parser.add_argument('blkdat2',     type=int, help='')
   #parser.add_argument('lon1',        type=float, help='')
   #parser.add_argument('lat1',        type=float, help='')
   #parser.add_argument('lon2',        type=float, help='')
   #parser.add_argument('lat2',        type=float, help='')
   main(blkdat1,blkdat2)
