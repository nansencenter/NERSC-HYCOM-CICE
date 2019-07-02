#import scipy.io.netcdf
import netCDF4
import numpy
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


#default_threshold=-5.


class GEBCO2014(object) :

   #def __init__(self,filename="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D_median2km.nc") :
   def __init__(self,filename="/work/shared/nersc/msc/ModelInput/bathymetry/GEBCO_2014/GEBCO_2014_2D") :
      self._filename = filename
      #self._nc = scipy.io.netcdf.netcdf_file(self._filename)
      self._nc = netCDF4.Dataset(self._filename,"r")

      self._lon = self._nc.variables["lon"][:]
      self._lat = self._nc.variables["lat"][:]
      #self._elevation = self._nc.variables["elevation"][:]
      self._elevation = self._nc.variables["z"][:]


      self._dlon= ( self._lon[-1] - self._lon[0] ) / (self._lon.size-1)
      self._dlat= ( self._lat[-1] - self._lat[0] ) / (self._lat.size-1)


   def regrid(self,lon,lat,width=None,myfilter="median"):

      # Pivot points
      wi= numpy.mod(lon-self._lon[0],360.) / self._dlon
      wj= (lat-self._lat[0]) / self._dlat

      i = numpy.floor(wi).astype("int")
      j = numpy.floor(wj).astype("int")
      ip1 = numpy.mod(i+1,self._lon.size)
      jp1 = numpy.minimum(j+1,self._lat.size-1)

      #print i
      #print wi
      wi = wi - i
      wj = wj - j



      if width is not None :

         numpts = int(numpy.round(width / (self._dlat*111*1000 )))
         if numpts < 3 :
            logger.debug("Too few points to do filtering. Using raw bathymetry")
            w = (1.-wi)*(1.-wj)*self._elevation[j  ,i  ] + \
                (   wi)*(1.-wj)*self._elevation[j  ,ip1] + \
                (   wi)*(   wj)*self._elevation[jp1,ip1] + \
                (1.-wi)*(   wj)*self._elevation[jp1,i  ]

         else :
            numpts = numpts/2 # -> half width
            logger.debug("Filtering width half-width = %d grid cells"%numpts)
            wtmp = numpy.zeros( ( (2*numpts)**2,lon.shape[0],lon.shape[1] ) )
            cnt=0
            for k in range(-numpts,numpts) :
               for l in range(-numpts,numpts) :
                  tmpi = i + k
                  tmpi = numpy.mod(tmpi,self._lon.shape[0])
                  tmpj = j + l
                  I = numpy.where(tmpj < 0)
                  if I[0].size>0 : 
                     tmpj[I] = numpy.abs(tmpj[I])
                     tmpi[I] = numpy.mod(tmpi[I]+ self._lon.size/2,self._lon.shape[0])
                  I = numpy.where(tmpj >= self._lat.size)
                  if I[0].size>0 : 
                     tmpj[I] = self._lat.size - 1 -  (tmpj[I] - (self._lat.size - 1))
                     tmpi[I] = numpy.mod(tmpi[I]+self._lon.size/2,self._lon.shape[0])
                  #print cnt,wtmp.shape
                  #print k,l,numpts,cnt,wtmp.shape
                  wtmp[cnt,:] = self._elevation[tmpj,tmpi]
                  cnt+=1
            if myfilter == "median" :
               w = numpy.median(wtmp,axis=0)
            elif myfilter == "average" :
               w = numpy.mean(wtmp,axis=0)
            else :
               w = numpy.mean(wtmp,axis=0)


      else :
         logger.debug("Using raw bathymetry")
         w = (1.-wi)*(1.-wj)*self._elevation[j  ,i  ] +\
             (   wi)*(1.-wj)*self._elevation[j  ,ip1] +\
             (   wi)*(   wj)*self._elevation[jp1,ip1] +\
             (1.-wi)*(   wj)*self._elevation[jp1,i  ]


      return w
