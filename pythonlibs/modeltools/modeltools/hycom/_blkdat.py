""" MNodule for parsing hycom blkdat """
import numpy 
import sys
import logging
import re


# TODO: Add writer

class BlkdatError(Exception):
    """Base class for exceptions in this module."""
    pass

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch) 



class BlkdatParser(object) :
   """ Class for doing binary input/output on hycom .a files """
   _integer_fields = [ "iversn","iexpt","idm","jdm","itest","jtest","kdm",
      "nhybrd","nsigma","kapref","thflag","vsigma","iniflg","jerlv0","yrflag",
      "sshflg","incflg","incstp","incupf","hybmap","hybflg","advflg",
      "advtyp","momtyp","ishelf","ntracr","trcflg","tsofrq","mlflag","pensol",
      "dypflg","bblkpp","shinst","dbdiff","nonloc","bodiw","difout","difsmo","hblflg",
      "niter","langmr","clmflg","fltflg","nfladv","wndflg","ustflg","flxflg","empflg",
      "dswflg","albflg","sssflg","sstflg","lwflag","icmflg",",slprs","stroff","flxoff",
      "flxsmo","relax","trcrlx","priver","epmass"]
   _float_fields = [ "dp00", "dp00x", "dp00f", "ds00", "ds00x", "ds00f", "isotop", "saln0",
                     "cplifq", "slip", "baclin", "batrop","dp00i","nestfq","bnstfq"]
      
   def __init__(self,filename) :

      self._sigma    = []
      self._dp0k     = []
      self._ds0k     = []
      self._datadict = {}
      self._descdict = {}
      fid=open(filename,"r")

      # 4 first line is header
      self.header=[]
      self.header.append(fid.readline())
      self.header.append(fid.readline())
      self.header.append(fid.readline())
      self.header.append(fid.readline())


      for line in fid :

         # general pattern to search for
         #m = re.match("^(.*)'([a-z0-9]{6})'[ ]*=[ ]*(.*$)",line.strip())
         m = re.match("^[\s]*(.*)[\s]*'([a-z0-9 _]{6})'[\s]*=[\s]*(.*$)",line.strip())
         if not m :
            print m,line
            raise BlkdatError,"Error when parsing file. Line does not match required pattern"

         value=m.group(1).strip()
         key  =m.group(2).strip()
         desc =m.group(3).strip()

         #print key
         if key == "sigma" :
            self._sigma.append(float(value))
         elif key == "ds0k" :
            self._ds0k.append(float(value))
         elif key == "dp0k" :
            self._dp0k.append(float(value))
         else  :
            self._datadict[key] = value

         if key in  self._integer_fields :
            self._datadict[key] = int(self._datadict[key])
         elif key in  self._float_fields :
            self._datadict[key] = float(self._datadict[key])

      # Create ds0k and dp0k from info in blkdat.input
      self._ds0k= self._ds0k_profile()
      self._dp0k= self._dp0k_profile()


   def __getitem__(self,i) :
      if i=="sigma" :
         return self._sigma
      elif i=="dp0k" :
         return self._dp0k
      elif i=="ds0k" :
         return self._ds0k
      elif i in self._integer_fields :
         return int(self._datadict[i])
      else :
         return self._datadict[i]


   def _dp0k_profile(self) :
      if self["dp0k"] :
         dp0k=self["dp0k"]
      else  :
         # Create dp0k (deep z-level) from parameters
         dp00=self["dp00"]
         dp00x=self["dp00x"]
         dp00f=self["dp00f"]
         dp0k=[]
         for k in range(self["kdm"]) :
            dp0k.append(dp00*dp00f**k)
         dp0k=[min(elem,dp00x) for elem in dp0k]
         dp0k=numpy.array(dp0k)
      dp0k[self["nhybrd"]:]=self["dp00i"]
      return dp0k


   def _ds0k_profile(self) :
      # Check for "ds0k" 
      if self["ds0k"] :
         ds0k=self["ds0k"]
      else  :
         # Create ds0k (shallow z-level)  from parameters
         ds00=self["ds00"]
         ds00x=self["ds00x"]
         ds00f=self["ds00f"]
         ds0k=[]
         for k in range(self["kdm"]) :
            ds0k.append(ds00*ds00f**k)
         ds0k=[min(elem,ds00x) for elem in ds0k]
         ds0k=numpy.array(ds0k)
      ds0k[self["nsigma"]:]=self["dp00i"]
      ds0k[self["nhybrd"]:]=self["dp00i"]
      return ds0k

   @property
   def intf_deep(self) :
      """ Min layer thickness for deep z-level """
      dp0k=self["dp0k"]
      intf0k=numpy.zeros((self["kdm"]+1))
      for i in range(self["nhybrd"]) :
         intf0k[i+1] = intf0k[i] + dp0k[i]
      return intf0k

   @property
   def intf_shallow(self) :
      """ Min layer thickness for shallow z-level """
      ds0k=self["ds0k"]
      intf0s=numpy.zeros((self["kdm"]+1))
      for i in range(self["nsigma"]) :
         intf0s[i+1] = intf0s[i] + ds0k[i]
      intf0s[self["nsigma"]:] = intf0s[self["nsigma"]-1]
      return intf0s


   def intf_min_profile(self,bottom) :
      """ From an input layer profile, set up min thickness interfaces appropriate for pathymetry """
      kdm   =self["kdm"]
      nsigma=self["nsigma"]
      nhybrd=self["nhybrd"]
      intf_shallow = self.intf_shallow
      intf_deep    = self.intf_deep
      intf=numpy.zeros(tuple(list(bottom.shape,)+[self["kdm"]+1]))
      logger.info("shallow z-level interface[nsigma+1]: %5.1f"%intf_shallow[nsigma])
      logger.info("Deep    z-level interface[nsigma+1]: %5.1f"%intf_deep[nsigma])

      # Depth of deepest sigma interface, compared to bottom
      ideep        = intf_deep   [self["nsigma"]]
      ishallow     = intf_shallow[self["nsigma"]]
      f_ishallow   = ishallow/numpy.maximum(bottom,1e-4)
      f_ideep      = ideep/numpy.maximum(bottom,1e-4)
      msk0=bottom<1e-4

      # Case 1) Fixed shallow z_level (nsigma shallow z levels extend beyond ocean floor)
      Imask=f_ishallow>=1. 
      I=numpy.where(Imask)
      intf[Imask,:nsigma+1] = intf_shallow[:nsigma+1]
      intf[Imask,nsigma+1:] = intf_shallow[nsigma]

      # Case 2) Fixed deep z_level (nsigma deep z levels above ocean floor)
      Jmask=f_ideep<=1.
      J=numpy.where(Jmask)
      intf[Jmask,:nhybrd+1] = intf_deep[:nhybrd+1]
      intf[Jmask,nhybrd+1:] = intf_deep[nhybrd]

      # Case 3) Sigma coordinates where sigma-th deep z levels is below ocean floor, and sigma-th shallow z level is above ocean floor
      itest=52
      ktest=10

      # Squeeze deep layers. Ooook... Thought this should work, but apparently not always?! 
      #Kmask=numpy.logical_and(~Jmask,~Imask)
      #intf[Kmask,:nsigma+1] = intf_deep[:nsigma+1]
      #tmp        = numpy.transpose(intf[Kmask,:])/f_ideep[Kmask]
      #tmp[nsigma+1:,] = tmp[nsigma,:]
      #intf[Kmask,:] = tmp.transpose()
      #
      # Stretch shallow layers. This works 
      Kmask=numpy.logical_and(~Jmask,~Imask)
      intf[Kmask,:nsigma+1] = intf_shallow[:nsigma+1]
      tmp        = numpy.transpose(intf[Kmask,:])/f_ishallow[Kmask]
      tmp[nsigma+1:,] = tmp[nsigma,:]
      intf[Kmask,:] = tmp.transpose()

      # Interfaces are 0 where depth is 0
      intf[msk0,:]=0.

      #intf = numpy.transpose(numpy.minimum(numpy.transpose(intf),bottom))
      intf = numpy.transpose(numpy.minimum(numpy.transpose(intf),bottom.transpose()))
      return intf,{"Shallow z":Imask,"Deep z":Jmask,"Sigma":Kmask}




          











