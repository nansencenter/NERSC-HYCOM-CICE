#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
import abfile
import modeltools.hycom
import modeltools.tools
import argparse
import matplotlib
import numpy
import logging
import datetime
import re
matplotlib.use('Agg')
import logging
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot
import cmocean
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
create spatio-varying thkdf4.[ab], veldf4.[ab], thkdf4.a
Usage:
   python ./hycom_creat_thkdf4_veldf4.py 10
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




def main(rmuwidth,exptid="") :

   
    #fh = Dataset('mld_dr003_l3.nc', mode='r')
    #fh_lons = fh.variables['lon'][:]
    #fh_lats = fh.variables['lat'][:]
    #fh_time = fh.variables['time'][:]
    #fh_mld_clim = fh.variables['mld_dr003_rmoutliers_smth_okrg'][:]
    #fh.close()
    
    grdpath='./'
    #gfile = abfile.ABFileGrid(jip+"regional.grid","r")
    gfile = abfile.ABFileGrid("regional.grid","r")
    idm=gfile.idm
    jdm=gfile.jdm
    #gfile = abfile.ABFileGrid("regional.grid","r")
    
    plon=gfile.read_field("plon")
    plat=gfile.read_field("plat")
    gfile.close()
    dpfile = abfile.ABFileBathy("regional.depth","r",idm=gfile.idm,jdm=gfile.jdm)
    depth=dpfile.read_field("depth")
    dpfile.close()
    # read from nc file
    
    #wdth=40
    wdth=rmuwidth
    mnn=0.02
    mxx=0.125
    thkdf4_array=numpy.zeros(depth.shape)
    thkdf4_array[:,:]=mnn
    ddd_mnn_mxx=np.linspace(mnn, mxx, num=wdth)
    ddd_mxx_mnn=numpy.flip(ddd_mnn_mxx[:])
    #print("ddd_E=",ddd_mnn_mxx[:])
    ##AA thkdf4_array[:,:]=0.015
    ##AA ddd=[(0.015+ii*0.0055) for ii in range(20)]
    #print "dd=",ddd[:]
    print("thkdf4_array.shape=",thkdf4_array.shape)
    for jjj in range(thkdf4_array.shape[0]) :
       thkdf4_array[jjj, 0:wdth] =ddd_mxx_mnn[:]
       thkdf4_array[jjj,-1*wdth:]=ddd_mnn_mxx[:]
    for iii in range(thkdf4_array.shape[1]):
       thkdf4_array[0:wdth  ,iii]=ddd_mxx_mnn[:]
       thkdf4_array[-1*wdth:,iii]=ddd_mnn_mxx[:]
    #   
    #fix lower left corner
    for iii in range(wdth):
       for jjj in  range(iii,wdth):
          thkdf4_array[jjj,iii]=ddd_mxx_mnn[iii]
          #print(thkdf4_array[jjj,iii])
    #fix upper left corner
    for iii in range(wdth):
       for jjj in  range(jdm-iii-1,jdm-wdth-1,-1):
          thkdf4_array[jjj,iii]=ddd_mxx_mnn[iii]
    #fix lower right corner
    LL = range(idm-wdth,idm)
    LL_r=LL[::-1]
    for iii in LL_r[:]:
       for jjj in  range(idm-iii-1,wdth):
          #print("jjj,iii=",jjj,iii)
          thkdf4_array[jjj,iii]=ddd_mxx_mnn[idm-iii-1]     
    #fix upper right corner
    for iii in LL_r[:]:
       for jjj in  range(jdm-1-(idm-iii-1),jdm-1-wdth,-1):
          thkdf4_array[jjj,iii]=ddd_mxx_mnn[idm-iii-1]

    #print '-----------------'
    #write thkdf4 to ab file
    af = abfile.AFile(plon.shape[1],plon.shape[0],"thkdf4.a","w")
    hmin,hmax = af.writerecord(thkdf4_array,None,record=0)
    af.close()
    print("thkdf4: range = %14.6e%14.6e\n"%(hmin,hmax))
    bf = open("thkdf4.b","w")
    bf.write("thkdf4: range = %14.6e%14.6e\n"%(hmin,hmax))
    bf.close()
    
    #print '-----------------'
    #write veldf4 to ab file
    af1 = abfile.AFile(plon.shape[1],plon.shape[0],"veldf4.a","w")
    hmin,hmax = af1.writerecord(thkdf4_array,None,record=0)
    af1.close()
    print("veldf4: range = %14.6e%14.6e\n"%(hmin,hmax))
    bf1 = open("veldf4.b","w")
    bf1.write("veldf4: range = %14.6e%14.6e\n"%(hmin,hmax))
    bf1.close()
   
   #print '-----------------'
    print("write to NC file------")
    ncfilename='thkdf4.nc'
    rootgrp = Dataset(ncfilename, "w", format="NETCDF4")
    logger.info("output to ncfile in  %s"%ncfilename)
    #dimension
    lat = rootgrp.createDimension("lat", gfile.jdm)
    lon = rootgrp.createDimension("lon", gfile.idm)
    #variable
    thkdf4_v = rootgrp.createVariable("thkdf4","f4",("lat","lon",))
    print('thkdf4_v.shape=',thkdf4_v.shape)
    thkdf4_v[:,:]=thkdf4_array[:,:]
    rootgrp.close()



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
   parser.add_argument('rmuwidth',     type=int, help='rum zone width,  example  10 []')
   parser.add_argument('--exptid'      , type=str,default="")
   args = parser.parse_args()

   main(exptid=args.exptid,rmuwidth=args.rmuwidth)



