""" Module for parsing hycom blkdat """
import datetime
import numpy 
import sys
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

def leapyear(iyr) :
   if ( iyr%4==0 and 1901+iyr%400==0) or iyr%4!=0  :
     leapyear=False
   else :
     leapyear=True
   return leapyear


def datetime_to_ordinal(dt,yrflag) :
   if yrflag != 3 :
      raise ValueError("Only yearflag=3 supported for datetime_to_ordinal")
   else :
      td = dt - datetime.datetime(dt.year,1,1,0,0,0) 
      tdd = td.days + 1
      #tdh = td.seconds/86400.
      tdh = td.seconds/3600.
      tds = td.seconds - (tdh*3600)
   return tdd,tdh,tds


def dayfor(iyear,iday,ihour,yrflag) :
   """ Converts hycom ordinal day (year, day of year, hour) to model day """

   # 360 days per model year, starting Jan 16
   if  yrflag == 0 :
     dtime =  (iyear-1) * 360.0 +  \
              iday              +  \
              ihour/24.0        -  \
              16.0
   # 366 days per model year, starting Jan 16
   elif yrflag == 1 :
     dtime =  (iyear-1) * 366.0 + \
              iday              + \
              ihour/24.0        - \
              16.0
   # 366 days per model year, starting Jan 01
   elif yrflag == 2 :
     dtime =  (iyear-1) * 366.0 + \
              iday              + \
              ihour/24.0        - \
              1.0
   # model day is calendar days since 01/01/1901
   elif yrflag == 3 :

     if iyear < 1901 :
        raise ValueError('error in forday - yrflag value %d must have year >=1901. year is %d'%(yrflag,iyear))


     iyr=iyear-1901
     dtime=0.
     for k in range(1,iyr+1) :
        if leapyear(k) :
           dtime=dtime+366.0
        else :
           dtime=dtime+365.0
     iyr=iyr+1

     iday2=iday
     while iday2+ihour/24.>0.0 :

        if leapyear(iyr) : 
           diy=366.0
        else :
           diy=365.0
        
        if ((iday2+(ihour/24.))/diy>=1.0):
           dtime=dtime+diy
        else :
           dtime=dtime+iday2+ihour/24.0

        iday2=iday2-diy

   else :
      raise ValueError('error in forday - unsupported yrflag value: %d'%yrflag)

   return dtime


def dayfor_datetime(iyear,iday,ihour,yrflag) :
   dtime=dayfor(iyear,iday,ihour,yrflag) 
   if yrflag == 3 :
      tmp =datetime.datetime(1901,1,1,0,0,0)
      return tmp + datetime.timedelta(days=dtime)
   else :
      raise NotImplementedError("yrflag = %d"%yrflag)



def forday(dtime,yrflag) :
   """ converts model day to "calendar" date (year,ordinal-day,hour). """

   # ---   360 days per model year, starting Jan 16
   if  yrflag == 0 :
      iyear =  int((dtime+15.0010)/360.0) + 1
      iday  =  int(numpy.mod( dtime+15.0010 ,360.0) + 1)
      ihour = int((numpy.mod( dtime+15.0010 ,360.0) + 1.0 - iday)*24.0)
   # ---   366 days per model year, starting Jan 16
   elif yrflag == 1 :
      iyear =  int((dtime+15.001)/366.) + 1
      iday  =  int(numpy.mod( dtime+15.001 ,366.) + 1)
      ihour = int((numpy.mod( dtime+15.001 ,366.) + 1. - iday)*24.)
   # ---   366 days per model year, starting Jan 01
   elif yrflag == 2 :
      iyear =  int((dtime+ 0.001)/366.) + 1
      iday  =  int(numpy.mod( dtime+ 0.001 ,366.) + 1)
      ihour =  int((numpy.mod( dtime+ 0.001 ,366.) + 1. - iday)*24.)
   # ---   model day is calendar days since 01/01/1901
   elif yrflag ==3 :

      iyr   = int((dtime-1.)/365.25)
      nleap = int(iyr/4)
      dtim1 = 365.*iyr + nleap + 1.
      day   = dtime - dtim1 + 1.
      if  dtim1 > dtime :
        iyr = iyr - 1
      elif day>=367.0 :
        iyr = iyr + 1
      elif day>=366.0 and iyr%4!=3 :
        iyr = iyr + 1


      nleap = int(iyr/4)
      dtim1 = 365.0*iyr + nleap + 1.0

      iyear =  1901 + iyr
      iday  =  int(dtime - dtim1 + 1.001)
      ihour = int((dtime - dtim1 + 1.001 - iday)*24.0)

   else :
      raise ValueError('error in forday - unsupported yrflag value: %d'%yrflag)

   return iyear, iday, ihour



def forday_datetime(dtime,yrflag) :
   print(dtime)
   iy,id,ih=forday(dtime,yrflag) 
   print(iy,id,ih)
   if yrflag == 3 :
      tmp =datetime.datetime(iy,1,1,0,0,0)
      return tmp + datetime.timedelta(days=id-1) + datetime.timedelta(hours=ih)
   else :
      raise NotImplementedError("yrflag = %d"%yrflag)








