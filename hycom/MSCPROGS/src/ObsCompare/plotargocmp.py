#!/local/python-2.5.1-login-gnu/bin/python
import pylab as p
import re
from numpy import *
import os
import glob
#import matplotlib
#matplotlib.use('Agg')

#for fname in ( 'profile_date20080619_id0006900494_pos-18.21Ex44.67N' ) :
#if (True) :
#for fname in os.listdir("./*N") :
for fname in glob.glob("./*N") :
   #fname = 'profile_date20080620_id0004900438_pos-72.57Ex26.84N'
   print fname
   subs=fname.split("_")
   subs2=subs[3].split("pos")
   pos=subs2[1]
   title =  subs[1] + " Position:" + subs2[1] + "Float ID: " + subs[2]


   depths=[];
   msaln=[];
   mtemp=[];
   osaln=[];
   otemp=[];
   file = open(fname)
   for line in file:
      subs=line.split()
      depths.append(-float(subs[0]))
      osaln.append (float(subs[1]))
      msaln.append (float(subs[2]))
      otemp.append (float(subs[3]))
      mtemp.append (float(subs[4]))
   file.close()
   msaln=array(msaln)
   osaln=array(osaln)
   mtemp=array(mtemp)
   otemp=array(otemp)
   depths=array(depths)

   I=where(msaln != float(-999))
   J=where(osaln != float(-999))
   #print prod(I[0].shape), prod(J[0].shape), prod(msaln.shape), prod(osaln.shape)
   if ( prod(I[0].shape) > 1 and prod(J[0].shape) > 1 ) :
      p.hold(True)
      leg=[]
      if ( prod(I[0].shape) > 1 ) :
         p.plot(msaln[I],depths[I], color='r', linewidth=2)
         leg.append("Model")
      if ( prod(J[0].shape) > 1 ) :
         p.plot(osaln[J],depths[J], color='b', linewidth=2)
         leg.append("Argo")
      p.grid(True)
      p.ylabel("Depths[m]")
      p.xlabel("Salinity[psu]")
      p.title(title)
      p.legend(leg)
      p.savefig("saltprofile_" + pos + ".png")
      print "---saltprofile_" + pos + ".png"
      p.close()

   I=where(mtemp != -999)
   J=where(otemp != -999)
   #print prod(I[0].shape), prod(J[0].shape), prod(mtemp.shape), prod(otemp.shape)
   if ( prod(I[0].shape) > 1 and prod(J[0].shape) > 1 ) :
      p.hold(True)
      leg=[]
      if ( prod(I[0].shape) > 1 ) :
         p.plot(mtemp[I],depths[I], color='r', linewidth=2)
         leg.append("Model")
      if ( prod(J[0].shape) > 1 ) :
         p.plot(otemp[J],depths[J], color='b', linewidth=2)
         leg.append("Argo")
      p.grid(True)
      p.ylabel("Depths[m]")
      p.xlabel("Temperature[C]")
      p.title(title)
      p.legend(leg)
      p.savefig("tempprofile_" + pos + ".png")
      print "---tempprofile_" + pos + ".png"
      p.close()
