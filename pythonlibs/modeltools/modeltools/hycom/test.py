#!/usr/bin/env python
import unittest
#import timetools
import random
import numpy
import modeltools.hycom

class TimeTest(unittest.TestCase):

    def test_time_conversion(self):
       """ Tests that forday and dayfor are consistent """
       numtest=20000
       for k in range(numtest) :
          yrflag = random.randrange(0,4)
          selector = random.randrange(0,3)
          dtime = random.randrange(1,200000)
          
          # Tests wnen close to integer, and normal
          if selector==0 :
             dtime = dtime + random.random()/100.
          elif selector==1 :
             dtime = dtime + 1.-random.random()/100.
          elif selector==2 :
             dtime = dtime + random.random()   

          if k%50000 == 0 : 
              print "Doing test %10d of %10d. yrflag=%d, dtime=%20.16g"%(k,numtest,yrflag,dtime)

          try :
              #print "dtime,yrflag=",dtime,yrflag
              iyear,iday,ihour = modeltools.hycom.forday(dtime,yrflag)
              #print "iyear,iday,ihour=",iyear,iday,ihour
              dtimeout=modeltools.hycom.dayfor(iyear,iday,ihour,yrflag)
              #print "dtime, reversed:",dtime
          except:
              unittest.fail("Caught error. dtime = %20.10g,yrflag=%d"%(dtime,yrflag))
          
          #print yrflag,dtime,dtimeout-dtime,3600./86400.
          
          #dtimeout is "floored" to nearest hour. So error should never be greater than this:
          if (dtime-dtimeout)>3600./86400:
              unittest.fail("Error: inn forday or dayfor, dtime=%20.10g, yrflag%d="%(dtime,yrflag))

class HycomIOTest(unittest.TestCase) :

   def test_afile_writeread_nomask(self) :

      idm = random.randrange(10,5000)
      jdm = random.randrange(10,5000)
      print "Creating %dX%d random arrays"%(idm,jdm)
      scale=1e4
      wfld1=numpy.random.rand(jdm,idm)*scale
      wfld2=numpy.random.rand(jdm,idm)*scale
      wfld3=numpy.random.rand(jdm,idm)*scale

      afile=modeltools.hycom.AFile(idm,jdm,"test.a","w")
      afile.writerecord(wfld1,None,record=None)
      afile.writerecord(wfld2,None,record=None)
      afile.writerecord(wfld3,None,record=None)
      afile.close()

      afile=modeltools.hycom.AFile(idm,jdm,"test.a","r")
      rfld1=afile.readrecord(None,record=0)
      rfld2=afile.readrecord(None,record=1)
      rfld3=afile.readrecord(None,record=2)
      afile.close()

      maxdiff1=numpy.abs(numpy.amax(rfld1-wfld1))/scale
      maxdiff2=numpy.abs(numpy.amax(rfld2-wfld2))/scale
      maxdiff3=numpy.abs(numpy.amax(rfld3-wfld3))/scale

      if all([elem < 1e-7 for elem in [maxdiff1,maxdiff2,maxdiff3]]) :
         print "afile io test passed"
      else :
         unittest.fail("AFile IO failed. MAx diff between read/written: %14.7g"%max([maxdiff1,maxdiff2,maxdiff3]))
      #print "end"

   def test_abfilebathy_writeread_nomask(self) :
      idm = random.randrange(10,5000)
      jdm = random.randrange(10,5000)
      scale=1e4

      print "Creating %dX%d bathy array"%(idm,jdm)
      wfldout=scale + numpy.random.rand(jdm,idm)*scale
      bathyfile=modeltools.hycom.ABFileBathy("testbathy","w")
      bathyfile.write_field(wfldout,None)
      bathyfile.close()

      print "Reading %dX%d bathy array"%(idm,jdm)
      bathyfile=modeltools.hycom.ABFileBathy("testbathy","r",idm=idm,jdm=jdm)
      wfldin=bathyfile.read_field("depth",None)
      bathyfile.close()

      fldmaxdiff=numpy.abs(numpy.amax(wfldin-wfldout))/scale
      bmin,bmax = bathyfile.bminmax("depth")
      amax = numpy.amax(wfldin)
      amin = numpy.amin(wfldin)
      abmindiff=numpy.abs(amin-bmin)/scale
      abmaxdiff=numpy.abs(amax-bmax)/scale
      #print amin,amax
      #print bmin,bmax

      if fldmaxdiff > 1e-7 and abmindiff >1e-5 and abmaxdiff > 1e-5:
         print "test_Abfilebathy_writeread_nomasl test passed"
      else :
         unittest.fail("AFile IO failed. MAx diff between read/written: %14.7g"%max([fldmaxdiff,abmaxdiff,abmindiff]))


   def test_abfilegrid_read(self) :
      regfile=modeltools.hycom.ABFileGrid("regional.grid","r")
      plon=regfile.read_field("plon")
      print regfile.fieldnames
      print plon
      #regfile.writefield(wfld,None)
      #regfile.close()

   def test_abfilegrid_write(self) :
      idm = random.randrange(10,5000)
      jdm = random.randrange(10,5000)
      scale=1e4
      print "Creating %dX%d bathy array"%(idm,jdm)
      wfldout=scale + numpy.random.rand(jdm,idm)*scale
      regfile=modeltools.hycom.ABFileGrid("test.grid","w")
      plon=regfile.write_field(wfldout,None,"plon")
      regfile.close()

   def test_abfilegrid_write(self) :
      idm = random.randrange(10,5000)
      jdm = random.randrange(10,5000)
      scale=1e4
      print "Creating %dX%d bathy array"%(idm,jdm)
      wfldout=scale + numpy.random.rand(jdm,idm)*scale
      regfile=modeltools.hycom.ABFileGrid("test.grid","w")
      plon=regfile.write_field(wfldout,None,"plon")
      regfile.close()


   def test_abfilegrid_writeread(self) :
      import modeltools.grid
      grid = modeltools.grid.Proj4Grid("+proj=stere +lat_ts=80 +lat_0=90 +lon_0=0",-45,60,20000,20000,200,300)
      modeltools.hycom.write_regional_grid(grid)
      inflds = modeltools.hycom.read_regional_grid()
      scale=1e-4

      flds={}
      flds["plon"],flds["plat"] = grid.pgrid()
      flds["ulon"],flds["ulat"] = grid.ugrid()
      flds["vlon"],flds["vlat"] = grid.vgrid()
      flds["qlon"],flds["qlat"] = grid.qgrid()
      flds["scpx"] = grid.scpx()
      flds["scux"] = grid.scux()
      flds["scvx"] = grid.scvx()
      flds["scqx"] = grid.scqx()
      flds["scpy"] = grid.scpy()
      flds["scuy"] = grid.scuy()
      flds["scvy"] = grid.scvy()
      flds["scqy"] = grid.scqy()
      flds["cori"] = grid.corio()
      flds["pang"] = grid.p_azimuth()
      flds["pasp"] = grid.aspect_ratio()

      for fldname in flds.keys() : 
         print fldname
         maxdiff=max_diff(flds[fldname],inflds[fldname])
         print fldname, maxdiff
         if maxdiff > 1e-7 :
            unittest.fail("AFile IO failed. MAx diff between %s read/written: %14.7g"%(fldname,maxdiff))


      #bmin,bmax = bathyfile.bminmax("depth")
      #amax = numpy.amax(wfldin)
      #amin = numpy.amin(wfldin)
      #abmindiff=numpy.abs(amin-bmin)/scale
      #abmaxdiff=numpy.abs(amax-bmax)/scale
      #print amin,amax
      #print bmin,bmax




      raise NameError,"test"


def max_diff(fld1,fld2) :
   absmax1=numpy.abs(numpy.amax(fld1))
   absmax2=numpy.abs(numpy.amax(fld2))
   absmax=numpy.max([absmax1,absmax2])
   absmaxdiff=numpy.abs(numpy.amax(fld1-fld2))
   return absmaxdiff/absmax
             
if __name__ == "__main__" :
   unittest.main()
