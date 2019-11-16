#!/usr/bin/env python
import unittest
#import timetools
import random
import numpy
import abfile

class ABFileTest(unittest.TestCase) :

   def test_afile_writeread_nomask(self) :

      idm = random.randrange(10,5000)
      jdm = random.randrange(10,5000)
      #print "Creating %dX%d random arrays"%(idm,jdm)
      scale=1e4
      wfld1=numpy.random.rand(jdm,idm)*scale
      wfld2=numpy.random.rand(jdm,idm)*scale
      wfld3=numpy.random.rand(jdm,idm)*scale

      afile=abfile.AFile(idm,jdm,"test.a","w")
      afile.writerecord(wfld1,None,record=None)
      afile.writerecord(wfld2,None,record=None)
      afile.writerecord(wfld3,None,record=None)
      afile.close()

      afile=abfile.AFile(idm,jdm,"test.a","r")
      rfld1=afile.read_record(record=0)
      rfld2=afile.read_record(record=1)
      rfld3=afile.read_record(record=2)
      afile.close()

      maxdiff1=numpy.abs(numpy.amax(rfld1-wfld1))/scale
      maxdiff2=numpy.abs(numpy.amax(rfld2-wfld2))/scale
      maxdiff3=numpy.abs(numpy.amax(rfld3-wfld3))/scale
      #print maxdiff1,maxdiff2,maxdiff3

      if all([elem < 1e-7 for elem in [maxdiff1,maxdiff2,maxdiff3]]) :
         #print "afile io test passed for %d x %d array " %(idm,jdm)
         pass
      else :
         self.fail("AFile IO failed. MAx diff between read/written: %14.7g"%max([maxdiff1,maxdiff2,maxdiff3]))

   def test_abfilebathy_writeread_nomask(self) :
      idm = random.randrange(10,5000)
      jdm = random.randrange(10,5000)
      scale=1e4

      #print "Creating %dX%d bathy array"%(idm,jdm)
      wfldout=scale + numpy.random.rand(jdm,idm)*scale
      bathyfile=abfile.ABFileBathy("testbathy","w")
      bathyfile.write_field(wfldout,None)
      bathyfile.close()

      #print "Reading %dX%d bathy array"%(idm,jdm)
      bathyfile=abfile.ABFileBathy("testbathy","r",idm=idm,jdm=jdm)
      wfldin=bathyfile.read_field("depth",None)
      bathyfile.close()

      fldmaxdiff=numpy.abs(numpy.amax(wfldin-wfldout))/scale
      bmin,bmax = bathyfile.bminmax("depth")
      amax = numpy.amax(wfldin)
      amin = numpy.amin(wfldin)
      abmindiff=numpy.abs(amin-bmin)/scale
      abmaxdiff=numpy.abs(amax-bmax)/scale
      if not (fldmaxdiff > 1e-7 and abmindiff >1e-5 and abmaxdiff > 1e-5):
         pass
      else :
         self.fail("AFile IO failed. MAx diff between read/written: %14.7g"%max([fldmaxdiff,abmaxdiff,abmindiff]))


#   def test_abfilegrid_read(self) :
#      regfile=abfile.ABFileGrid("regional.grid","r")
#      plon=regfile.read_field("plon")

   def test_abfilegrid_write(self) :
      idm = random.randrange(10,5000)
      jdm = random.randrange(10,5000)
      scale=1e4
      #print "Creating %dX%d bathy array"%(idm,jdm)
      wfldout=scale + numpy.random.rand(jdm,idm)*scale
      regfile=abfile.ABFileGrid("test.grid","w")
      plon=regfile.write_field(wfldout,None,"plon")
      regfile.close()


#   def test_abfilegrid_contents(self) :
#      regfile=abfile.ABFileGrid("regional.grid","r")
#      for k 
      

   def test_abfilearchv_read(self) :
      archv=abfile.ABFileArchv("archv.2013_152_00","r")
      fld=archv.read_field("salin",1)
      #print fld.min(),fld.max()
      #print archv.fieldnames
      #print archv.fieldlevels
      archv.close()


   def test_abfileforcing_write(self) :
      archv=abfile.ABFileForcing("test.airtmp","w",cline1="test",cline2="test")
      fld=numpy.random.rand(100,100)
      archv.write_field(fld,fld,"airtmp",39814.0000,0.25)
      archv.close()

   def test_abfileforcing_writeread(self) :
      scale=1e2
      archv=abfile.ABFileForcing("test.airtmp","w",cline1="test",cline2="test")
      fld1=numpy.random.rand(100,100)
      fld2=numpy.random.rand(100,100)
      fld3=numpy.random.rand(100,100)
      archv.write_field(fld1,fld1,"airtmp",39814.0000,0.25)
      archv.write_field(fld2,fld2,"airtmp",39814.5000,0.25)
      archv.write_field(fld3,fld3,"airtmp",39815.0000,0.25)
      archv.close()

      archv=abfile.ABFileForcing("test.airtmp","r")
      fldin = archv.read_field("airtmp",39814.5000)
      bmin,bmax = archv.bminmax("airtmp",39814.5)

      fldmaxdiff=numpy.abs(numpy.amax(fldin-fld2))/scale
      amax = numpy.amax(fldin)
      amin = numpy.amin(fldin)
      abmindiff=numpy.abs(amin-bmin)/scale
      abmaxdiff=numpy.abs(amax-bmax)/scale
      #print abmaxdiff
      if not (fldmaxdiff > 1e-7 and abmindiff >1e-5 and abmaxdiff > 1e-5):
         pass
      else :
         self.fail("AFile IO failed. MAx diff between read/written: %14.7g"%max([fldmaxdiff,abmaxdiff,abmindiff]))


def max_diff(fld1,fld2) :
   absmax1=numpy.abs(numpy.amax(fld1))
   absmax2=numpy.abs(numpy.amax(fld2))
   absmax=numpy.max([absmax1,absmax2])
   absmaxdiff=numpy.abs(numpy.amax(fld1-fld2))
   return absmaxdiff/absmax
             
if __name__ == "__main__" :
   unittest.main(verbosity=2)
