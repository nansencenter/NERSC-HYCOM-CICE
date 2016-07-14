#!/usr/bin/env python
import modelgrid
import numpy
import abfile

def main() :
   print "ok"
   cg = modelgrid.ConformalMapping(
      -40.,140.,
      -50.,140.,
      178.2,181.42,800,
      0.,80,760,
      True,
      0.365,False)

   la_n,lo_n = cg.oldtonew(70,-45)
   #print la_n,lo_n
   ip,jp =  cg.pivotp(la_n,lo_n)
   #print ip,jp
   #print numpy.mod(-180.,360.)
   i,j= cg.ll2gind(70,-45)
   #print i,j
   lo,la = cg.gind2ll(i,j)
   #print lo,la


   tmp =  modelgrid.ConformalGrid(
      -40.,140.,
      -50.,140.,
      178.2,181.42,1000,
      0.,80,950,
      True,
      0.365,False)
   plat,plon=tmp.pgrid()
   #print plat.max()
   #print tmp.scpx()
   print plat.shape

   fig=tmp.plotgrid(1.5)
   fig.savefig("tst.png")
   tmp.save_to_scrip("tst.nc")

   datadict=tmp.create_datadict_hycom()
   abfile.write_regional_grid(datadict,endian="big")


if __name__ == "__main__":
   print "hei"
   main()
