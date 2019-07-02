#!/usr/bin/env python
import argparse
import numpy

def micomsteps(batrop_lo,batrop_hi,baclin_lo,baclin_hi) :
   # []=MICOMSTEPS(BACLIN,BATROP) displays valid time steps for the model MICOM.
   #    BACLIN and BATROP are two element vectors containing the lower and upper
   #    baroclinic and barotropic time steps, respectively.
   sd=86400;
   bc=numpy.arange(baclin_lo,baclin_hi+1,1)
   bc=bc[sd%bc == 0]
   for n in bc :
      bt=numpy.arange(batrop_lo,batrop_hi+1,1)
      bt=bt[(n/2)%bt==0]
      print "Allowed baroclinic time step %d"%n
      for m in bt :
         print "**Allowed barotropic time step %4d for baroclinic time step %d (ratio %d)"%(n,m,n/m)


if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('batrop_lo', type=int, help='')
   parser.add_argument('batrop_hi', type=int, help='')
   parser.add_argument('baclin_lo', type=int, help='')
   parser.add_argument('baclin_hi', type=int, help='')

   args = parser.parse_args()
   micomsteps(
         args.batrop_lo,
         args.batrop_hi,
         args.baclin_lo,
         args.baclin_hi
         )

