#!/usr/bin/env python
import modeltools.hycom
import argparse
import numpy
import logging


_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # Dont propagate to parent in hierarchy (determined by "." in __name__)


def micomsteps(blkdatfile,baclin_lo,baclin_hi) :
   # []=MICOMSTEPS(BACLIN,BATROP) displays valid time steps for the model MICOM.
   #    BACLIN and BATROP are two element vectors containing the lower and upper
   #    baroclinic and barotropic time steps, respectively.


   tmp = modeltools.hycom.BlkdatParser(blkdatfile)

   # Low frequency forcing
   if tmp["yrflag"] <= 1 :
      sd=86400;
   # High frequency forcing
   else :
      sd=6*3600;
   logger.info("Restriction: %d seconds must be divisible by baclin"%sd)


   # NB: Only integers
   bcdiv1 = int(numpy.floor(sd/baclin_lo))
   bcdiv2 = int(numpy.ceil(sd/baclin_hi))
   for i in range(bcdiv2,bcdiv1) :
      baclin = sd/float(i)
      logger.info("Allowed baroclinic time step %14.6f (%8d / %3d)"%(baclin,sd,i))

      # baclin_lo, baclin_hi not used
      for i in range(10,32,2) :
         batrop=baclin/float(i)
         logger.info("--->Allowed barotropic time step %14.6f (%8d / %3d)"%(batrop,baclin,i))


if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('blkdatfile', type=str, help='')
   parser.add_argument('baclin_lo', type=int, help='')
   parser.add_argument('baclin_hi', type=int, help='')

   args = parser.parse_args()
   micomsteps(
         args.blkdatfile,
         args.baclin_lo,
         args.baclin_hi
         )


