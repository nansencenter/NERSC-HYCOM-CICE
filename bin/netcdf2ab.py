#!/usr/bin/env python
import netCDF4
import numpy as np
import numpy
import logging
import netcdftime
import argparse
import abfile


# This script convert netcdf file into [ab] files for the purpose of creating
# offlux files. 
#
# Inputs
# filename : input netcdf file name
# fieldname: field which will be saved in [ab] file
# OPTIONAL
# maskfield: It possible you decide to use for example depth as masking file
#            you then need to specify the mask field's name here.
#
# Mostafa Bakhoday-Paskyabi: 4 Feb. 2019
#
#

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # TODO: :Not sure why this is needed...




if __name__ == "__main__":
   class ClimParseAction(argparse.Action) :
      def __call__(self, parser, args, values, option_string=None):
         tmp = values.split(",")
         tmp = [float(elem) for elem in tmp[0:2]]
         setattr(args, self.dest, tmp)
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('filename', help="",nargs='+')
   parser.add_argument('fieldname', type=str, help="")
   parser.add_argument('--maskfield',default=None, type=str, help="")
   parser.add_argument('--ofieldname',default=None, type=str, help="")
   args = parser.parse_args()


   for myfile0 in args.filename :
      logger.info("Opening netcdf file %s"%myfile0)
      ncid=netCDF4.Dataset(myfile0,"r")
      ncid.set_auto_mask(False) 
      field = ncid.variables[args.fieldname][:,:].squeeze()
      field=np.where(np.abs(field)>1e+2,0,field)
      print field.min(),field.max()
      if args.maskfield is not None:
         dummy = ncid.variables[args.maskfield][:,:]
         mask=np.where(np.abs(dummy)>1e+5,0,1) # here depth
      else:
         print "Field data will be used to diagnose dry/wet cells."
         mask=np.where(np.abs(field)>1e+2,0,1)
      jdm,idm=field.shape
      print idm,jdm
      outfile=abfile.ABFileForcing(myfile0[:-3],"w",idm=idm, jdm=jdm)
      if args.ofieldname is not None: 
         outfile.write_field(field,mask,args.ofieldname,0,0)
      else:
         outfile.write_field(field,mask,"sst",0,0)


