"""
Replace NaN values to Zeros for the Montgomery potential AO
"""
#!/usr/bin/env python
import argparse
import abfile
import numpy as np
import modeltools.hycom
import logging
import re
import os
import os.path
import pickle
#############
# this code takes the nesting archive files that were already interpolated
# both vertically and horizontally (these files are  possibly located in nest folder),
# converts ssh (anomaly) to montg1 variable using the slope and intercept from 
# regress_file. Use the scripts as (in experiment folder):
#
# source ../REGION.src
# python ../bin/fix_montg1.py ../nest/011/archv.2016_001_00.a ${INPUTDIR}montg_regress.pckl ./  ../topo/regional.depth.a
#
# this code works standalone, so you can use it on existing archive nest files
# no need to do this if you're creating a new archive file since the code is already embedded to bin/archvz2hycom_biophys.sh
#############
# Set up logger
_loglevel=logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False

cline1="HYCOM nested archive files"
cline2="NEMO fields on hycom grid with corrected montg1 variable"
cline3="Expt 01.1  nhybrd=50 nsigma= 0 ds00= 1.00 dp00= 1.00 dp00x= 450.0 dp00f=1.150"


def main(archv_files,regress_file,opath,depth_file,header_line1=cline1,header_line2=cline2,header_line3=cline3 ) :

   regr = open(regress_file,'rb')
   slope,intercept,nemomeandt = pickle.load(regr,encoding="latin1")

   bp=modeltools.hycom.BlkdatParser("blkdat.input")
   idm    = bp["idm"]
   jdm    = bp["jdm"]
   kdm    = bp["kdm"]
   iversn = bp["iversn"]
   yrflag = bp["yrflag"]
   iexpt  = bp["iexpt"]

   # get domain depth
   abdepth = abfile.ABFileBathy(depth_file, \
            "r",idm=idm,jdm=jdm)
   depthm=abdepth.read_field("depth")   

   for archv_file in archv_files :

      logger.debug("Processing %s"%(archv_file))
      arcfile=abfile.ABFileArchv(archv_file,"r")
      
      srfhgt=arcfile.read_field("srfhgt",0)
      montg1=arcfile.read_field("montg1",0)
      montg1=(srfhgt-nemomeandt) * slope + intercept
      montg1=np.where(np.isnan(montg1), 0, montg1)
      dummy = arcfile.read_field("temp",1)
      dummy[~depthm.mask] = montg1[~depthm.mask] # your fixed montg1
      logger.info("Montg Minimum")
      print(dummy.min())
      logger.info("Montg Maximum")
      print(dummy.max())
#      logger.info("Estimated montg1 ")
      fnameout = opath+str(archv_file[-19:])
      arcfile_out=abfile.ABFileArchv(fnameout,"w",iversn=iversn,yrflag=yrflag,iexpt=iexpt,mask=False,cline1=header_line1,cline2=header_line2,cline3=header_line3)

      for keys in sorted( arcfile.fields.keys() ) :
          fieldname = arcfile.fields[keys]["field"]
          time_step = arcfile.fields[keys]["step"]
          model_day = arcfile.fields[keys]["day"]
          k         = arcfile.fields[keys]["k"]
          dens      = arcfile.fields[keys]["dens"]
          fld       = arcfile.read_field(fieldname,k)

          if fieldname == "montg1" :
            logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,fnameout))
            arcfile_out.write_field(dummy,None,fieldname,time_step,model_day,k,dens)
          else :
            arcfile_out.write_field(fld   ,None,fieldname,time_step,model_day,k,dens)

      logger.info("Finished writing to %s"%fnameout)
      arcfile_out.close()
      arcfile.close()

      archv_file_b = archv_file[0:-1] + 'b'
      fnameout_b = fnameout[0:-1] + 'b'

      os.rename(fnameout,archv_file) # replaces the new archive files with the original ones in nest folder
      os.rename(fnameout_b,archv_file_b)

if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='Calculate montgomery potential in 1st layer (used in archv files)')
   parser.add_argument('--header_line1',          type=str, help="First header line, applicable only for archm & archv",default=cline1)
   parser.add_argument('--header_line2',          type=str, help="Montgomery with or without barotropic pressure",default=cline2)
   parser.add_argument('--header_line3',          type=str, help="Montgomery with or without barotropic pressure",default=cline3)

   parser.add_argument('archive_file',     type=str, help='Files to add montg1 to',nargs="+")
   parser.add_argument('opath',     type=str, help='Output files path')
   parser.add_argument('regress_file', type=str, help='regression file path')
   parser.add_argument('depth_file', type=str, help='path to regional_depth.a file')
   args = parser.parse_args()
   main(args.archive_file,args.opath,args.regress_file,args.depth_file,header_line1=args.header_line1,header_line2=args.header_line2,header_line3=args.header_line3)

