#!/usr/bin/env python
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
## AO from mpl_toolkits.basemap import Basemap
import modeltools.hycom
import logging
import re
import os
import os.path

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

# Note that this part can be easily implemented into the python code 
# that aims to generate nesting [ab] files. I may not have time to manage it 
# but evrything has been provided here or other earlier developed NERSC scripts.
# By Mostafa Bakhoday-Paskyabi
# July 4th 2019 (almost at last day of my work at NERSC)
#
# ../bin/calc_montg1.py ../nest/015/archv.2016_001_00.a ./
# ../bin/calc_montg1.py YOUR_INPUT_AFILE OUTPUT_PATH
#  It saves with same name as a-file to outpput path


cline1="HYCOM archive files from climatology generated from"
cline2="NEMO Relaxation fields on native grid, MBP2019"
cline3="Expt 01.5  nhybrd=50 nsigma= 0 ds00= 1.00 dp00= 1.00 dp00x= 450.0 dp00f=1.150"

def montg1_pb(thstar,p) :
   kdm=thstar.shape[0]
   thref=1e-3
   # Part that depends on pbavg :
   montgpb=-thstar[kdm-1,:,:]*thref**2
   for k in reversed(range(kdm-1)) :
      # PArt that depends on pbavg (through oneta) 
      montgpb[:,:]=montgpb[:,:]+ p[k+1,:,:]*(thstar[k+1,:,:]-thstar[k,:])*thref**2/p[kdm,:,:]
   return montgpb




# Get component of montgomery potential at surface that does not depend on pbavg
def montg1_no_pb(psikk,thkk,thstar,p) :
   kdm=thstar.shape[0]
   thref=1e-3
   montgc=psikk+(p[kdm,:,:]*(thkk-thstar[kdm-1,:,:]))*thref**2
   for k in reversed(range(kdm-1)) :
      montgc [:,:]=montgc [:,:]+ p[k+1,:,:]*(thstar[k+1,:,:]-thstar[k,:])*thref**2 
   return montgc



def main(archv_files,opath,header_line1=cline1,header_line2=cline2,header_line3=cline3 ) :


   # Open blkdat files. Get some properties
   bp=modeltools.hycom.BlkdatParser("blkdat.input")
   idm    = bp["idm"]
   jdm    = bp["jdm"]
   kdm    = bp["kdm"]
   thflag = bp["thflag"]
   thbase = numpy.float(bp["thbase"])
   kapref = bp["kapref"]
   iversn = bp["iversn"]
   iexpt  = bp["iexpt"]
   yrflag = bp["yrflag"]
   #thref=modeltools.hycom.thref
   thref=1e-3
   if kapref == -1 : 
      kapnum = 2
      msg="Only kapref>=0 is implemented for now"
      logger.error(msg)
      raise ValueError,msg
   else :
      kapnum = 1 

   if kapnum > 1 :
      msg="Only kapnum=1 is implemented for now"
      logger.error(msg)
      raise ValueError,msg


   # hycom sigma and kappa, written in python. NB: sigver is not used here.
   # Modify to use other equations of state. For now we assume sigver is:
   #    1 (7-term eqs referenced to    0 bar)
   #    2 (7-term eqs referenced to 2000 bar)
   if thflag == 0 :
      sigver=1
   else :
      sigver=2
   sig  = modeltools.hycom.Sigma(thflag)
   if kapref > 0  : kappa = modeltools.hycom.Kappa(kapref,thflag*1000.0e4) # 



   # Loop through archive files
   for archv_file in archv_files :
      archv_file=archv_file[:-2]
      logger.debug("Processing %s"%(archv_file))
      arcfile=abfile.ABFileArchv(archv_file,"r")
      
      temp    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      saln    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      th3d  =numpy.ma.zeros((kdm,jdm,idm))
      thstar=numpy.ma.zeros((kdm,jdm,idm))
      dp    =numpy.ma.zeros((jdm,idm))
      p     =numpy.ma.zeros((kdm+1,jdm,idm))
      pup  =numpy.ma.zeros((jdm,idm))
      montg=numpy.ma.zeros((jdm,idm))
      for k in range(kdm) :
         logger.info("Reading layer %d from %s"%(k,archv_file))
         temp  =arcfile.read_field("temp",k+1)
         saln  =arcfile.read_field("salin",k+1)
         dp    [:,:]=arcfile.read_field("thknss",k+1)
         th3d  [k  ,:,:]=sig.sig(temp,saln) - thbase
         p     [k+1,:,:]= p[k,:,:] + dp[:,:]
         if k>0 : thstarup   =numpy.ma.copy(thstar[k  ,:,:])
         thstar[k  ,:,:]=numpy.ma.copy(th3d  [k  ,:,:])
         if kapref > 0 :
            thstar[k  ,:,:]=thstar  [k  ,:,:] + kappa.kappaf(temp[:,:],
                                                           saln[:,:], 
                                                           th3d[k,:,:]+thbase,
                                                           p[k,:,:])
         elif kapref < 0 :
            msg="Only kapref>=0 is implemented for now"
            logger.error(msg)
            raise ValueError,msg
         if k > 0 : 
            montg = montg - p[k,:,:] * (thstar[k  ,:,:] - thstarup) * thref**2

      thkk = thstar[k  ,:,:]
      psikk = montg

      # This part of montg1 does nto depend on pbavg
      montg1c  = montg1_pb(thstar,p)
      # This part depends linearly on pbavg
      montg1pb = montg1_no_pb(psikk,thkk,thstar,p) 
      # Barotropic pressure. NB: srfhgt is in geopotential height : 
      srfhgt=arcfile.read_field("srfhgt",0)
      pbavg  = (srfhgt - montg1c) / (montg1pb+thref)
      montg1 = montg1pb*pbavg + montg1c
      logger.info("Estimated montg1 ")
 
      # Apply mask
      #montg1=numpy.ma.masked_where(srfhgt.mask,montg1) 
      #montg1[~msk]=2**100
      #montg1 = numpy.ma.masked_where(srfhgt<2.0**99,montg1)
      msk=numpy.where(srfhgt>2.0**99,1,0)


      # Open new archive file, and write montg1 to it.
      fnameout = opath+str(archv_file[-17:])
      arcfile_out=abfile.ABFileArchv(fnameout,"w",iversn=iversn,yrflag=yrflag,iexpt=iexpt,mask=False,cline1=header_line1,cline2=header_line2,cline3=header_line3)

      fields=arcfile.get_fields()
      for key in sorted(fields.keys()) :
         fieldname = fields[key]["field"]
         time_step = fields[key]["step"]
         model_day = fields[key]["day"]
         k         = fields[key]["k"]
         dens      = fields[key]["dens"]
         fld       =arcfile.read_field(fieldname,k)

         if fieldname == "montg1" :
            logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,fnameout))
            arcfile_out.write_field(montg1,None,fieldname,time_step,model_day,k,dens) 
         else :
            arcfile_out.write_field(fld   ,None,fieldname,time_step,model_day,k,dens) 

      logger.info("Finished writing to %s"%fnameout)
      arcfile_out.close()
      arcfile.close()



if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='Calculate montgomery potential in 1st layer (used in archv files)')
   parser.add_argument('--header_line1',          type=str, help="First header line, applicable only for archm & archv",default=cline1)
   parser.add_argument('--header_line2',          type=str, help="Montgomery with or without barotropic pressure",default=cline2)
   parser.add_argument('--header_line3',          type=str, help="Montgomery with or without barotropic pressure",default=cline3)

   parser.add_argument('archive_file',     type=str, help='Files to add montg1 to',nargs="+")
   parser.add_argument('opath',     type=str, help='Output files path')
   args = parser.parse_args()
   main(args.archive_file,args.opath,header_line1=args.header_line1,header_line2=args.header_line2,header_line3=args.header_line3)
