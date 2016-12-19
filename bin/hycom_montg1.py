#!/usr/bin/env python
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import abfile
import numpy
from mpl_toolkits.basemap import Basemap
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


# Northern line 
northern_limit_longitudes = [-35 ,-16, -6.6, -4.2, -0.3,   4., 31.6, 97., 172.,-172,
   -150., -129, -101.,-88, -83.84, -82.7, -75., -46.]
northern_limit_latitudes  = [67.5, 65,  62.,  57.7, 51.0, 49.4, 51.5, 57., 66.7, 66.1,
   65.,    64.,  64., 66.6, 67.81,  73.7, 81.9, 80.5]


def main(psikk_file,archv_files, psikk_file_type="restart",month=1) :

   logger.info("psikk file type is %s"%psikk_file_type)
   from_relax_archv = psikk_file_type=="relax_archv"
   from_relax       = psikk_file_type=="relax"
   from_restart     = psikk_file_type=="restart"
   if not from_relax_archv and not from_relax and not from_restart :
      msg="psikk_file_type must be either restart relax or relax_archv"
      logger.error(msg)
      raise ValueError,msg

   # Open blkdat files. Get some properties
   bp=modeltools.hycom.BlkdatParser("blkdat.input")
   idm    = bp["idm"]
   jdm    = bp["jdm"]
   kdm    = bp["kdm"]
   thflag = bp["thflag"]
   thbase = bp["thbase"]
   kapref = bp["kapref"]
   iversn = bp["iversn"]
   iexpt  = bp["iexpt"]
   yrflag = bp["yrflag"]
   thref=modeltools.hycom.thref
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

   path0=os.path.join(".","montg1")
   if os.path.exists(path0) and os.path.isdir(path0) :
      pass
   else : 
      os.mkdir(path0)

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

   # Option 1: Get psikk and thkk directly from a restart file
   if from_restart :

      # Read input variables from lowest layer of a restart file
      m=re.match( "^(.*)(\.[ab])", psikk_file)
      if m : psikk_file=m.group(1)
      rfile=abfile.ABFileRestart(psikk_file,"r",idm=idm,jdm=jdm)
      psikk=rfile.read_field("psikk",0)
      thkk=rfile.read_field("thkk",0)
      rfile.close()

   # Option 2: Get temp, salinity and dp from relaxaiton fields. Estimate thkk and psikk
   elif from_relax :
      pattern="^(.*)_(tem|int|sal)"
      psikk_file = abfile.ABFile.strip_ab_ending(psikk_file)
      m=re.match(pattern,psikk_file)
      if m :  
         relaxtem=abfile.ABFileRelax(m.group(1)+"_tem","r")
         relaxsal=abfile.ABFileRelax(m.group(1)+"_sal","r")
         relaxint=abfile.ABFileRelax(m.group(1)+"_int","r")
      else :
         msg="""Input hycom relaxation file %s does not match pattern
         %s"""%(psikk_file,pattern)
         logger.error(msg)
         raise ValueError,msg

   # Option 3: Get temp, salinity and dp from relaxation fields (archive files generated during relaxation setup). Estimate thkk and psikk
   elif from_relax_archv :
      arcfile0=abfile.ABFileArchv(psikk_file,"r")

   else :
      msg="No way of estimating psikk and thkk. Aborting ..."
      logger(msg)
      raise ValueError, msg

   # Estimate psikk, thkk from climatology. (Option 1 or option 2)
   if from_relax or from_relax_archv :
      logger.info("Estimating psikk and thkk from climatology")
      p    =numpy.ma.zeros((jdm,idm))
      pup  =numpy.ma.zeros((jdm,idm))
      montg=numpy.ma.zeros((jdm,idm))
      for k in range(kdm) :
         logger.debug("Reading layer %d from climatology"%k)
         if k>0 : thstarup   =numpy.ma.copy(thstar)

         if from_relax :
            saln       =relaxsal.read_field("sal",k+1,month)
            temp       =relaxtem.read_field("tem",k+1,month)
            plo        =relaxint.read_field("int",k+1,month) # lowest interface of layer
            dp         = plo - pup
            pup        = numpy.copy(plo)              # Next loop, pup = plo of this loop
         elif from_relax_archv :
            saln       =arcfile0.read_field("salin",k+1)
            temp       =arcfile0.read_field("temp",k+1)
            dp         =arcfile0.read_field("thknss",k+1)
         else :
            msg="No way of estimating sal, tem and dp. Aborting ..."
            logger(msg)
            raise ValueError, msg

         th3d       =sig.sig(temp,saln) - thbase
         thstar     =numpy.ma.copy(th3d)
         if kapref > 0 :
            thstar=thstar + kappa.kappaf(temp[:,:],
                                         saln[:,:], 
                                         th3d[:,:]+thbase,
                                         p[:,:])
         elif kapref < 0 :
            msg="Only kapref>=0 is implemented for now"
            logger.error(msg)
            raise ValueError,msg
         # From hycom inicon.f
         #        montg(i,j,1,kkap)=0.0
         #        do k=1,kk-1
         #          montg(i,j,k+1,kkap)=montg(i,j,k,kkap)-
         #     &    p(i,j,k+1)*(thstar(i,j,k+1,kkap)-thstar(i,j,k,kkap))*thref**2
         #        enddo
         #c
         #        thkk( i,j,kkap)=thstar(i,j,kk,kkap)
         #        psikk(i,j,kkap)=montg( i,j,kk,kkap)
         if k > 0 : 
            #print (thstar - thstarup).min(), (thstar - thstarup).min()
            montg = montg - p * (thstar - thstarup) * thref**2
         p   = p +dp
      thkk = thstar
      psikk = montg
      if from_relax :
         relaxtem.close()
         relaxsal.close()
         relaxint.close()
      elif from_relax_archv :
         arcfile0.close()
   else :
      pass
   logger.info("Min max of thkk: %12.6g %12.6g"%(thkk.min(),thkk.max()))
   logger.info("Min max of psikk: %12.6g %12.6g"%(psikk.min(),psikk.max()))






   # Loop through archive files
   for archv_file in archv_files :

      logger.debug("Processing %s"%(archv_file))
      arcfile=abfile.ABFileArchv(archv_file,"r")

      # Read all layers .. (TODO: If there is memory problems, read and estimate sequentially)
      temp    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      saln    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      th3d  =numpy.ma.zeros((kdm,jdm,idm))
      thstar=numpy.ma.zeros((kdm,jdm,idm))
      dp    =numpy.ma.zeros((jdm,idm))
      p     =numpy.ma.zeros((kdm+1,jdm,idm))
      for k in range(kdm) :
         logger.info("Reading layer %d from %s"%(k,archv_file))
         temp  =arcfile.read_field("temp",k+1)
         saln  =arcfile.read_field("salin",k+1)
         #dp    [k  ,:,:]=arcfile.read_field("thknss",k+1)
         dp    [:,:]=arcfile.read_field("thknss",k+1)
         th3d  [k  ,:,:]=sig.sig(temp,saln) - thbase
         p     [k+1,:,:]= p[k,:,:] + dp[:,:]
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

      # This part of montg1 does nto depend on pbavg
      montg1c  = modeltools.hycom.montg1_pb(thstar,p) 

      # This part depends linearly on pbavg
      montg1pb = modeltools.hycom.montg1_no_pb(psikk,thkk,thstar,p) 
      print montg1c.min(),montg1c.max()
      print montg1pb.min(),montg1pb.max()


# ... we have ...
#     montg1 = montgc + montgpb * pbavg
#     srfhgt = montg1 + thref*pbavg = montgc + montgpb * pbavg + thref*pbavg
#  ... which gives new montg1n for srfhgtn
#     pbavgn  = (srfhgtn-montgc) / (montgpb+thref)
#     montg1n = montgc + montgpb*pbavgn
# Systematic differences due to choice of sigma and or eq of state is not taken into account ...


      # Barotropic pressure. NB: srfhgt is in geopotential height : 
      srfhgt=arcfile.read_field("srfhgt",0)
      pbavg  = (srfhgt - montg1c) / (montg1pb+thref)
      print pbavg.min(),pbavg.max()
      montg1 = montg1pb*pbavg + montg1c
      logger.info("Estimated montg1 ")
      print montg1.min(),montg1.max()

      # Open new archive file, and write montg1 to it.
      fnameout = os.path.join(path0,os.path.basename(arcfile.basename))
      arcfile_out=abfile.ABFileArchv(fnameout,"w",iversn=iversn,yrflag=yrflag,iexpt=iexpt,mask=False)



      for key in sorted(arcfile.fields.keys()) :
         fieldname = arcfile.fields[key]["field"]
         time_step = arcfile.fields[key]["step"]
         model_day = arcfile.fields[key]["day"]
         k         = arcfile.fields[key]["k"]
         dens      = arcfile.fields[key]["dens"]
         fld       =arcfile.read_field(fieldname,k)

         if fieldname == "montg1" :
            logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,fnameout))
            arcfile_out.write_field(montg1,None,fieldname,time_step,model_day,sigver,thbase) 
         else :
            arcfile_out.write_field(fld   ,None,fieldname,time_step,model_day,k,dens) 
            #logger.info("Writing field %10s at level %3d to %s (copy from original)"%(fieldname,k,fnameout))

      logger.info("Finished writing to %s"%fnameout)
      arcfile_out.close()
      arcfile.close()

   logger.warning("Sigver assumed to be those of 7 term eqs")
   logger.warning("    1 for sigma-0/thflag=0, 2 for sigma-2/thflag=2")
   logger.warning("For other eqs, you need to modify the code so that sigver is set correctly in the archv file""") 
   logger.warning("psikk and thkk from relaxation fields is beta. Preferred method is restart")

   logger.info("Finito")

if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='Calculate montgomery potential in 1st layer (used in archv files)')
   parser.add_argument('--month',          type=int, help="Month to process if psikkfile is relaxation fields. Defaults to 1",default=1)
   parser.add_argument('--psikk-file-type',       type=str, help="""File to get psikk and
         thkk from. Accepted values are restart (a restart file containing psikk
         and thkk, relax_archv (archive files created when running reanalysis),
         relax (relaxation files)""", default="restart") 
   parser.add_argument('psikk_file',       type=str, help='file to get thkk and psikk from ')
   parser.add_argument('archive_file',     type=str, help='Files to add montg1 to',nargs="+")
   args = parser.parse_args()
   main(args.psikk_file,args.archive_file,month=args.month,psikk_file_type=args.psikk_file_type)
