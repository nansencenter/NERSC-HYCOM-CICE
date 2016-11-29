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
_loglevel=logging.DEBUG
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


def main(restart_file,archv_files,sigver) :


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

   path0=os.path.join(".","montg1")
   if os.path.exists(path0) and os.path.isdir(path0) :
      pass
   else :
      msg="This routine assumes there is a directory called %s in the current directory. Create that first"%path0
      logger.error(msg)
      raise NameError,msg
         

   # Initialize relevant sigma class



   ## Read input bathymetry
   #bfile=abfile.ABFileBathy("regional.depth","r",idm=gfile.idm,jdm=gfile.jdm,mask=True)
   #in_depth_m=bfile.read_field("depth")
   #bfile.close()

   # Read input variables from lowest layer
   m=re.match( "^(.*)(\.[ab])", restart_file)
   if m : restart_file=m.group(1)
   rfile=abfile.ABFileRestart(restart_file,"r",idm=idm,jdm=jdm)
   psikk=rfile.read_field("psikk",0)
   thkk=rfile.read_field("thkk",0)
   rfile.close()

   # hycom sigma and kappa, written in python
   sig  = modeltools.hycom.Sigma(thflag)
   if kapref > 0  : kappa = modeltools.hycom.Kappa(kapref,thflag*1000.0e4) # 



   # Loop through archive files
   for archv_file in archv_files :

      logger.info("Processing %s"%(archv_file))
      arcfile=abfile.ABFileArchv(archv_file,"r")

      # Read all layers .. (TODO: If there is memory problems, read and estimate sequentially)
      temp    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      saln    = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      montgpb = numpy.ma.zeros((jdm,idm))    # Only need top layer
      montgc  = numpy.ma.zeros((jdm,idm))    # Only need top layer
      th3d  =numpy.ma.zeros((kdm,jdm,idm))
      thstar=numpy.ma.zeros((kdm,jdm,idm))
      dp    =numpy.ma.zeros((kdm,jdm,idm))
      p     =numpy.ma.zeros((kdm+1,jdm,idm))
      for k in range(kdm) :
         temp  =arcfile.read_field("temp",k+1)
         saln  =arcfile.read_field("salin",k+1)
         dp    [k  ,:,:]=arcfile.read_field("thknss",k+1)
         th3d  [k  ,:,:]=sig.sig(temp,saln) - thbase
         p     [k+1,:,:]= p[k,:,:] + dp[k,:,:]
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


#     Original fortran code (mod_momtum)
#     # m_prime in lowest layer
#     montg(i,j,kk,1)=psikk(i,j,1)+
#    &        ( p(i,j,kk+1)*(thkk(i,j,1)-thstar(i,j,kk,1))
#    &          -pbavg(i,j,m)*thstar(i,j,kk,1) )*thref**2
#     if     (kapref.eq.-1) then
#             montg(i,j,kk,2)=psikk(i,j,2)+
#    &          ( p(i,j,kk+1)*(thkk(i,j,2)-thstar(i,j,kk,2))
#    &            -pbavg(i,j,m)*thstar(i,j,kk,2) )*thref**2
#           endif !kapref.eq.-1
      if kapref >= 0 :

         # Part that depends on pbavg :
         montgpb[:,:]=-thstar[kdm-1,:,:]*thref**2

         # The rest
         montgc[:,:]=psikk+(p[kdm,:,:]*(thkk-thstar[kdm-1,:,:]))*thref**2

      else :
         msg="Only kapref>=0 is implemented for now"
         logger.error(msg)
         raise ValueError,msg

      #print "montgc:",montgc[itest,jtest],montgc[itest,jtest],thstar[kdm-1,itest,jtest],thkk[itest,jtest]
      #print "montgpb:",montgpb[itest,jtest],montgpb[itest,jtest],thkk[itest,jtest],psikk[itest,jtest]
      #raise ValueError,"test"


#     Original fortran code (mod_momtum)
#c ---       m_prime in remaining layers:
#            do k=kk-1,1,-1
#              montg(i,j,k,1)=montg(i,j,k+1,1)+p(i,j,k+1)*oneta(i,j)
#     &            *(thstar(i,j,k+1,1)-thstar(i,j,k,1))*thref**2
#              if     (kapref.eq.-1) then
#                montg(i,j,k,2)=montg(i,j,k+1,2)+p(i,j,k+1)*oneta(i,j)
#     &              *(thstar(i,j,k+1,2)-thstar(i,j,k,2))*thref**2
#              endif !kapref.eq.-1
#            enddo !k



      if kapref >= 0 :
         for k in reversed(range(kdm-1)) :
            #print "montgc:",montgc[itest,jtest],montgc[itest,jtest],(thstar[k+1,itest,jtest]-thstar[k,itest,jtest])
            #print "montgpb:",montgpb[itest,jtest],montgpb[itest,jtest]
            #print k,kdm
            # PArt that depends on pbavg (through oneta) 
            montgpb[:,:]=montgpb[:,:]+ p[k+1,:,:]*(thstar[k+1,:,:]-thstar[k,:])*thref**2/p[kdm,:,:]

            # The rest
            montgc [:,:]=montgc [:,:]+ p[k+1,:,:]*(thstar[k+1,:,:]-thstar[k,:])*thref**2 
      else :
         msg="Only kapref>=0 is implemented for now"
         logger.error(msg)
         raise ValueError,msg

# ... we have ...
#     montg1 = montgc + montgpb * pbavg
#     srfhgt = montg1 + thref*pbavg
#  ... which gives ...
#     pbavg  =(srfhgt - montgc) /(montgpb + thref)
#     srfhgt = montg1 + thref*pbavg
# Systematic differences due to choice of sigma and or eq of state is not taken into account ...


      # Barotropic pressure. NB: srfhgt is in geopotential height : 
      srfhgt=arcfile.read_field("srfhgt",0)
      pbavg  = (srfhgt - montgc[:,:]) / (montgpb[:,:]+thref)
      montg1 = montgpb[:,:]*pbavg + montgc[:,:]
      #print "srfhgt",srfhgt[itest,jtest]
      #print "pbavg",pbavg[itest,jtest]
      #print "montg1",montg1[itest,jtest]
      logger.info("Estimated montg1 ")

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

   logger.info("Finito")

if __name__ == "__main__" :

   parser = argparse.ArgumentParser(description='Calculate montgomery potential in 1st layer (used in archv files)')
   parser.add_argument('sigver',           type=int, help='Version of equation of state')
   parser.add_argument('restart_file',     type=str, help='Restart file to get thkk and psikk from ')
   parser.add_argument('archive_file',     type=str, help='Files to add montg1 to',nargs="+")
   args = parser.parse_args()
   main(args.restart_file,args.archive_file,args.sigver)
