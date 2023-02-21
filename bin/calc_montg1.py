#!/usr/bin/env python
import modeltools.hycom
import subprocess
import argparse
import datetime
import os.path
import sys
import re
import logging
import abfile
import pickle
import numpy
# Set up logger
##_loglevel=logging.DEBUG
_loglevel=logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False


"""
usage:
python ../bin/topaz6_calc_montg1.py archv_file restart_file  output_path
python ../bin/topaz6_calc_montg1.py ../nest/070/archv.2022_004_00.a ./data/TP6restart.2021_365_00_0000.a  ./

python /cluster/home/achoth/pyScripts_AO_Betzy/ValidationFreq/Alf_calc_montg1.py /cluster/work/users/achoth/TP5a0.06/nest/080_NewMontg/archv.1997_001_00.b /cluster/work/users/achoth/TP5a0.06/expt_08.1/data/restart.1997_001_00_0000.b /cluster/work/users/achoth/TP5a0.06/nest/test_dir/


"""
# Get component of montgomery potential at surface that depends on pbavg
def montg1_pb(thstar,p) :
   kdm=thstar.shape[0]
   #AO thref=_constants.thref
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
   ##AO thref=_constants.thref
   thref=1e-3
   montgc=psikk+(p[kdm,:,:]*(thkk-thstar[kdm-1,:,:]))*thref**2
   for k in reversed(range(kdm-1)) :
      montgc [:,:]=montgc [:,:]+ p[k+1,:,:]*(thstar[k+1,:,:]-thstar[k,:])*thref**2
   return montgc


def approx_montg1(thstar,p,srfhgt,idm,jdm,rstr_file):
         # estimate montg1---------------------------------------------
         # Read input variables from lowest layer of a restart file
         thref=1e-3
         #psikk_file="./TP6restart.2021_041_00_0000.a"
         psikk_file=rstr_file#[0]#[-30:]
         logger.info("----->rstr_file :%s"%psikk_file)
         m=re.match( "^(.*)(\.[ab])", psikk_file)
         if m : psikk_file=m.group(1)
         logger.info("----->rstr_file :%s"%psikk_file)
         rfile=abfile.ABFileRestart(psikk_file,"r",idm=idm,jdm=jdm)
         psikk=rfile.read_field("psikk",0)
         thkk=rfile.read_field("thkk",0)
         rfile.close()
         #print(( "psikk.min(),psikk.max()=",psikk.min(),psikk.max())
         #print( "thkk.min(),thkkk.max()=",thkk.min(),thkk.max()
         #
         #montg1pb = modeltools.hycom._montg_tools.montg1_pb(thstar,p)
         #montg1c  = modeltools.hycom._montg_tools.montg1_no_pb(psikk,thkk,thstar,p) 
         montg1pb = montg1_pb(thstar,p)
         montg1c  = montg1_no_pb(psikk,thkk,thstar,p)
         print( "montg1c.min(),montg1c.max()=",montg1c.min(),montg1c.max())
         print( "montg1pb.min(),montg1pb.max()=",montg1pb.min(),montg1pb.max() )
         pbavg  = (srfhgt - montg1c) / (montg1pb+thref)
         print( "pbavg.min(),pbavg.max()=",pbavg.min(),pbavg.max() )
         montg1 = montg1pb*pbavg + montg1c
         logger.info("Estimated montg1 ")
         print( " min max montg1=",montg1.min(),montg1.max() )
         # end estimate montg1---------------------------------------------
         return montg1



def main(archv_files,restart_file,out_path) :
#def main(dtlist,archv_files,restart_file,archv_basedir="") :
#    main([],args.archive_file,archv_basedir=args.archv_basedir)

   # Set paths using REGION.src environment (must be available at python invocation)
   
   # Open blkdat files. Get some properties
   header_line1="HYCOM nested archive files"
   header_line2="Compute montg1 from Nemo ssh"
   header_line3="Expt 01.7  nhybrd=50 nsigma= 0 ds00= 1.00 dp00= 1.00 dp00x= 450.0 dp00f=1.150"
   #
   bp=modeltools.hycom.BlkdatParser("blkdat.input")
   idm    = bp["idm"]
   jdm    = bp["jdm"]
   kdm    = bp["kdm"]
   thflag = bp["thflag"]
   thbase = bp["thbase"]
   kapref = bp["kapref"]
   iversn = bp["iversn"]
   yrflag = bp["yrflag"]
   iexpt  = bp["iexpt"]
   thref=1e-3
   # First read nemo ssh from the meain files
   #archvf_00=archv_files[0][-19:-4]+"00.a"

   for archv_file in archv_files :
      archvf_00=archv_file 
      logger.info("----->Processing archv file :%s"%archvf_00)
      logger.info("----->Restart file _pthkk :%s"%restart_file)
      logger.info("Daily mean nemo file:%s"%archvf_00)
      if os.path.isfile(archvf_00):
         logger.info("file  %s  OK:"%archvf_00)
      else:
         logger.info("Ooobs file does not exist %s:"%archvf_00)
      print( "" )
      ##logger.debug("Reading nemo_srfhgt from %s"%(archvf_00))
      arcfile_in=abfile.ABFileArchv(archvf_00,"r")
      nemo_srfhgt=arcfile_in.read_field("srfhgt",0)
      #nemo_montg1=arcfile_in.read_field("montg1",thflag)
      # ---------------------------------------------------
      logger.info("compute thstar and p for montg1")
      sig  = modeltools.hycom.Sigma(thflag)
      # Read all layers .. (TODO: If there are memory problems, read and estimate sequentially)
      temp  = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      saln  = numpy.ma.zeros((jdm,idm))    # Only needed when calculating density
      th3d  =numpy.ma.zeros((kdm,jdm,idm))
      thstar=numpy.ma.zeros((kdm,jdm,idm))
      dp    =numpy.ma.zeros((jdm,idm))
      p     =numpy.ma.zeros((kdm+1,jdm,idm))
      logger.info("---------------")
      logger.info("-----Reading layers to get thstar and p----------")
      for k in range(kdm) :
         ##logger.debug("Reading layer %d from %s"%(k,archvf_00))
         temp  =arcfile_in.read_field("temp",k+1)
         saln  =arcfile_in.read_field("salin",k+1)
         #dp    [k  ,:,:]=arcfile.read_field("thknss",k+1)
         dp    [:,:]=arcfile_in.read_field("thknss",k+1)
         th3d  [k  ,:,:]=sig.sig(temp,saln)-float(thbase)
         #dens_sig=sig.sig(temp,saln) 
         p     [k+1,:,:]= p[k,:,:] + dp[:,:]
         thstar[k  ,:,:]=numpy.ma.copy(th3d  [k  ,:,:])
         if kapref > 0 :
            thstar[k  ,:,:]=thstar  [k  ,:,:] + kappa.kappaf(
                  temp[:,:], saln[:,:], th3d[k,:,:]+float(thbase), p[k,:,:])
         elif kapref < 0 :
            msg="Only kapref>=0 is implemented for now"
            logger.error(msg)
            raise ValueError(msg)
      # ---------------------------------------------------
      #logger.info("Adding nemo ssh to the tide elevation:%s"%archv_file_in_string)
      logger.info("---------------------------------")
      archv_file_out = out_path+str(archvf_00[-19:])
      logger.info("----->Adding montg1 archv_file_out:%s"%archv_file_out)
      logger.info("Output file %s"%(archv_file_out))
      arcfile_out=abfile.ABFileArchv(archv_file_out,"w",iversn=iversn,yrflag=yrflag,iexpt=iexpt,mask=False,cline1=header_line1,cline2=header_line2,cline3=header_line3)
      #
      ###AA regr = open("montg_regress_tide.pckl",'rb')
      ###AA a,b,nemomeandt = pickle.load(regr)
      ###AA tide_montg1 = (nemo_srfhgt - nemomeandt) * a + b
      estimated_montg1=approx_montg1(thstar,p,nemo_srfhgt,idm,jdm,restart_file)
      #---#
      #read all keys and write updated version
      for keys in sorted(arcfile_in.fields.keys()) :
          fieldname = arcfile_in.fields[keys]["field"]
          time_step = arcfile_in.fields[keys]["step"]
          model_day = arcfile_in.fields[keys]["day"]
          k         = arcfile_in.fields[keys]["k"]
          dens      = arcfile_in.fields[keys]["dens"]
          fld       = arcfile_in.read_field(fieldname,k)
          #logger.info("Current field %10s at level %3d, and dens %s"%(fieldname,k,dens))
          if fieldname == "montg1" :
             logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,archv_file_out))
             arcfile_out.write_field(estimated_montg1,None,fieldname,time_step,model_day,k,dens)
          #elif fieldname == "srfhgt" :
          #   logger.info("Writing field %10s at level %3d to %s (modified)"%(fieldname,k,archv_file_out))
          #   #arcfile_out.write_field(tide_srfhgt,None,fieldname,time_step,model_day,k,dens)
          else :
             arcfile_out.write_field(fld   ,None,fieldname,time_step,model_day,k,dens)
      
      arcfile_in.close()
      arcfile_out.close()
      
      # Remove archive files lying around
   #logger.info("Finished - files placed in location %s"%os.path.abspath("./"))





if __name__ == "__main__" :
   basecmd=os.path.basename(sys.argv[0])

   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)


   parser = argparse.ArgumentParser(description='')
   #parser.add_argument("--archv-basedir",type = str, help = "",default="")
   parser.add_argument('archive_file', type=str, nargs="+")
   parser.add_argument('restart_file', type=str, help="Restart file produced by model run")
   parser.add_argument('out_path',     type=str, help='Output files path')
   #parser_opt2.set_defaults(subparser_name="2")

   args = parser.parse_args()

   main(args.archive_file,args.restart_file,args.out_path)
