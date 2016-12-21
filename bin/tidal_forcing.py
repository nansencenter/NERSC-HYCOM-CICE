#!/usr/bin/env python
import hycom_add_tides
import modeltools.hycom
import subprocess
import argparse
import datetime
import os.path
import sys
import re
import logging

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


def main(dtlist,archv_files,include_uv=False,tidal_database="FES2014",archv_basedir="") :

   # Set paths using REGION.src environment (must be available at python invocation)
   if tidal_database == "FES2014" :
      tide_cmd = os.path.join(os.environ["MSCPROGS"],"bin_setup","fes2014hycom")
      tide_file = "fes2014_h.nc"
      tide_file_u = "fes2014_u.nc"
      tide_file_v = "fes2014_v.nc"
   elif tidal_database == "FES2004" :
      tide_cmd = os.path.join(os.environ["MSCPROGS"],"bin_setup","fes2004hycom")
      tide_file = "fes2004_h.nc"
   else :
      msg = "Unknown tidal database %s"%tidal_database
      logger.error(msg)
      raise ValueError,msg

   # currents only for FES2104
   if tidal_database == "FES2004" and include_uv :
      msg = "FES 2004 does not have currents. Ignoring option include_uv"
      include_uv=False
      logger.warning(msg)

   # If archv_basedir is provided, any paths are relative to archv_basedir (unless full path is provided)
   for i in range(len(archv_files)) :

      if archv_files[i][0] == "/" : 
         pass
      elif archv_basedir :
         archv_files[i] = os.path.join(archv_basedir,archv_files[i])


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

   # Extract time from archive file name if provided
   a_datetimes=[]
   for a in archv_files :
      m=re.match("archv.([0-9]{4})_([0-9]{3})_([0-9]{2}).*",os.path.basename(a))
      if m :
         year=int(m.group(1))
         oday=int(m.group(2))
         hour=int(m.group(3))
         dtime = modeltools.hycom.dayfor(year,oday,hour,yrflag)
         mydatetime=modeltools.hycom.forday_datetime(dtime,yrflag)
         a_datetimes.append(mydatetime)
      else  :
         msg="Unable to extract time from archive file name %s"%archv_file
         logger.error(msg)
         raise ValueError,msg


   #1) A datetime list is provided. Use that for calculating tides
   if dtlist :
      tide_forcing_dt=dtlist
   #1) A archive file list is provided. Use that for calculating tides
   elif archv_files :
      tide_forcing_dt=a_datetimes

   # Transform times to be relative to 1950-01-01 00:00:00 UTC. This is used by the FES C API
   reftime=datetime.datetime(1950,1,1,0,0,0)
   jdays=[ (elem-reftime).total_seconds()/86400. for elem in tide_forcing_dt ]
   logger.info("Calculating %s tides for julian days relative to %s : %s"%(tidal_database,str(reftime),str(jdays)))

   # Remove old tide files lying around
   if os.path.exists(tide_file) : os.unlink(tide_file)

   # Now call tidal generation routine
   #print ["%16.8f"%elem for elem in jdays]
   rc = subprocess.call([tide_cmd,"h"] + ["%16.8f"%elem for elem in jdays])
   if rc <> 0 :
      msg = "An error occured when running %s for h. Return code=%d"%(tide_cmd,rc)
      logger.error(msg)
      raise ValueError,msg


   if include_uv :
      if os.path.exists(tide_file_u) : os.unlink(tide_file_u)
      if os.path.exists(tide_file_v) : os.unlink(tide_file_v)
      rc=subprocess.call([tide_cmd,"u"] + ["%16.8f"%elem for elem in jdays])
      if rc <> 0 :
         msg = "An error occured when running %s for u. Return code=%d"%(tide_cmd,rc)
         logger.error(msg)
         raise ValueError,msg
      rc=subprocess.call([tide_cmd,"v"] + ["%16.8f"%elem for elem in jdays])
      if rc <> 0 :
         msg = "An error occured when running %s for v. Return code=%d"%(tide_cmd,rc)
         logger.error(msg)
         raise ValueError,msg

   # Remove archive files lying around
   path0="./archv_with_tide/"
   if os.path.exists(path0) :
      for file in os.listdir(path0) :
         fullpath= os.path.join(path0,file)
         if os.path.isfile(fullpath) : os.unlink(fullpath)


   # Tidal forcing generated - Now call hycom_add_tides, which will create archive files, or add to existing 
   # archive file
   hycom_add_tides.main(tide_file,archv_files,include_uv=include_uv)

   # Move files into parent dir
   for file in os.listdir(path0) :
      fullpath= os.path.join(path0,file)
      if os.path.isfile(fullpath) : 
         newpath= os.path.join("..",file)
         os.rename(fullpath,newpath)









if __name__ == "__main__" :
   basecmd=os.path.basename(sys.argv[0])
   usage="""
   This routine will use FES to calculate the didal forcing on the boundary of your
   model domain. The tidal forcing will then be used to calculate hycom archive
   files suitable for use as nesting files.

   There are two ways to call this routine:

   1)Input to this routine is a bunch of archive files. The script will extract the
   times of the archive files, provide tidal forcing for those, and add to the
   archive files.

   2) Input to this routine is a start time, an end time, and a time step in days.
   Tidal forcing will be generated from start time to end time with the provided
   time step. Any archive files given on the command line will be modified to take
   the tidal forcing into account IF the times match those of the tidal data



   Usage 1:
      {0} [--include-uv] 1 archive_file1 [archive_file2 ...]

   Usage 2:
      {0} [--include-uv] 2 start_time end_time delta_t [archive_file1 archive_file2 ..]


   Arguments:
      archive_file : a HYCOM archv file
      tidal_dataset: FES2004 or FES2014
      start_time   : (Only usage 2) Start of tide generation
      end_time     : (Only usage 2) End of tide generation
      delta_t      : (Only usage 2) Time step in hours

   Optional arguments
      --tidal-database dbase  : Tidal database to use; FES2004 or FES2014(default)
      --include-uv            : Include tidal currents (only FES2014)


   Examples:
      {0} 1 2013-01-01T00:00:00 2013-01-10T00:00:00 1.0 
      {0} 1 2013-01-01T00:00:00 2013-01-10T00:00:00 1.0  archv.2013_003_12.a archv.2013_004.a 
      {0} 2 archv.2013_003_12.a archv.2013_004.a 
   """.format(basecmd)

   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)


   parser = argparse.ArgumentParser(description='')
   parser.add_argument("--include_uv",action="store_true", default=False, help="Include tidal currents")
   parser.add_argument("--archv-basedir",type = str, help = "",default="")
   parser.add_argument("--tidal-database",type = str, default="FES2014", help = "tidal dataset to use  FES2004 or FES2014(default)")
   subparsers = parser.add_subparsers(help='sub-command help')

   parser_opt1 = subparsers.add_parser('1', help='')
   parser_opt1.set_defaults(subparser_name="1")
   parser_opt1.add_argument('start_time', action=DateTimeParseAction, help='Start time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser_opt1.add_argument('end_time', action=DateTimeParseAction, help='End time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser_opt1.add_argument('deltah', type=float, help='Time step in hours')
   parser_opt1.add_argument('archive_file', type=str, nargs="*")

   parser_opt2 = subparsers.add_parser('2', help='')
   parser_opt2.set_defaults(subparser_name="2")
   parser_opt2.add_argument('archive_file', type=str, nargs="+")
   args = parser.parse_args()


   if args.subparser_name == "1" :

      dtlist=[]
      mytime=args.start_time
      while mytime <= args.end_time :
         dtlist.append(mytime)
         mytime=mytime+ datetime.timedelta(hours=args.deltah)
      main(dtlist,args.archive_file,
            include_uv=args.include_uv,
            tidal_database=args.tidal_database,
            archv_basedir=args.archv_basedir)
   elif args.subparser_name == "2" :
      main([],args.archive_file,
            include_uv=args.include_uv,
            tidal_database=args.tidal_database,
            archv_basedir=args.archv_basedir)
