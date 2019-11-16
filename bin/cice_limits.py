#!/usr/bin/env python
##!/usr/bin/python -E
import modeltools.hycom
import f90nml
import argparse
import datetime
import numpy
import os


def main(start_time,end_time,init,nmpi,fnml) :
   #fnml = "ice_in"
   nml  = f90nml.read(fnml)
   dt   = nml["setup_nml"]["dt"]
   year_init   = nml["setup_nml"]["year_init"]

   # Integration time, initial year
   tint      = end_time - start_time
   tint_secs = tint.days*86400. + tint.seconds
   #year_init = start_time.year
   #year_init = 1958
   npt=int(numpy.floor(tint_secs/dt))+1
   
   # Seconds lapsed into this year
   sec0=start_time - datetime.datetime(year_init,1,1,0,0,0)
   sec0=sec0.days*86400. + sec0.seconds
   istep0=int(numpy.floor(sec0/dt))-1

   if istep0 < 0 :
      #year_init=year_init-1
      #sec0=start_time - datetime.datetime(year_init,1,1,0,0,0)
      #sec0=sec0.days*86400. + sec0.seconds
      #istep0=int(numpy.floor(sec0/dt))-1
      raise NameError,"negative istep0, adjust year_init"


   nml["setup_nml"]["year_init"] = year_init
   nml["setup_nml"]["istep0"] = istep0
   nml["setup_nml"]["npt"]=npt
   print "year_init = ",nml["setup_nml"]["year_init"]
   print "istep0    = ",nml["setup_nml"]["istep0"]
   print "npt       = ",nml["setup_nml"]["npt"]

   nml["domain_nml"]["nprocs"]=nmpi

   if init :
      nml["setup_nml"]["runtype"]="initial"
      nml["setup_nml"]["restart"]=False
      nml["setup_nml"]["use_restart_time"]=False
      nml["setup_nml"]["ice_ic"]="default"
   else :
      nml["setup_nml"]["runtype"]="continue"
      nml["setup_nml"]["restart"]=True
      nml["setup_nml"]["use_restart_time"]=True

   nml.write("ice_in",force=True)

   path_restart= os.path.join("./", nml["setup_nml"]["restart_dir"])
   if not os.path.isdir(path_restart) :
      os.mkdir(path_restart)

   path_pointer= os.path.join("./", os.path.dirname(nml["setup_nml"]["pointer_file"]))
   if not os.path.isdir(path_pointer) :
      os.mkdir(path_pointer)

if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)
   parser = argparse.ArgumentParser(description='Sets up a new fortran namelist from CICE infile based on provided arguments')
   parser.add_argument('--init',       action="store_true", default=False)
   parser.add_argument('start_time', action=DateTimeParseAction, help='Start time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser.add_argument('end_time',   action=DateTimeParseAction, help='Stop  time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser.add_argument('nmpi',       type=int,help="Number of mpi tasks to use")
   parser.add_argument('infile',     help='CICE namelist')
   args = parser.parse_args()
   main(args.start_time,args.end_time,args.init,args.nmpi,args.infile)



