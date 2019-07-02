#!/usr/bin/env python
##!/usr/bin/python -E
import modeltools.hycom
import argparse
import datetime
import os


def main(start_time,end_time,restart) :
   bp=modeltools.hycom.BlkdatParser("blkdat.input")
   yrflag=bp["yrflag"]

   if restart : 
      rfac=-1
   else :
      rfac=1

   fday,fhour,fsecond = modeltools.hycom.datetime_to_ordinal(start_time,bp["yrflag"])
   lday,lhour,lsecond = modeltools.hycom.datetime_to_ordinal(end_time,  bp["yrflag"])
   #print start_time,fday,fhour,fsecond,start_time.year

   # NB Add 1 to fday, ldat - they start from 0 but hycom
   # starts from 1
   fdtime=modeltools.hycom.dayfor(start_time.year,fday,fhour,yrflag)+fsecond/86400.
   ldtime=modeltools.hycom.dayfor(end_time.year  ,lday,lhour,yrflag)+lsecond/86400.
   #print fdtime, start_time,fday,fhour,fsecond,start_time.year

   f = open("limits", "w")
   f.write( "%14.5f %14.5f" % ( (rfac*fdtime), ldtime ))
   f.close

if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--init',       action="store_true", default=False)
   parser.add_argument('start_time', action=DateTimeParseAction, help='Start time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser.add_argument('end_time',   action=DateTimeParseAction, help='Stop  time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   args = parser.parse_args()

   main(args.start_time,args.end_time,args.init)



