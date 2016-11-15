#!/usr/bin/env python
##!/usr/bin/python -E
import modeltools.hycom
import argparse
import datetime
import os


def main(input,output,value) :

   if input == output :
      logger.warning("input format == output format ")
      return value

   # Convert input to dtime
   if input == "dtime" :
      pass
      value=float(value[0])
   elif input == "datetime" :
      value=datetime.datetime.strptime( value[0], "%Y-%m-%dT%H:%M:%S")
      iy,id,ih = modeltools.hycom.datetime_to_ordinal(value,3)
      value=modeltools.hycom.dayfor(iy,id,ih,3)
   elif input == "ordinal" :
      iy = int(value[0])
      id = int(value[1])
      ih = int(value[2])
      value=modeltools.hycom.dayfor(iy,id,ih,3)
   else :
      raise NotImplementedError,input

   if output == "dtime" :
      return value 
   elif output == "datetime" :
      return modeltools.hycom.forday_datetime(value,3).strftime("%Y-%m-%dT%H:%M:%S")

   elif output == "ordinal" :
      return modeltools.hycom.forday(value,3)
   else :
      raise NotImplementedError,output





if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('input',       type=str)
   parser.add_argument('output',     type=str)
   parser.add_argument('value',nargs="+")
   args = parser.parse_args()

   value = main(args.input,args.output,args.value)
   print value



