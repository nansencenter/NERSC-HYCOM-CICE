#!/usr/bin/env python
##!/usr/bin/python -E
import modeltools.hycom
import argparse
import datetime
import os
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
"""
Ex:
Converting ['2020-03-01T00:00:00'] from datetime to ordinal
python ../bin/hycom_date.py datetime ordinal 2020-03-01T00:00:00
(2020, 61, 0)

Ex:
Converting ['2020', '299', '00'] from ordinal to datetime
python ../bin/hycom_date.py ordinal datetime 2020 299 00
43763.0
2020 299 0
2020-10-25T00:00:00

Ex:
Converting ['43726'] from dtime to datetime   
python ../bin/hycom_date.py dtime datetime  43726
43726.0
2020 262 0
2020-09-18T00:00:00
"""

def main(input,output,value) :

   logger.info("Converting %s from %s to %s"%(str(value),input,output))

   if input == output :
      logger.warning("input format == output format ")
      return value

   # Convert input to dtime
   if input == "dtime" :
      pass
      value=float(value[0])
   elif input == "datetime" :
      value=datetime.datetime.strptime( value[0], "%Y-%m-%dT%H:%M:%S")
      iy=value.year
      id,ih,isec = modeltools.hycom.datetime_to_ordinal(value,3)
      value=modeltools.hycom.dayfor(iy,id,ih,3)
   elif input == "ordinal" :
      iy = int(value[0])
      id = int(value[1])
      ih = int(value[2])
      value=modeltools.hycom.dayfor(iy,id,ih,3)
   else :
      logger.error("input must be either dtime, datetime or ordinal")
      raise NotImplementedError(input)

   if output == "dtime" :
      return value 
   elif output == "datetime" :
      return modeltools.hycom.forday_datetime(value,3).strftime("%Y-%m-%dT%H:%M:%S")

   elif output == "ordinal" :
      return modeltools.hycom.forday(value,3)
   else :
      logger.error("output must be either dtime, datetime or ordinal")
      raise NotImplementedError(input)





if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='')
   parser.add_argument('input',       type=str, 
      help="Input format; dtime (hycom), datetime (actual date), ordinal (year ddd hh)")
   parser.add_argument('output',     type=str,
      help="Output format; dtime (hycom), datetime (actual date), ordinal (year ddd hh)")
   parser.add_argument('value',nargs="+")
   args = parser.parse_args()

   value = main(args.input,args.output,args.value)
   print(value)
