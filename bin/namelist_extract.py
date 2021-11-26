#!/usr/bin/env python
##!/usr/bin/python -E
import modeltools.hycom
import f90nml
import argparse
import datetime
import numpy
import os


def main(fnml,block,param) :
   nml  = f90nml.read(fnml)

   print(nml[block][param])

if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   parser.add_argument('fnml',        type=str)
   parser.add_argument('block',       type=str)
   parser.add_argument('param',       type=str)
   args = parser.parse_args()
   main(args.fnml,args.block,args.param)



