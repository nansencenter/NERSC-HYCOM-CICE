#!/usr/bin/env python
##!/usr/bin/python -E
import matplotlib.pyplot as plt
import modeltools.hycom
import logging
import argparse
import datetime
import numpy
import os

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


def main(s1,s2,s3,Nlin,Nexp):





if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   subparsers = parser.add_subparsers(help='sub-command help')

   # Option 1 : 
   parser_1 = subparser.add_parser("lin+exp") 
   parser_1.add_argument('s1')
   parser_1.add_argument('s2')
   parser_1.add_argument('s3')
   parser_1.add_argument('Nlin', help="Number of points in linear range[s1,s2)")
   parser_1.add_argument('Nexp', help="Number of points in exponential range[s2,s3]")

   parser.add_argument('file' , help="blkdat file")
   args = parser.parse_args()
   main(args.file)



