#!/usr/bin/env python
import matplotlib 
matplotlib.use("Agg")
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


def main_1(s1,s2,s3,Nlin=None,Nexp=None,Ntot=None,dfac=4.0):

   if Ntot is not None :
      tmp = (s2 - s1)*(Ntot) / (dfac * (s3-s2) + s2 - s1)
      print tmp
      Nlin=int(numpy.floor(tmp))
      Nexp=Ntot-Nlin
   else :
      Ntot=Nlin+Nexp

   s=numpy.zeros((Nlin+Nexp))

   # derivative in linear range
   t2=(s2-s1)/(Nlin)

   # Linear derivative in exp range
   t3=(s3-s2)/(Nexp)


   logger.info("Total number of layers            :%d"%Ntot)
   logger.info("Number of layers in linear segment:%d"%Nlin)
   logger.info("Number of layers in exp    segment:%d"%Nexp)


   # 
   print t2,t3
   if t3 >= t2 :
      logger.info("For this routine to work, this conditions must be met:")
      #logger.info("   s3 > s2 > s1")
      logger.info("   (s2-s1)/(Nlin)   >  (s3-s2)/(Nexp)")
      logger.info("Either increase Nexp or reduce Nlin (or modify s1, s2, s3)")
      raise ValueError




   #Linear range 
   s[:Nlin] = numpy.linspace(s1,s1+t2*Nlin,Nlin)

   # The derivative should match that of the exponential range at s2
   #F(N) = A + B*exp(k(N-Nlin))
   #F(0)  = s2
   #F(Ne) = s3
   #df/dN (s2) = t2
   # Yields function for k solved by Newtons method
   def g(k) :
      return (s2-s3) * k/t2 - 1 + numpy.exp(k*Nexp)
   def dgdN(k) :
      return (s2-s3)/t2  + k*numpy.exp(k*Nexp)
   k=-.5 
   for i in range(100):
      tmp0=g(k)
      tmp1=dgdN(k)
      k = k - tmp0/tmp1
      #print k
   B=t2/k
   A=s2-B

   ## derivative of F at Nlin+Nexp :
   #print 
   #print B*k*numpy.exp(k*(Nlin))
   #print B*k*numpy.exp(k*(Nlin+Nexp))

   for i in range(Nlin) :
      print "l",i,s[i]
   for i in range(Nlin,Nexp+Nlin) :
      s[i]=A+B*numpy.exp(k*(i+1-Nlin))
      print "e",i,
      print s[i]


   return s





def main_plot(s) :
   # Plot vertical profile of sig-2 vs depth amd index ?
   f,(ax1,ax2) = plt.subplots(1,2,sharey=False)
   colt="b"
   #ax1.plot(sigma,-intfmid,lw=2,color=colt)
   #ax1.set_ylabel("Min interface depth",color=colt)
   #ax1.set_xlabel("sigma-%d"%thflag,color=colt)
   ax2.plot(s,-numpy.arange(s.size),lw=2,color=colt)
   ax2.set_ylabel("Layer index",color=colt)
   ax2.set_xlabel("sigma",color=colt)
   #ax1.set_ylabel("Temperature[C]",color=colt)
   #ax1.grid(True)
   #for t in ax1.get_xticklabels() :
   #   t.set_color(colt)
   #   t.set_size(8)
   #   t.set_rotation(-45)
   plt.gcf().savefig("tst.png",dpi=180)
#   raise NameError,"test"
     















if __name__ == "__main__" :
   parser = argparse.ArgumentParser(description='')
   subparsers = parser.add_subparsers(help='sub-command help')

   # Option 1 : 
   parser_1 = subparsers.add_parser("lin+exp") 
   parser_1.add_argument('--sigma', default=0)
   parser_1.add_argument('s1', type=float)
   parser_1.add_argument('s2', type=float)
   parser_1.add_argument('s3', type=float)
   parser_1.add_argument('Nlin', type=int,  help="Number of points in linear range[s1,s2]")
   parser_1.add_argument('Nexp', type=int,  help="Number of points in exponential range(s2,s3]")
   parser_1.set_defaults(subparser_name="lin+exp")

   parser_2 = subparsers.add_parser("lin+expv2",description=
         "Linear + exp profile. Assumption: (s2-s1)/Nlin = dfac*(s3-s2)/Nexp. Nlin and Nexp are estimated from input" ) 
   parser_2.add_argument('--dfac', type=float, help="Ratio of (s2-s1)*Nlin  : (s3-s1)*Nexp",default=4.0)
   parser_2.add_argument('--sigma', default=0)
   parser_2.add_argument('s1', type=float, help="Surface density")
   parser_2.add_argument('s2', type=float, help="density (end of linear range)")
   parser_2.add_argument('s3', type=float, help="Bottom density")
   parser_2.add_argument('Ntot', type=int,  help="Number of points in total range")
   parser_2.set_defaults(subparser_name="lin+expv2")

   args = parser.parse_args()
   if args.subparser_name == "lin+exp" :
      s=main_1(args.s1,args.s2,args.s3,Nlin=args.Nlin,Nexp=args.Nexp)
   elif args.subparser_name == "lin+expv2" :
      s=main_1(args.s1,args.s2,args.s3,Ntot=args.Ntot,dfac=args.dfac)
   else :
      raise ValueError,"Unknown parser "
   main_plot(s)



