#!/usr/bin/python

def main():
   import sys
   if len(sys.argv)!= 6 :
      print "Routine interpolates between two longitude, latitude points to get "
      print "station positions between them. Interpolation is along rhumb lines  "
      print ""
      print "Script requires 5 arguments - first_lon  first_lat last_lon last_lat number_of_stations"
      sys.exit(-1)
   tmp,flon, flat, llon, llat, numstat =sys.argv
   numstat=int(numstat)
   flon=float(flon); llon=float(llon); 
   flat=float(flat); llat=float(llat); 
   for i in range(int(numstat)) :
      print  "%10.3f %10.3f" % (
      flon*(numstat-i-1)/(numstat-1) + llon*i/(numstat-1) ,
      flat*(numstat-i-1)/(numstat-1) + llat*i/(numstat-1) 
      )

main()
