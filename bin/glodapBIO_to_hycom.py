import numpy
import os
import math
import argparse

'''
CAGLAR AUG 2020 ----
HOW TO USE:
locate the river flux file and assign to "volumefile"
locate the GLODAPv2 nutrient load file and assign to "nutfile"
save the output file in git folder as " NERSC-HYCOM-CICE/input/biorivers.dat"
river scripts look for that exact output location and biorivers.dat
works in any directory

script execute:
python glodapBIO_to_hycom.py volumefile nutfile hycompath radius

example:
python glodapBIO_to_hycom.py /cluster/home/cagyum/DEVEL/NERSC-HYCOM-CICE/input/rivers_ahype-ehype_clim_rev2.dat /cluster/home/cagyum/DEVEL/NERSC-HYCOM-CICE/input/globalNEWS_2012.dat /cluster/home/cagyum/DEVEL/NERSC-HYCOM-CICE/ 20
'''

def distance_on_unit_sphere(lat1, long1, lat2, long2):

# Convert latitude and longitude to
# spherical coordinates in radians. Assumes the Earth is a perfect sphere.
  degrees_to_radians = math.pi/180.0
# phi = 90 - latitude
  phi1 = (90.0 - lat1)*degrees_to_radians
  phi2 = (90.0 - lat2)*degrees_to_radians
# theta = longitude
  theta1 = long1*degrees_to_radians
  theta2 = long2*degrees_to_radians

  cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
  math.cos(phi1)*math.cos(phi2))
  arc = math.acos( cos )
  arc = arc * 6378.137 # convert to km

  return arc

def main(volumefile,nutfile,hycompath,radius):
   searchradius = radius # kilometers
   river_volume_flux_file   = open(volumefile) 
   river_nutrient_flux_file = open(nutfile) 

   fflux      = river_volume_flux_file.readlines()
   nfluxlines = len(fflux)/3 # calculates number of rivers in flux file
   river_volume_flux_file   = open(volumefile)  # open again for colocation

   class DataFlux:
      name    = ["" for x in range(nfluxlines)]
      onoff   = ["" for x in range(nfluxlines)]
      lat     = numpy.zeros(nfluxlines)
      lon     = numpy.zeros(nfluxlines)
      annual  = numpy.zeros(nfluxlines)
      monthly = numpy.zeros((nfluxlines,12))


   k=0 # read in flux information
   for line in river_volume_flux_file: 
    if   ( line.split()[0] == 'T' ) or ( line.split()[0] == 'F' ) :
              DataFlux.onoff[k] = line.split()[0]
              DataFlux.name[k] = line.split()[1]
    elif ( len(line.split()) > 2 ) :
              DataFlux.annual[k] = float( line.split()[0] )
              for x in range(1,len(line.split())) :
                       DataFlux.monthly[k,x-1] = float( line.split()[x] )
    else :
              DataFlux.lat[k] = float( line.split()[0] )
              DataFlux.lon[k] = float( line.split()[1] )
              k=k+1
   river_volume_flux_file.close()


   fload      = river_nutrient_flux_file.readlines() # open again for colocation
   nloadlines = len(fload)-1

   class DataLoad:
      name = ["" for x in range(nloadlines)]
      lat  = numpy.zeros(nloadlines)
      lon  = numpy.zeros(nloadlines)
      nit  = numpy.zeros(nloadlines)
      pho  = numpy.zeros(nloadlines)
      sil  = numpy.zeros(nloadlines)
      don  = numpy.zeros(nloadlines)

   k=0 # read in nutrient information
   with open(nutfile) as river_nutrient_flux_file:
     first_line = river_nutrient_flux_file.readline()
     for line in river_nutrient_flux_file:
         DataLoad.name[k] = line.split()[1]
         DataLoad.lon[k]  = float( line.split()[4] )
         DataLoad.lat[k]  = float( line.split()[5] )
         DataLoad.nit[k]  = float( line.split()[17] )
         DataLoad.pho[k]  = float( line.split()[18] )
         DataLoad.don[k]  = float( line.split()[19] )
         DataLoad.sil[k]  = float( line.split()[30] ) / 1000.
         k=k+1

   class Colocated:
      fluxname = ["" for x in range(nfluxlines)]
      loadname = ["" for x in range(nfluxlines)]
      onoff    = ["" for x in range(nfluxlines)]
      fluxlon  = numpy.zeros(nfluxlines)
      fluxlat  = numpy.zeros(nfluxlines)
      fluxann  = numpy.zeros(nfluxlines)
      loadlon  = numpy.zeros(nfluxlines)
      loadlat  = numpy.zeros(nfluxlines)
      loadnit  = numpy.zeros(nfluxlines)
      loadpho  = numpy.zeros(nfluxlines)
      loadsil  = numpy.zeros(nfluxlines)
      loaddon  = numpy.zeros(nfluxlines)

   outfile = open(hycompath + "input/biorivers.dat","w")
   for x in range(nfluxlines):
      distance = numpy.zeros(nloadlines)
      for y in range(nloadlines):
            distance[y] = distance_on_unit_sphere( DataFlux.lat[x], DataFlux.lon[x], DataLoad.lat[y], DataLoad.lon[y] ) 
      found = numpy.where(distance<=searchradius)
      if numpy.size(found) != 0:
         message = 'For: (Lat:%s,Lon:%s) found %s points: Lat%s, Lon%s within %skm radius'
         print( message%(DataFlux.lat[x],DataFlux.lon[x],numpy.size(found),DataLoad.lat[found],DataLoad.lon[found],searchradius) )
         Colocated.fluxlat[x]  = DataFlux.lat[x]
         Colocated.fluxlon[x]  = DataFlux.lon[x]
         Colocated.fluxname[x] = DataFlux.name[x]
         Colocated.onoff[x]    = DataFlux.onoff[x]
         Colocated.fluxann[x]    = DataFlux.annual[x]
         Colocated.loadlon[x]  = numpy.mean( DataLoad.lon[found] )
         Colocated.loadlat[x]  = numpy.mean( DataLoad.lat[found] )
         Colocated.loadnit[x]  = numpy.max( DataLoad.nit[found] ) * 1E9 * (6.625*12.01/14.01) / (365.25*24.0*60.0*60.0) # MgN/yr  --> mgC/s 
         Colocated.loadpho[x]  = numpy.max( DataLoad.pho[found] ) * 1E9 * (106.0*12.01/30.97) / (365.25*24.0*60.0*60.0) # MgP/yr  --> mgC/s
         Colocated.loadsil[x]  = numpy.max( DataLoad.sil[found] ) * 1E9 * (6.625*12.01/28.09) / (365.25*24.0*60.0*60.0) # MgSi/yr --> mgC/s
         Colocated.loaddon[x]  = numpy.max( DataLoad.don[found] ) * 1E9 * (6.625*12.01/14.01) / (365.25*24.0*60.0*60.0) # MgN/yr  --> mgC/s
      else :
         Colocated.fluxlat[x]  = DataFlux.lat[x]
         Colocated.fluxlon[x]  = DataFlux.lon[x]
         Colocated.fluxname[x] = DataFlux.name[x]
         Colocated.onoff[x]    = DataFlux.onoff[x]
         Colocated.fluxann[x]   = DataFlux.annual[x]
         Colocated.loadlon[x]  = 0.#-999. 
         Colocated.loadlat[x]  = 0.#-999.
         Colocated.loadnit[x]  = 0.#-999.
         Colocated.loadpho[x]  = 0.#-999.
         Colocated.loadsil[x]  = 0.#-999.
         Colocated.loaddon[x]  = 0.#-999.

   for x in range(nfluxlines):
      found2 = numpy.where((Colocated.loadlon == Colocated.loadlon[x]) & (Colocated.loadlat == Colocated.loadlat[x]))
      if (Colocated.loadlon[x] != 0.) & (Colocated.loadlat[x] != 0.) :
        if numpy.size(found2)>1 : 

           print('found ' + numpy.str(numpy.size(found2)) + ' locations for the same bioriver point @indexes ' + numpy.str(found2))

           nloc = numpy.size(found2)
           distloc = numpy.zeros((nloc))
           for dd in range(nloc):
               distloc[dd] = 1./distance_on_unit_sphere( Colocated.fluxlat[found2[0][dd]], Colocated.fluxlon[found2[0][dd]], Colocated.loadlat[found2[0][dd]], Colocated.loadlon[found2[0][dd]] )
           Colocated.loadnit[found2[0][found2[0]==x]] = distloc[found2[0]==x]/distloc.min() / ( numpy.sum( distloc/distloc.min() ) ) * Colocated.loadnit[found2[0][found2[0]==x]]
           Colocated.loadpho[found2[0][found2[0]==x]] = distloc[found2[0]==x]/distloc.min() / ( numpy.sum( distloc/distloc.min() ) ) * Colocated.loadpho[found2[0][found2[0]==x]]
           Colocated.loadsil[found2[0][found2[0]==x]] = distloc[found2[0]==x]/distloc.min() / ( numpy.sum( distloc/distloc.min() ) ) * Colocated.loadsil[found2[0][found2[0]==x]]

      outfile.write('{:3}'.format(Colocated.onoff[x]) + Colocated.fluxname[x]  + os.linesep)
      outfile.write(str( '{:10.1f}'.format(Colocated.fluxann[x])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][0])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][1])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][2])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][3])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][4])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][5])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][6])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][7])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][8])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][9])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][10])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][11])  ) \
                                                                     + os.linesep)
      outfile.write(str( '{:11.2f}'.format(Colocated.loadnit[x])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][0])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][1])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][2])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][3])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][4])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][5])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][6])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][7])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][8])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][9])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][10])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][11])  ) \
                                                                     + os.linesep)
      outfile.write(str( '{:11.2f}'.format(Colocated.loadpho[x])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][0])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][1])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][2])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][3])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][4])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][5])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][6])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][7])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][8])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][9])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][10])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][11])  ) \
                                                                     + os.linesep)
      outfile.write(str( '{:11.2f}'.format(Colocated.loadsil[x])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][0])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][1])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][2])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][3])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][4])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][5])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][6])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][7])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][8])  ) \
                                                                     + str( '{:8.5f}'.format(DataFlux.monthly[x][9])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][10])  ) + str( '{:8.5f}'.format(DataFlux.monthly[x][11])  ) \
                                                                     + os.linesep)
      outfile.write( str( '{:9.2f}'.format(Colocated.fluxlat[x])  ) + str( '{:9.2f}'.format(Colocated.fluxlon[x]) ) + os.linesep)
   outfile.close()

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('volumefile',  type=str)
    parser.add_argument('nutfile',  type=str)
    parser.add_argument('hycompath',  type=str)
    parser.add_argument('radius',  type=float)
    args = parser.parse_args()

    main(args.volumefile,args.nutfile,args.hycompath,args.radius)
