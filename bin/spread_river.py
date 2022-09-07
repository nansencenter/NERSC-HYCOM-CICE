import abfile
import numpy as np
import math
import argparse
import os

'''
Author: Veli Çağlar Yumruktepe
Date  : June 21st, 2022

This script spreads the Ob river nutrient fluxes towards the outer bay.
It spreads with an inverse of distance to source point, as such,
   there is a gradient with a decrease towards the outer bay.
The mass is conserved. Follow the print screen to see the old and new mass.


Usage:
python spread_river.py /cluster/work/users/cagyum/TP5a0.06/force/rivers/010/ 

Output:
Filenames will not be changed but a copy with a suffix: _orig will be preserved.

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

def main(path):
    print(path)


    # coordinates for the river mouth based on the grid location on TP5 domain
    river_lat = 66.37281
    river_lon = 71.15888 
    # coordinates for the original river discharge effective area
    #(lat,lon)_1 = 71.496376,72.00326
    #(lat,lon)_2 = 70.26682,82.36286 
    #(lat,lon)_3 = 62.861946,67.83625
    #(lat,lon)_4 = 62.23391,75.26198

    abgrid = abfile.ABFileGrid(path + "../../../topo/regional.grid","r")
    plon=abgrid.read_field("plon")
    plat=abgrid.read_field("plat")
    jdm,idm=plon.shape

    abdepth = abfile.ABFileBathy(path+"../../../topo/regional.depth.b","r",idm=idm,jdm=jdm)
    depthm=abdepth.read_field("depth")

    #
    cooINDEX = abs( plat-river_lat ) + abs( plon-river_lon )
    la,lo = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
    #
    cooINDEX = abs( plat-71.496376 ) + abs( plon-72.00326 )
    la1,lo1 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
    cooINDEX = abs( plat-70.26682 ) + abs( plon-82.36286 )
    la2,lo2 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
    cooINDEX = abs( plat-62.861946 ) + abs( plon-67.83625 )
    la3,lo3 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)
    cooINDEX = abs( plat-62.23391 ) + abs( plon-75.26198 )
    la4,lo4 = np.unravel_index(cooINDEX.argmin(), cooINDEX.shape)


    dist = np.zeros((plon.shape))
    for j in range(plon.shape[0]-1):
       for i in range(plon.shape[1]-1):
           dist[j,i] = distance_on_unit_sphere(plat[j,i],plon[j,i],plat[la,lo],plon[la,lo])
    dist = np.ma.masked_where(depthm<0.,dist)


    dist = np.ma.masked_where(plat>=74.,dist)
    dist = np.ma.masked_where(plon<=70.,dist)
    dist = np.ma.masked_where(plon>=90.,dist)
    dist = dist/dist.max()
    dist = np.abs(dist-1.)
    totaldist = np.sum(dist)
    # now modify rivers for 3 nutrients
    for varib in ['no3','pho','sil']:
       if varib == 'no3': longname = 'nitrate'
       if varib == 'pho': longname = 'phosphate'
       if varib == 'no3': longname = 'silicate'
       outfile=abfile.ABFileRiver(path + "ECO_"+varib+"_new.a","w",idm=idm,jdm=jdm,\
                   cline1='River '+longname+' fluxes (Ob River dispersed out the bay)',\
                   cline2='mgC m-2 s-1')
       outfile.write_header()
       abriver = abfile.AFile(idm,jdm,path + "ECO_"+varib+".a","r")
       for month in range(12):
           river  = abriver.read_record(month)
           river_modified  = abriver.read_record(month)
           river = np.ma.masked_where(depthm<0.,river)
           river_modified = np.ma.masked_where(depthm<0.,river_modified)

           unit = np.sum(river[la1:la2,lo2:lo4])/totaldist
           final = unit*dist

           river_modified[la1:la2,lo2:lo4][ ~river_modified[la1:la2,lo2:lo4].mask ] = 0.
           river_modified[~final.mask] = river_modified[~final.mask] + final[~final.mask]
           print("original total mass: "+str(np.sum(river))+" modified total mass: "+str(np.sum(river_modified)))
           outfile.write_field(river_modified,None,"river "+longname,month+1)

       abriver.close()
       outfile.close()

       origAfile = path + "ECO_"+varib+".a"
       origBfile = path + "ECO_"+varib+".b"
       oldAfile  = path + "ECO_"+varib+"_orig.a"
       oldBfile  = path + "ECO_"+varib+"_orig.b"
       newAfile  = path + "ECO_"+varib+"_new.a"
       newBfile  = path + "ECO_"+varib+"_new.b"

       os.rename(origAfile,oldAfile)
       os.rename(origBfile,oldBfile)
       os.rename(newAfile,origAfile)
       os.rename(newBfile,origBfile)

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('path',  type=str, help="path to river forcing files")
    args = parser.parse_args()

    main(args.path)

