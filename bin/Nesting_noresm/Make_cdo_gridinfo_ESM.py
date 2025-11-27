#!/usr/bin/env python
# This python routine generates a grid information for the ESM so that NetCDF-files can be used with CDO.
# For example to make spatial aveages from netcdf-files or remapping:
#
#cdo remapbil,cdogridinfofile satelite_fld.nc sat_fld_mapped_to_esm.nc 
#
# usage: python Make_cdo_gridinfo_ESM.py ESM_grid_file(full path)
# eg.g.
#./Make_cdo_gridinfo_ESM.py areacello_Ofx_CNRM-ESM2-1_historical_r1i1p1f2_gn.nc --esm_id=CNRM-ESM2-1
# output: cdogrid_CNRM-ESM2-1


import numpy as np
from netCDF4 import Dataset
import argparse


def main(meshfile,esm_id):
    # read in the necessary grid variables
    
    fh = Dataset(meshfile, mode='r')
    print(fh.variables)
    print(esm_id)
    if esm_id=="CNRM-ESM2-1":
        plon=fh.variables["lon"][:]
        plat=fh.variables["lat"][:]
        scpx=fh.variables["bounds_lon"][:]
        scpy=fh.variables["bounds_lat"][:]
    else:
        plon=fh.variables["longitude"][:] # variable names not always consisten between ESMS
        plat=fh.variables["latitude"][:]
        scpx=fh.variables["vertices_longitude"][:]
        scpy=fh.variables["vertices_latitude"][:]
    
    fh.close()

    print(np.shape(scpy))

    if esm_id=="":
        outfile="cdogrid_ESM"
    else:
        outfile="cdogrid_"+esm_id
    

    #write the grid information to the gridinfofile
    fileout = open(outfile,"w")

    idm,jdm=np.shape(np.transpose(plon))
    print(idm,jdm,idm*jdm,(idm-1)*(jdm-1))
    fileout.write('gridtype = curvilinear\n');
    fileout.write('gridsize = ' + str((idm)*(jdm)) + '\n');
    fileout.write('xsize    = ' + str(idm) + '\n');
    fileout.write('ysize    = ' + str(jdm) + '\n');

    fileout.write('xvals    = ');
    for i in range(np.size(plon,0)):
        for j in range(np.size(plon,1)):
            fval = "{:.3f}".format(plon[i,j])
            fileout.write(fval + " ") 
        fileout.write('\n');   
    fileout.write('\n'); 

    fileout.write('yvals    = ');
    for i in range(np.size(plon,0)):
        for j in range(np.size(plon,1)):
            fval = "{:.3f}".format(plat[i,j])
            fileout.write(fval + " ") 
        fileout.write('\n'); 
    fileout.write('\n');

    fileout.write('nvertex  = 4\n');

    fileout.write('xbounds    = ');
    for i in range(np.size(plon,0)):
        for j in range(np.size(plon,1)):

            fval1 = "{:.3f}".format(scpx[i,j,0])
            fval2 = "{:.3f}".format(scpx[i,j,1]) 
            fval3 = "{:.3f}".format(scpx[i,j,2]) 
            fval4 = "{:.3f}".format(scpx[i,j,3]) 
            fileout.write(fval1 + " " + fval2 + " " + fval3 + " " + fval4 + "\n")

    fileout.write('\n');        

    fileout.write('ybounds    = ');
    for i in range(np.size(plon,0)):
        for j in range(np.size(plon,1)):

            fval1 = "{:.3f}".format(scpy[i,j,0])
            fval2 = "{:.3f}".format(scpy[i,j,1]) 
            fval3 = "{:.3f}".format(scpy[i,j,2]) 
            fval4 = "{:.3f}".format(scpy[i,j,3]) 
            fileout.write(fval1 + " " + fval2 + " " + fval3 + " " + fval4 + "\n")  
    fileout.write('\n');

    fileout.close()


if __name__ == "__main__" :

   parser = argparse.ArgumentParser(
         description='This python routine generates a grid information for the ESM so that NetCDF-files can be used with CDO'
         )
   parser.add_argument('meshfile', type=str,help="ESM mesh file in netcdf format - full path")
   parser.add_argument('--esm_id', type=str,default="",help="ESM model identifier, eg: CNRM-ESM2-1")
   args = parser.parse_args()
   main(args.meshfile,args.esm_id)
