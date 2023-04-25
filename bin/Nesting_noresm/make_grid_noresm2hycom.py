#!/usr/bin/env python
import modeltools.nemo
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import modeltools.forcing.bathy
import abfile
import numpy
import netCDF4
import logging
import re
import pyproj
import os.path
import time
import cfunits
import sys
import shutil
import glob
import numpy as np
import geopy.distance as gpd


# Code reads mesh filename and then generate both regional and bathy [ab] files.
# Note that information about native grid netcdffiles (their paths ,etc) can be found in REGION.src
# I have used hard-coded grid size which you need to change in the case using different subdomains.
#
# Usage:
#
# ../bin/Nesting_noresm/make_grid_noresm2hycom.py /nird/projects/nird/NS9481K/annettes/O2OCEAN/grid/areacello_Ofx_NorESM2-MM_historical_r1i1p1f1_gn.nc
#
# Mostafa Bakhoday-Paskyabi, 16 May 2019.
# History:
# Annette Samuelsen: modified from NEMO to NORESM, April, 2022
#
#


def u_to_hycom_u(field2d)  :
   return numpy.roll(field2d,1,axis=1)


# nemo V at i,j -> HYCOM U at i,j+1. 
def v_to_hycom_v(field2d,extrapolate="none")  :
   myfield=numpy.copy(field2d)
   # For now the bottom row is just replicated
   myfield[1:,:] = myfield[:-1,:]
   
   return myfield

def f_to_hycom_q(field,extrapolate="none")  :
   myfield=u_to_hycom_u(field)
   myfield=v_to_hycom_v(myfield)
   return myfield



def make_grid(meshfile ,idm=360,jdm=385):
    
    ncidg=netCDF4.Dataset(meshfile,"r")
    # AS: since the data are on a regular grid only plon and plat are needed.
    # Now acquire the data. P-cell data
    plon = np.squeeze(ncidg.variables["longitude"][:,:])
    plat = np.squeeze(ncidg.variables["latitude"][:,:])
    vertices_plon = np.squeeze(ncidg.variables["vertices_longitude"][:,:,:])
    vertices_plat = np.squeeze(ncidg.variables["vertices_latitude"][:,:,:])
    ncidg.close()
    print("PLON",plon[0,0:2],vertices_plon[0,0,:])
    print("PLAT",plat[0,0:2],vertices_plat[0,0,:])

    # Compute the grid length:
    scpx=np.zeros((jdm,idm))
    scpy=np.zeros((jdm,idm))
    for i in range(idm-1):
         for j in range(jdm-1):
              scpx[j,i]=gpd.geodesic((plat[j,i],plon[j,i]),(plat[j,i+1],plon[j,i+1])).m
              scpy[j,i]=gpd.geodesic((plat[j,i],plon[j,i]),(plat[j+1,i],plon[j+1,i])).m

    # U-cell data.                                                                                                                                 
    ulon = np.squeeze(vertices_plon[:,:,0])
    ulat = np.squeeze(plat)

    # Compute the grid length:
    scux=np.zeros((jdm,idm))
    scuy=np.zeros((jdm,idm))                                                                                                                     
    for i in range(idm-1):
         for j in range(jdm-1):
              scux[j,i]=gpd.geodesic((ulat[j,i],ulon[j,i]),(ulat[j,i+1],ulon[j,i+1])).m
              scuy[j,i]=gpd.geodesic((ulat[j,i],ulon[j,i]),(ulat[j+1,i],ulon[j+1,i])).m

    # V-cell data.                                                                                                                      
    # TODO: Proper extrapolation of data on grid edges (mainly bottom row)                                                    
    vlon = np.squeeze(plon)
    vlat = np.squeeze(vertices_plat[:,:,2])

    # Compute the grid length:                                                                                                                     
    scvx=np.zeros((jdm,idm))
    scvy=np.zeros((jdm,idm))
    for i in range(idm-1):
         for j in range(jdm-1):
              scvx[j,i]=gpd.geodesic((vlat[j,i],vlon[j,i]),(vlat[j,i+1],vlon[j,i+1])).m
              scvy[j,i]=gpd.geodesic((vlat[j,i],vlon[j,i]),(vlat[j+1,i],vlon[j+1,i])).m

    # Q-cell data                                                                                              
    # TODO: Proper extrapolation of data on grid edges (mainly bottom row)                                         
    qlon = np.squeeze(np.squeeze(vertices_plon[:,:,2]))
    qlat = np.squeeze(np.squeeze(vertices_plat[:,:,2]))
    # Compute the grid length:                                                                                                                     
    scqx=np.zeros((jdm,idm))
    scqy=np.zeros((jdm,idm))
    for i in range(idm-1):
         for j in range(jdm-1):
              scqx[j,i]=gpd.geodesic((qlat[j,i],qlon[j,i]),(qlat[j,i+1],qlon[j,i+1])).m
              scqy[j,i]=gpd.geodesic((qlat[j,i],qlon[j,i]),(qlat[j+1,i],qlon[j+1,i])).m

    print("MINX",np.min(scqx),"MAXX",np.max(scqx))
    print("MINY",np.min(scqy),"MAXY",np.max(scqy))

    # Angle used for rotation                                                                                                                       
    ulon_rgt = numpy.copy(ulon)
    ulat_rgt = numpy.copy(ulat)
    ulon_lft = numpy.roll(ulon,1,axis=1)
    ulat_lft = numpy.roll(ulat,1,axis=1)
    pang = modeltools.tools.p_azimuth(ulon_lft,ulat_lft,ulon_rgt,ulat_rgt)

    # Aspect ratio                                                                                                                                  
    asp = numpy.where(scpy==0.,99.0,scpx/scpy)

    # Coriolis                                                                                                                                      
    corio = numpy.sin(numpy.radians(qlat)) * 4. * numpy.pi / 86164.0 # Sidereal day                                                                 

    # Put inside datadict for abfile writing                                                                                                        
    ddict={}
    ddict["plon"]=plon
    ddict["plat"]=plat
    ddict["ulon"]=ulon
    ddict["ulat"]=ulat
    ddict["vlon"]=vlon
    ddict["vlat"]=vlat
    ddict["qlon"]=qlon
    ddict["qlat"]=qlat
    #                                                                                                                                               
    ddict["scpx"]=scpx
    ddict["scpy"]=scpy
    ddict["scux"]=scux
    ddict["scuy"]=scuy
    ddict["scvx"]=scvx
    ddict["scvy"]=scvy
    ddict["scqx"]=scqx
    ddict["scqy"]=scqy
    #                                                                                                                                               
    ddict["cori"]=corio
    ddict["pang"] =pang
    ddict["pasp"] =asp
    abfile.write_regional_grid(ddict)

    # Bathymetry
    NORESM_depth="/nird/projects/nird/NS9481K/annettes/O2OCEAN/grid/deptho_Ofx_NorESM2-MM_historical_r3i1p1f1_gn.nc"
    dncid=netCDF4.Dataset(NORESM_depth,"r")
    hdepw  = np.squeeze(dncid.variables["deptho"][:,:])
    mbathy = np.squeeze(dncid.variables["deptho"][:,:])
    dncid.close()

#    tmask = np.squeeze(ncidg.variables["pmask"][:,:])
#    print(hdepw.mask[0,:])
#    print(hdepw.mask[:,0])
#    print(hdepw[0,:])
#    print(hdepw[:,0])
    tmp2 = numpy.where( mbathy>0.  ,hdepw,0.)
    tmp2 = numpy.where( mbathy>1.0e19,0.0 ,tmp2)
    abfile.write_bathymetry("bathy",1,tmp2,0.)
    
    shutil.move("depth_bathy_01.a",'./dummped_depth_ESMa1.00_01.a')
    shutil.move("depth_bathy_01.b",'./dummped_depth_ESMa1.00_01.b')
    shutil.move("regional.grid.a","./dummped_regional.grid.a")
    shutil.move("regional.grid.b","./dummped_regional.grid.b")
    
if __name__ == "__main__" :

    parser = argparse.ArgumentParser(
          description='This tool will convert NORESM netcdf files to hycom archive files. It will also create grid and topo files for hycom.'
          )
    parser.add_argument('meshfile',   type=str,help="NEMO mesh file in netcdf format")
    parser.add_argument('--idm',    type=int,default=360,  help="    ")
    parser.add_argument('--jdm',    type=int,default=385,  help="    ")

    args = parser.parse_args()
    
    make_grid(args.meshfile,idm=args.idm,jdm=args.jdm )
