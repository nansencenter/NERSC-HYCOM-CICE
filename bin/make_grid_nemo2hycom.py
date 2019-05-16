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


# Code reads mesh filename and then generate both regional and bathy [ab] files.
# Note that information about native grid netcdffiles (their paths ,etc) can be found in REGION.src
# I have used hard-coded grid size which you need to change in the case using different subdomains.
# Mostafa Bakhoday-Paskyabi, 16 May 2019.
# History:
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



def make_grid(filemesh ,idm=4320,jdm=1188):
    
    meshfile_hgr=filemesh[:-6]+"hgr.nc"
    meshfile_zgr=filemesh[:-6]+"zgr.nc"
    maskfile=filemesh[:-11]+"mask.nc"

    ncidh=netCDF4.Dataset(meshfile_hgr,"r")
    ncidz=netCDF4.Dataset(meshfile_zgr,"r")
    ncidm=netCDF4.Dataset(maskfile,"r")

    # Now acquire the data. P-cell data
    plon = np.squeeze(ncidh.variables["glamt"][0,:,:])
    plat = np.squeeze(ncidh.variables["gphit"][0,:,:])
    scpx = np.squeeze(ncidh.variables["e1t"]  [0,:,:])
    scpy = np.squeeze(ncidh.variables["e2t"]  [0,:,:])

    # U-cell data. 
    ulon = np.squeeze(u_to_hycom_u( ncidh.variables["glamu"][0,:,:]))
    ulat = np.squeeze(u_to_hycom_u( ncidh.variables["gphiu"][0,:,:]))
    scux = np.squeeze(u_to_hycom_u( ncidh.variables["e1u"  ][0,:,:]))
    scuy = np.squeeze(u_to_hycom_u( ncidh.variables["e2u"  ][0,:,:]))

    # V-cell data.
    # TODO: Proper extrapolation of data on grid edges (mainly bottom row)
    vlon = np.squeeze(v_to_hycom_v( ncidh.variables["glamv"][0,:,:]))
    vlat = np.squeeze(v_to_hycom_v( ncidh.variables["gphiu"][0,:,:]))
    scvx = np.squeeze(v_to_hycom_v( ncidh.variables["e1v"]  [0,:,:]))
    scvy = np.squeeze(v_to_hycom_v( ncidh.variables["e2v"]  [0,:,:]))

    # Q-cell data
    # TODO: Proper extrapolation of data on grid edges (mainly bottom row)
    qlon = np.squeeze(f_to_hycom_q( ncidh.variables["glamf"][0,:,:] ))
    qlat = np.squeeze(f_to_hycom_q( ncidh.variables["gphif"][0,:,:] ))
    scqx = np.squeeze(f_to_hycom_q( ncidh.variables["e1f"  ][0,:,:] ))
    scqy = np.squeeze(f_to_hycom_q( ncidh.variables["e2f"  ][0,:,:] ))

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
    hdepw  = np.squeeze(ncidz.variables["hdepw"] [0,:,:])
    mbathy = np.squeeze(ncidz.variables["mbathy"][0,:,:])

    ncidt=netCDF4.Dataset(maskfile,"r")
    tmask = np.squeeze(ncidt.variables["tmaskutil"][0,:,:])
    tmp2 = numpy.where( mbathy>0.  ,hdepw,0.)
    abfile.write_bathymetry("bathy",1,tmp2,0.)
    
    shutil.move("depth_bathy_01.a",'./dummped_depth_NMOa0.08_01.a')
    shutil.move("depth_bathy_01.b",'./dummped_depth_NMOa0.08_01.b')
    shutil.move("regional.grid.a","./dummped_regional.grid.a")
    shutil.move("regional.grid.b","./dummped_regional.grid.b")
    
if __name__ == "__main__" :

    parser = argparse.ArgumentParser(
          description='This tool will convert NEMO netcdf files to hycom archive files. It will also create grid and topo files for hycom.'
          )
    parser.add_argument('meshfile',   type=str,help="NEMO mesh file in netcdf format")
    parser.add_argument('--idm',    type=int,default=4320,  help="    ")
    parser.add_argument('--jdm',    type=int,default=1188,  help="    ")

    args = parser.parse_args()
    
    make_grid(args.meshfile,idm=args.idm,jdm=args.jdm )
