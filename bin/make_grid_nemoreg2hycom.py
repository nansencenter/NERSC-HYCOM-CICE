#!/usr/bin/env python
import modeltools.nemo
import argparse
import modeltools.forcing.bathy
import abfile
import netCDF4
import shutil
import numpy.matlib as matlib
import numpy as np
import numpy.ma as ma
import geopy.distance as gpd


# Code reads mesh filename and then generate both regional and bathy [ab] files.
# Note that information about native grid netcdffiles (their paths ,etc) can be found in REGION.src
# I have used hard-coded grid size which you need to change in the case using different subdomains.
#
# Usage:
#
# ../bin/make_grid_nemoreg2hycom.py /cluster/work/users/annettes/TestDownload/Grid/GLO-MFC_001_030_mask_bathy_sub.nc
#
# Mostafa Bakhoday-Paskyabi, 16 May 2019.
# History:
# Annette Samuelsen: modified from NEMO native to NEMO regular , March, 2022
#
#

def make_grid(meshfile):

    ncid0=netCDF4.Dataset(meshfile[:-17]+"coordinates_sub.nc","r")
    np.seterr(invalid='ignore')
    e3t=ncid0.variables["e3t"][:]  #AS: Is not spatially dependent.              
    ncid0.close()
    lev_bnds=np.cumsum(e3t)

    # Now acquire the data. P-cell data
    ncid0=netCDF4.Dataset(meshfile,"r")
    plon1d=ncid0.variables["longitude"][:]
    tmp=np.diff(plon1d)
    tmp=np.append(tmp,tmp[len(tmp)-1])
    ulon1d=plon1d-tmp
    vlon1d=plon1d
    qlon1d=plon1d-tmp

    plat1d=ncid0.variables["latitude"][:]
    tmp=np.diff(plat1d)
    tmp=np.append(tmp,tmp[len(tmp)-1])
    ulat1d=plat1d
    vlat1d=plat1d-tmp
    qlat1d=plat1d-tmp

    plon=matlib.repmat(plon1d,len(plat1d),1)
    plat=np.transpose(matlib.repmat(plat1d,len(plon1d),1))

    hdepw=ncid0.variables["deptho"][:,:]
    hdepw=ma.masked_where(np.isnan(hdepw),hdepw)
    mbathy=ncid0.variables["deptho"][:,:]
    mbathy=ma.masked_where(np.isnan(mbathy),mbathy)
#    tmp=hdepw[:,0]
#    print(np.shape(tmp),tmp.mask)

    jdm,idm=np.shape(hdepw)
    
    ncid0.close()

#C-grid
#--Q--V--Q--
#  |  |  |
#--U--P--U--
#  |  |  |
#--Q--V--Q--
#  |  |  |

    # Compute the grid length:
    scpx=np.zeros((jdm,idm))
    scpy=np.zeros((jdm,idm))
    for i in range(idm-1):
         for j in range(jdm-1):
              scpx[j,i]=gpd.geodesic((plat[j,i],plon[j,i]),(plat[j,i+1],plon[j,i+1])).m
              scpy[j,i]=gpd.geodesic((plat[j,i],plon[j,i]),(plat[j+1,i],plon[j+1,i])).m
    scpx[:,idm-1]=scpx[:,idm-2]
    scpy[:,idm-1]=scpy[:,idm-2]
    scpx[jdm-1,:]=scpx[jdm-2,:]
    scpy[jdm-1,:]=scpy[jdm-2,:]

    # U-cell data.          
    ulon=matlib.repmat(ulon1d,len(ulat1d),1)
    ulat=np.transpose(matlib.repmat(ulat1d,len(ulon1d),1))

    print("PLON",plon[0,0:2],ulon[0,0:2])
    print("PLAT",plat[0,0:2],ulat[0,0:2])

    # Compute the grid length:
    scux=np.zeros((jdm,idm))
    scuy=np.zeros((jdm,idm))               
    for i in range(idm-1):
         for j in range(jdm-1):
              scux[j,i]=gpd.geodesic((ulat[j,i],ulon[j,i]),(ulat[j,i+1],ulon[j,i+1])).m
              scuy[j,i]=gpd.geodesic((ulat[j,i],ulon[j,i]),(ulat[j+1,i],ulon[j+1,i])).m
    scux[:,idm-1]=scpx[:,idm-2]
    scuy[:,idm-1]=scpy[:,idm-2]
    scux[jdm-1,:]=scpx[jdm-2,:]
    scuy[jdm-1,:]=scpy[jdm-2,:]


    # V-cell data.                
    # TODO: Proper extrapolation of data on grid edges (mainly bottom row)                  
    vlon=matlib.repmat(vlon1d,len(vlat1d),1)
    vlat=np.transpose(matlib.repmat(vlat1d,len(vlon1d),1))

    # Compute the grid length:               
    scvx=np.zeros((jdm,idm))
    scvy=np.zeros((jdm,idm))
    for i in range(idm-1):
         for j in range(jdm-1):
              scvx[j,i]=gpd.geodesic((vlat[j,i],vlon[j,i]),(vlat[j,i+1],vlon[j,i+1])).m
              scvy[j,i]=gpd.geodesic((vlat[j,i],vlon[j,i]),(vlat[j+1,i],vlon[j+1,i])).m
    scvx[:,idm-1]=scpx[:,idm-2]
    scvy[:,idm-1]=scpy[:,idm-2]
    scvx[jdm-1,:]=scpx[jdm-2,:]
    scvy[jdm-1,:]=scpy[jdm-2,:]

    # Q-cell data         
    # TODO: Proper extrapolation of data on grid edges (mainly bottom row)       
    qlon=matlib.repmat(qlon1d,len(qlat1d),1)
    qlat=np.transpose(matlib.repmat(qlat1d,len(qlon1d),1))

    # Compute the grid length:               
    scqx=np.zeros((jdm,idm))
    scqy=np.zeros((jdm,idm))
    for i in range(idm-1):
         for j in range(jdm-1):
              scqx[j,i]=gpd.geodesic((qlat[j,i],qlon[j,i]),(qlat[j,i+1],qlon[j,i+1])).m
              scqy[j,i]=gpd.geodesic((qlat[j,i],qlon[j,i]),(qlat[j+1,i],qlon[j+1,i])).m
    scqx[:,idm-1]=scpx[:,idm-2]
    scqy[:,idm-1]=scpy[:,idm-2]
    scqx[jdm-1,:]=scpx[jdm-2,:]
    scqy[jdm-1,:]=scpy[jdm-2,:]

    print("MINX",np.min(scqx),"MAXX",np.max(scqx))
    print("MINY",np.min(scqy),"MAXY",np.max(scqy))

    # Angle used for rotation                 
#    ulon_rgt = np.copy(ulon)
#    ulat_rgt = np.copy(ulat)
#    ulon_lft = np.roll(ulon,1,axis=1)
#    ulat_lft = np.roll(ulat,1,axis=1)
    pang = modeltools.tools.p_azimuth(ulon,ulat,plon,plat)

    # Aspect ratio           
    asp = np.where(scpy==0.,99.0,scpx/scpy)

    # Coriolis               
    corio = np.sin(np.radians(qlat)) * 4. * np.pi / 86164.0 # Sidereal day              

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
    tmp2 = np.where( mbathy>0.  ,hdepw,0.)
    tmp2 = np.where( mbathy>1.0e19,0.0 ,tmp2)
    abfile.write_bathymetry("bathy",1,tmp2,0.)
#    hdepw=np.where(np.isnan(hdepw),-1.0,hdepw)
#    abfile.write_bathymetry("bathy",1,hdepw,0.)
    
    shutil.move("depth_bathy_01.a",'./dummped_depth_NMOc0.00_01.a')
    shutil.move("depth_bathy_01.b",'./dummped_depth_NMOc0.00_01.b')
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
    
    make_grid(args.meshfile)


