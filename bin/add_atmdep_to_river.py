import numpy as np
import abfile.abfile as abf
import argparse
from netCDF4 import Dataset as NetCDFFile
from scipy.interpolate import griddata
import os
'''
CAGLAR - November 2020

The script adds additional N and P on the existing river N and P forcing files.
Therefore use this script only if you have those in your force/rivers/expt_no/ folder
This script can be used standalone and need only be used once, or can be used in create_ref_case.sh
It will keep a copy of the original river files in the same directory.

how to use:
python add_atmdep_to_river.py PATH-TO-RIVER-FORCING-FOLDER  ATMOSPHERIC-DEPOSITIONS-NETCDF-FILE

python add_atmdep_to_river.py /cluster/work/users/cagyum/TP0a1.00/force/rivers/022/  /cluster/home/cagyum/NERSC-HYCOM-CICE/input/emep_2010_annual_1degree_rv4_17gfecl1p0.nc
'''
def main(path,source_file):
    print(path)
    print(source_file)
    nc = NetCDFFile(source_file,"r")
    wetNOy = nc.variables["WDEP_NOy"][:,:]
    wetNHx = nc.variables["WDEP_NHx"][:,:]
    dryNOy = nc.variables["DDEP_NOy_m2Grid"][:,:]
    dryNHx = nc.variables["DDEP_NHx_m2Grid"][:,:]
    lat = nc.variables["lat"][:]
    lon = nc.variables["lon"][:]
    Nread = wetNOy + wetNHx + dryNOy + dryNHx

    abgrid = abf.ABFileGrid(path + "../../../topo/regional.grid","r")
    plon=abgrid.read_field("plon")
    plat=abgrid.read_field("plat")
    jdm,idm=plon.shape

    xi,yi = np.meshgrid(lon,lat)
    N = griddata((xi.flatten(),yi.flatten()),Nread.flatten(),(plon,plat),method='linear') 
    N[np.isnan(N)]=0.
    N = N/365./86400. # year --> second
    N = N/14.01 * 6.625 * 12.01 # mgN m-2 s-1 --> mgC m-2 s-1


    outfile=abf.ABFileRiver(path + "ECO_no3_new.a","w",idm=idm,jdm=jdm,\
                   cline1='River nitrate fluxes + Atmospheric N deposition',\
                   cline2='mgC m-2 s-1')
    outfile.write_header()
    Nriver = abf.AFile(idm,jdm,path + "ECO_no3.a","r")

    for month in range(12):
        river  = Nriver.read_record(month)
        outfile.write_field(river+N,None,"river nitrate",month+1)

    Nriver.close()
    outfile.close()

    origAfile = path + "ECO_no3.a"
    origBfile = path + "ECO_no3.b"
    oldAfile  = path + "ECO_no3_noATM.a"
    oldBfile  = path + "ECO_no3_noATM.b"
    newAfile  = path + "ECO_no3_new.a"
    newBfile  = path + "ECO_no3_new.b"

    os.rename(origAfile,oldAfile) # river without atmpspheric deposition is kept as noATM file.
    os.rename(origBfile,oldBfile) # river without atmpspheric deposition is kept as noATM file.
    os.rename(newAfile,origAfile) # river with atmospheric deposition is renamed to the original file used by hycom
    os.rename(newBfile,origBfile) # river with atmospheric deposition is renamed to the original file used by hycom

    # according to Okin et al, 2011 - doi:10.1029/2010GB003858
    # for North Atlantic (northern section)
    # Ndep = 7.4 TgN/yr
    # Pdep = 0.02 TgP/yr
    # ratio --> ( 7.4 / 14.01 * 6.625 ) / ( 0.02 / 31. * 106. ) = 51.2
    P = N/51.2

    outfile=abf.ABFileRiver(path + "ECO_pho_new.a","w",idm=idm,jdm=jdm,\
                   cline1='River phosphate fluxes + Atmospheric P deposition',\
                   cline2='mgC m-2 s-1')
    outfile.write_header()
    Priver = abf.AFile(idm,jdm,path + "ECO_pho.a","r")

    for month in range(12):
        river  = Priver.read_record(month)
        outfile.write_field(river+P,None,"river phosphate",month+1)

    Priver.close()
    outfile.close()

    origAfile = path + "ECO_pho.a"
    origBfile = path + "ECO_pho.b"
    oldAfile  = path + "ECO_pho_noATM.a"
    oldBfile  = path + "ECO_pho_noATM.b"
    newAfile  = path + "ECO_pho_new.a"
    newBfile  = path + "ECO_pho_new.b"

    os.rename(origAfile,oldAfile) # river without atmpspheric deposition is kept as noATM file.
    os.rename(origBfile,oldBfile) # river without atmpspheric deposition is kept as noATM file.
    os.rename(newAfile,origAfile) # river with atmospheric deposition is renamed to the original file used by hycom
    os.rename(newBfile,origBfile) # river with atmospheric deposition is renamed to the original file used by hycom

    '''
    #in case you want to inspect the results
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    from matplotlib import ticker

    m = Basemap(width=8500000,height=8500000,
    resolution='c',projection='stere',\
    lat_ts=80,lat_0=80,lon_0=-40.,round='True')
    x,y=m(plon,plat)
    x1,y1 = m(xi,yi)

    cmin = 0.; cmax = 5E-4    
    warnings.filterwarnings('ignore')
    cmap = plt.get_cmap('Spectral_r')
    figure=plt.figure(figsize=(24,6))
    ax=figure.add_subplot(411)
    ax.set_position([0.0,0.025,0.25,0.95])

    levels = MaxNLocator(nbins=15).tick_values(cmin, cmax)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    pmesh = m.pcolormesh(x,y,N,cmap=cmap,norm=norm)

    m.drawcoastlines(linewidth=0.25)
    m.fillcontinents(color='lightgrey')
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='whitesmoke')

    ax2=figure.add_subplot(412)
    ax2.set_position([0.25,0.025,0.25,0.95])

    pmesh2 = m.pcolormesh(x1,y1,Nread/365./86400./14.01*6.625*12.01,cmap=cmap,norm=norm)
    m.drawcoastlines(linewidth=0.25)
    m.fillcontinents(color='lightgrey')
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='whitesmoke')

    ax3=figure.add_subplot(413)
    ax3.set_position([0.5,0.025,0.25,0.95])

    pmesh3 = m.pcolormesh(x,y,river,cmap=cmap,norm=norm)
    m.drawcoastlines(linewidth=0.25)
    m.fillcontinents(color='lightgrey')
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='whitesmoke')

    ax4=figure.add_subplot(414)
    ax4.set_position([0.75,0.025,0.25,0.95])

    pmesh4 = m.pcolormesh(x,y,river+N,cmap=cmap,norm=norm)
    m.drawcoastlines(linewidth=0.25)
    m.fillcontinents(color='lightgrey')
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='whitesmoke')


    cbaxes = figure.add_axes([0.95, 0.025, 0.01, 0.175])
    cb = ax.figure.colorbar(pmesh, cax = cbaxes)

    '''

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('path',  type=str, help="path to river forcing files")
    parser.add_argument('source_file',  type=str, help="N deposition file")
    args = parser.parse_args()

    main(args.path,args.source_file)

