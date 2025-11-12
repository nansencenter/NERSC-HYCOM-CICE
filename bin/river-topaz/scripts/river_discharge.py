import matplotlib.pyplot as plt
import netCDF4
import matplotlib as mpl
mpl.use('TkAgg')
import numpy as np
import matplotlib.colors as mcolors
mpl.rcParams.update({'font.size': 12}) #Pour changer taille des
import xarray as xr
import abfile
import pandas as pd





#--------------Config-------------------
#To select a region in particular
min_lon = 5
min_lat = 55
max_lon = 30
max_lat = 80

Plot_discharge=True
Plot_uparea=False
Plot_channel=False
Plot_climatology=False




#---------------------Grid file--------------------------
file='grid/regional.grid.a'
grid=abfile.ABFileGrid(file,'r')  #Chargement du ficher grid

# plon=grid.read_field('plon')#Grille de long
# plat=grid.read_field('plat')


# #read_fiel permet de lire un des champs du fichier a
plon=grid.read_field('plon').flatten() #Grille de long
plat=grid.read_field('plat').flatten()

mask=np.where((plon<max_lon) & (plon>min_lon) & (plat<max_lat) & (plat>min_lat))
plon=plon[mask]
plat=plat[mask]

#-----------------Depth file----------------------
#Pour bathymetrie

file_bathy='depth/depth_TP5a0.06_08.a'
depth_grid=abfile.ABFileBathy(file_bathy,'r',idm=800,jdm=760)

depth=depth_grid.read_field('depth').flatten()
depth=depth[mask]


#------------------------River data--------------------------
ds=xr.open_dataarray('river/river.nc')

mask_lon = (ds.longitude >= min_lon) & (ds.longitude <= max_lon)
mask_lat = (ds.latitude >= min_lat) & (ds.latitude <= max_lat)

ds_Norway = ds.where(mask_lon & mask_lat, drop=True)


dis24_March=ds_Norway.values[0,:,:]
dis24_September=ds_Norway.values[1,:,:]
dis24_December=ds_Norway.values[2,:,:]
x=ds_Norway.longitude
y=ds_Norway.latitude
x,y=np.meshgrid(x,y)

##---------------------Upstream Area---------------------------
da_uparea=xr.open_dataarray('river/uparea_glofas_v4_0.nc')

mask_lon = (da_uparea.longitude >= min_lon) & (da_uparea.longitude <= max_lon)
mask_lat = (da_uparea.latitude >= min_lat) & (da_uparea.latitude <= max_lat)

da_uparea_Norway = da_uparea.where(mask_lon & mask_lat, drop=True)


#----------------------Channel path-----------------
da_chan=xr.open_dataset('river/chan_Global_03min.nc')
da_chan=da_chan['Band1']

mask_lon = (da_chan.lon >= min_lon) & (da_chan.lon <= max_lon)
mask_lat = (da_chan.lat >= min_lat) & (da_chan.lat <= max_lat)

da_chan_Norway = da_chan.where(mask_lon & mask_lat, drop=True)


#----------------------Climatology data--------------------------
climatology=pd.read_csv('river/climatology.csv')



climatology_September=climatology.loc[climatology['month']=='September']

climatology_September_Norway=climatology.loc[(climatology['month']=='September') &
                                (climatology['lon']>min_lon) & (climatology['lon']<max_lon) &
                                (climatology['lat']>min_lat) & (climatology['lat']<max_lat)]

#------------------------------Local drainage direction-----------------------------------
ds_ldd=xr.open_dataset('river/ldd_glofas_v4_0.nc')

mask_lon = (ds_ldd.lon >= min_lon) & (ds_ldd.lon <= max_lon)
mask_lat = (ds_ldd.lat >= min_lat) & (ds_ldd.lat <= max_lat)

ds_ldd_Norway = ds_ldd.where(mask_lon & mask_lat, drop=True)


#For the local drainage direction, each number corresponds to a direction (angle).
#We make a dictionnary to link ldd number with corresponding angle
angles={1:225, 2:270, 3: 315, 4: 180, 5: np.nan, 6:0, 7: 135, 8: 90, 9: 45}
#TODO: Direction 5 could be interesting, it means that water doesn't flow anymore

def angles_func(i):
    if np.isnan(i):
        return np.nan
    return angles[i]

LDD_angles=xr.apply_ufunc(angles_func,ds_ldd_Norway['ldd'],vectorize=True)

LON,LAT=np.meshgrid(ds_ldd_Norway.lon,ds_ldd_Norway.lat)

U = np.cos(LDD_angles*np.pi/180)
V = np.sin(LDD_angles*np.pi/180)


#---------------------------------------Plot----------------------------------------------

if Plot_discharge:
    cmap = plt.get_cmap('rainbow')
    cmap.set_under('black') #Color for value less than vmin
    eps =15  #River flux threshold (considered river if above threshold)

    #Create mask for values under eps
    mask=np.ones_like(dis24_September)
    mask[np.where(dis24_September<eps)]=np.nan
    U*=mask
    V*=mask

    fig,ax=plt.subplots()
    plt.tight_layout()
    im=ax.pcolormesh(x,y,dis24_September,vmin=eps,cmap=cmap)
    sc=ax.scatter(plon,plat,s=25,c=depth,cmap='ocean')
    plt.quiver(LON, LAT, U, V, scale=200,angles='xy')
    ax.set_aspect('equal')
    cba=plt.colorbar(im,label='flux '+r'$m^3/s$')
    cbb=plt.colorbar(sc,label='depth [m]')
    plt.legend()
    plt.savefig('fig/Disharge_eps=' + str(eps))
    plt.show()


if Plot_uparea:
    cmap = plt.get_cmap('rainbow')
    cmap.set_under('black')  # Color for value less than vmin
    eps = 1e9

    fig, ax = plt.subplots(figsize=(15, 15))
    im = ax.pcolormesh(da_uparea_Norway.longitude, da_uparea_Norway.latitude, da_uparea_Norway.values, vmin=eps, cmap=cmap)
    sc = ax.scatter(plon, plat, s=25, c=depth, cmap='ocean')
    ax.set_aspect('equal')
    cba = plt.colorbar(im, label='Upstream area ' + r'$m^2$')
    cbb = plt.colorbar(sc, label='depth [m]')
    plt.legend()

    plt.show()

if Plot_channel:
    cmap = mpl.colors.ListedColormap(["navy", "red"])
    norm = mpl.colors.BoundaryNorm(np.arange(-0.5, 2), cmap.N)

    fig, ax = plt.subplots(figsize=(15, 15))
    im = ax.contourf(da_chan_Norway.lon, da_chan_Norway.lat, da_chan_Norway.values,  cmap=cmap)
    sc = ax.scatter(plon, plat, s=25, c=depth, cmap='ocean')
    ax.set_aspect('equal')
    cba = plt.colorbar(im, label='Upstream area ' + r'$m^2$')
    cbb = plt.colorbar(sc, label='depth [m]')
    plt.legend()

    plt.show()

if Plot_climatology:
    cmap = plt.get_cmap('rainbow')
    cmap.set_under('black') #Color for value less than vmin
    eps =25 #River flux threshold (considered river if above threshold)

    fig,ax=plt.subplots()
    im=ax.pcolormesh(x,y,dis24_September,vmin=eps,cmap=cmap)
    sc=ax.scatter(climatology_September_Norway.lon,
                  climatology_September_Norway.lat,s=100,
                  c=climatology_September_Norway.flux,cmap='rainbow')
    ax.set_aspect('equal')
    cba=plt.colorbar(im,label='flux '+r'$m^3/s$')
    cbb=plt.colorbar(sc,label='climatology')
    plt.legend()
    plt.tight_layout()
    plt.show()