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
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap, ListedColormap





#--------------------------------Config-------------------
#To select a region in particular
region_name='Vaasa' #Choose region name

#Dict--> 'Region_name':[min_lon,max_lon,min_lat,max_lat]
Region={'Norway':[5,55,30,80],'Bothnia':[18,28,60,70], 'Vaasa':[21,28,59,66], 'Iceland':[-40,-10,50,80]}
min_lon,max_lon,min_lat,max_lat=Region[region_name]

date_time='2024-07-02T00:00:00.000000000'
date_time_object=datetime.fromisoformat(date_time)#Choose date and time of River discharge data

river_threshold=1 #Threshold to consider flux is related to a river


#----------------------Commands------------------------
Apply_regional_mask=True
Plot_discharge=False
Plot_statistics=True




#---------------------Grid file--------------------------
file='grid/regional.grid.a'
grid=abfile.ABFileGrid(file,'r')  #Chargement du ficher grid

# plon=grid.read_field('plon')#Grille de long
# plat=grid.read_field('plat')


#read_field permet de lire un des champs du fichier a
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
ds_dis=xr.open_dataset('river/river.nc')
ds_dis=ds_dis.sel(valid_time=date_time) #Select date and time

#Create coast: NaN for sea, 1 for land
ds_dis['land']=ds_dis['dis24'].where(np.isnan(ds_dis['dis24']),1)


if Apply_regional_mask:
    mask_lon = (ds_dis.longitude >= min_lon) & (ds_dis.longitude <= max_lon)
    mask_lat = (ds_dis.latitude >= min_lat) & (ds_dis.latitude <= max_lat)
    ds_dis = ds_dis.where(mask_lon & mask_lat, drop=True)


x=ds_dis.longitude
y=ds_dis.latitude

x,y=np.meshgrid(x,y)

#----------------------Climatology data--------------------------
climatology=pd.read_csv('river/climatology.csv')

# climatology=climatology.loc[climatology['month']==date_time_object.strftime('%B')]
climatology=climatology.loc[climatology['month']=='Annual mean']

if Apply_regional_mask:
    climatology=climatology.loc[(climatology['lon']>min_lon) & (climatology['lon']<max_lon) &
                                (climatology['lat']>min_lat) & (climatology['lat']<max_lat)]

#I will assume that the annual mean is in km3/yr, so I will convert it to m3/s
climatology['flux']=climatology['flux']*1e9/31104000

#------------------------------Local drainage direction adn esturary localization-----------------------------------
ds_ldd=xr.open_dataset('river/ldd_glofas_v4_0.nc')


if Apply_regional_mask:
    mask_lon = (ds_ldd.lon >= min_lon) & (ds_ldd.lon <= max_lon)
    mask_lat = (ds_ldd.lat >= min_lat) & (ds_ldd.lat <= max_lat)
    ds_ldd = ds_ldd.where(mask_lon & mask_lat, drop=True)

# We create a mask to have: NaN if the cell can't be an estuary (ldd!=5)
# and 1 if the cell can be an estuary (ldd==5)
mask_estuaries = ds_ldd.where(ds_ldd == 5, 0)
mask_estuaries = mask_estuaries.where(ds_ldd != 5, 1)

# And now we apply mask to Discharge data:
# ds_dis*=mask_estuaries
ds_dis['Estuary_dis'] = ds_dis['dis24'] * mask_estuaries['ldd'].values


def find_estuaries(river_threshold):
    #We transform xarray dataset into a panda Dataframe because otherwise
    # It's impossible to drop NaN along multiple dimensions
    df_dis_Estuary=ds_dis['Estuary_dis'].to_dataframe('dis_Estuary').dropna() #Only keep the points that are estuaries
    df_dis_Estuary=df_dis_Estuary.reset_index() #Just a trick to make the plot easier
    df_dis_Estuary=df_dis_Estuary[df_dis_Estuary['dis_Estuary']>river_threshold] #We ignore the estuaries if the flux is uner the thresold
    return df_dis_Estuary

df_dis_Estuary=find_estuaries(river_threshold)


#For the local drainage direction, each number corresponds to a direction (angle)----> This entire sedtion is just for plotting the arrows.
#We make a dictionnary to link ldd number with corresponding angle
angles={1:225, 2:270, 3: 315, 4: 180, 5: np.nan, 6:0, 7: 135, 8: 90, 9: 45}

def angles_func(i):
    """
    Compute the angle corresponding to the LDD value, and returns NaN for input NaN
    :param i: LDD value
    :return: Corresponding angle (as float otherwise there's a type problem when applying the function)
    """
    if np.isnan(i):
        return np.nan
    return float(angles[i])

# # This can be used to plot the arrows for the LDD
LDD_angles=xr.apply_ufunc(angles_func,ds_ldd['ldd'],vectorize=True) #Vectorize= True because you can't apply a ufun to numpy array
LON,LAT=np.meshgrid(ds_ldd.lon,ds_ldd.lat)
U = np.cos(LDD_angles*np.pi/180)
V = np.sin(LDD_angles*np.pi/180)

mask = np.ones_like(ds_dis['dis24'])
mask[np.where(ds_dis['dis24'] < river_threshold)] = np.nan
U *= mask
V *= mask




#-----------------------------------------Statistics computation---------------------------------------


Nb_estuaries=df_dis_Estuary.dis_Estuary.count() #Number of estuaries
Fulx_sum=df_dis_Estuary.dis_Estuary.sum() #Sum of all fluxes from estuaries




if Plot_discharge:
    land_cmap=cmap = ListedColormap(["black"])

    cmap = plt.get_cmap('rainbow')
    cmap.set_under('black') #Color for value less than vmin
    eps=river_threshold
    flux_min=df_dis_Estuary.dis_Estuary.min()
    flux_max = df_dis_Estuary.dis_Estuary.max()



    fig,ax=plt.subplots()
    plt.tight_layout()
    # im=ax.pcolormesh(x,y,ds_dis['land'],cmap=land_cmap,label='GloFAS grid point')
    im = ax.pcolormesh(x, y, ds_dis['dis24'], vmin=eps,cmap=cmap)

    sc1=ax.scatter(df_dis_Estuary.longitude,df_dis_Estuary.latitude,
                   c=df_dis_Estuary.dis_Estuary,
                   s=150,cmap=cmap,label='GloFAS Estuary',
                   vmin=flux_min,vmax=flux_max,marker='v',edgecolors='black')

    sc2=ax.scatter(plon,plat,s=50,c=depth,cmap='ocean',
                   label='TOPAZ depth grid point')

    plt.quiver(LON, LAT, U, V, scale=200,angles='xy')

    # sc3=ax.scatter(climatology.lon,climatology.lat,s=20,
    #               c=climatology.flux,cmap='rainbow',marker='*',
    #                vmin=flux_min,vmax=flux_max,label='Climatology')

    ax.set_aspect('equal')
    cba=plt.colorbar(sc1,label='flux '+r'$m^3/s$')
    cbb=plt.colorbar(sc2,label='depth [m]')
    # cbc=plt.colorbar(sc3,label='climatology flux '+r'$m^3/s$')
    plt.legend(loc='upper right')
    # plt.savefig('fig/Discharge_eps=' + str(eps))
    plt.show()

if Plot_statistics:
    pass
