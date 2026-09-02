import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def find_closest_cells_to_river(ds,lon_river_gro,lat_river_gro):
    LON,LAT=np.meshgrid(ds.longitude.data,ds.latitude.data)
    dist=np.sqrt((LON-lon_river_gro)**2+(LAT-lat_river_gro)**2) #compute distance with river for every cell on GloFAS grid

    k=3
    ii_min=np.ones((3,2))
    for i in range(k):
        idx=np.unravel_index(np.argmin(dist),dist.shape)
        ii_min[i]=idx

        dist[idx]=1e9

    return ii_min.astype(int)

    # return ds.sel(longitude=LON[ii_min],latitude=LAT[ii_min]) #returns only cell clostest to river