import os
import sys
import time
import xarray as xr
from utils.grid_operations import refine_grid
from utils.grid_operations import transform_coordinates
import numpy as np
import pandas as pd
# import imageio.v3 as iio
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.use('TkAgg')
import warnings
import scipy.sparse as sp

# Shared with the AHYPE v4b river forcing: bin/river_spreading.py and bin/geus.py
_BINDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _BINDIR not in sys.path:
    sys.path.insert(0, _BINDIR)
from river_spreading import spreading_weights

class GloFAS_grid:
    """ Class containing the GloFAS data grid as a Xarray dataset """
    def __init__(self,file,region_name,region_list,field_name='dis24',x_coord='longitude',y_coord='latitude',
                 static=False,climatology=False, month=None,GloFAS_lon_shape=7200,GloFAS_lat_shape=3000,lazy_mode=False):
        """
        Constructor for the GloFAS_grid class. Opens the NetDCF4 daily discharge file with data on GloFAS grid and puts the data in a Xarray dataset.
        This class is used for daily discharge GloFAS data but also for auxiliary files: LDD, upstream area and elevation that are on the GloFAS grid.
        :param file: GloFAS NETCDF4 daily discharge file name [str]
        :param region_name: Name of the region of interest [str]
        :param region_list: Lists of available regions [dict]
        :param date: Date of the discharge data [str]
        :param field_name: Name of the field containing discharge data in NETCDF4 file [str]
        :param x_coord: Name of the x coordinate field in NETCDF4 file [str]
        :param y_coord: Name of the y coordinate field in NETCDF4 file [str]
        :param static: Whether the file a static auxiliary file (True) or a daily discharge file (False) [bool]
        """
        # self.ds = xr.open_dataset(file)
        self.ds=file



        self.__region_name=region_name
        self.__Region = region_list


        if climatology==False:
            try:
                self.__date=self.ds.valid_time.values[0]
            except:
                self.__date = self.ds.valid_time.values
            # self.__date_time_object=self.ds.valid_time.values[0]

            if not static:
                try:
                    self.ds = self.ds.sel(valid_time=self.ds.valid_time.values[
                    0])  # Select date and time, we take the first date in the datatset, normally there's only one, so that doesn't change much, it's mainly a shape problem
                except: pass

                # recreate the dataset to ensure right dim if the dataset isn't global (not all lon or lat)
                lon, lat = np.zeros(GloFAS_lon_shape), np.zeros(GloFAS_lat_shape)
                dis24 = np.zeros((GloFAS_lat_shape, GloFAS_lon_shape))

                lon[:self.ds.longitude.shape[0]] = self.ds.longitude
                lat[:self.ds.latitude.shape[0]] = self.ds.latitude
                # print(self.ds)
                # print('data:::::',self.ds.dis24)
                dis24[:self.ds.latitude.shape[0], :self.ds.longitude.shape[0]] = self.ds.dis24.data

                ds_tmp = xr.Dataset(coords={'longitude': lon, 'latitude': lat,
                                            'valid_time':self.__date, 'surface':0})
                ds_tmp['dis24'] = (['latitude', 'longitude'], dis24)

                self.ds=ds_tmp


                self.__lon = self.ds.longitude
                self.__lat = self.ds.latitude

        else:
            self.__date = month

            if not static:
                self.ds = self.ds.sel(valid_time=month) # Selecting right month
                # self.ds['dis24'] = self.ds['dis24'][0, :, :]  # dropping time dim on discharge values

                self.__lon = self.ds.longitude
                self.__lat = self.ds.latitude



        self.__field_name=field_name #The name of the field containing discharge values, that can be useful

        self.__xcoord=x_coord #Depending on wether they call it 'longitude' or 'lon' in the file
        self.__ycoord=y_coord


        self.apply_regional_mask()
        self.create_coast()

    def apply_regional_mask(self):
        """
        Applies the regional mask to the dataset by dropping with that are outside the lon/lat limits
        """
        min_lon,min_lat,max_lon,max_lat = self.__Region[self.__region_name]
        mask_lon = (self.ds[self.__xcoord] >= min_lon) & (self.ds[self.__xcoord] <= max_lon)
        mask_lat = (self.ds[self.__ycoord] >= min_lat) & (self.ds[self.__ycoord] <= max_lat)

        self.ds = self.ds.where(mask_lon & mask_lat, drop=True)

    def create_coast(self):
        # Create coast: NaN for sea, 1 for land
        self.ds['land'] = self.ds[self.__field_name].where(np.isnan(self.ds[self.__field_name]), 1)


    def add_variable(self,file,field_name,static=False,x_coord='longitude',y_coord='latitude'):
        """
        Add a variable to the dataset, this will be used to add LDD, uparea and elevation to the GloFAS grid dataset.
        :param file: Name of the NETCDF4 file to add to the dataset [str]
        :param field_name: Name of the field in NETCDF4 file [str]
        :param static: Whether the file is auxiliary static grid or not. Here we only add static grids, so it will always be True [bool]
        :param x_coord: Name of the x coordinate field in NETCDF4 file [str]
        :param y_coord: Name of the y coordinate field in NETCDF4 file [str]
        """
        var=xr.open_dataset(file) #Loading the new varible

        if x_coord=='lon':
            var=var.rename({'lon':'longitude','lat':'latitude'})
            x_coord,y_coord='longitude','latitude'


        var[x_coord]=self._GloFAS_grid__lon #Slightly extreme way to realign the datasets, we can probably do better
        var[y_coord]=self._GloFAS_grid__lat

        min_lon,min_lat,max_lon,max_lat = self._GloFAS_grid__Region[self._GloFAS_grid__region_name]
        mask_lon = (var[x_coord] >= min_lon) & (var[x_coord] <= max_lon) #Applying the regional mask, we can't use the one from the GloFAS grid because sometimes x_coord is longitude and sometimes lon
        mask_lat = (var[y_coord] >= min_lat) & (var[y_coord] <= max_lat)
        var = var.where(mask_lon & mask_lat, drop=True)

        self.ds[field_name]=(['latitude','longitude'],var[field_name].data) #We create a new field in the dataset

    # Approximate positions (lat, lon) of the rivers corrected in correct_largest_russian_rivers()
    LARGE_RIVERS = {'Lena': (72, 129), 'Yenisei': (70.90, 83.30), 'Ob': (66.5, 69.8)}

    def find_estuaries(self,lazy_mode=False,input_grids=None):
        """
        Finds the estuaries: all GloFAS river mouths (LDD = 5) that are close enough to the model grid
        (finite distance_to_topaz, see get_glofas_into_topaz). Every GloFAS catchment ends in exactly one
        river mouth, so each river is counted once.

        The result, self.df_dis_Estuary, has one row per river mouth. merge_estuaries() then combines the
        river mouths that are moved to the same model cell; correct_largest_russian_rivers() must be called
        before that.
        """

        if lazy_mode==False:
            mask_estuaries=np.where(self.ds['ldd']== 5 , 1,np.nan) #Mask based on LDD
            mask_estuaries=mask_estuaries*np.where(np.isfinite(self.ds['distance_to_topaz']), 1,np.nan) #Drop estuaries too far from the model grid (NaN distance, see get_glofas_into_topaz)

            self.mask_estuaries=mask_estuaries #We put the mask in GloGAS grid object so that it can be re used

        elif lazy_mode==True: #In lazy mode, we don't compute mask
            mask_estuaries=input_grids.dis_grid.mask_estuaries

        # And now we apply mask to Discharge data:
        self.ds['dis_Estuary'] = self.ds[self._GloFAS_grid__field_name] * mask_estuaries

        # We transform xarray dataset into a panda Dataframe because otherwise
        # It's impossible to drop NaN along multiple dimensions

        self.df_dis_Estuary = self.ds.to_dataframe().dropna(subset='dis_Estuary') # Only keep the points that are estuaries

        self.df_dis_Estuary = self.df_dis_Estuary.reset_index()  # Just a trick to make life easier when locating an index

        self.df_dis_Estuary = self.df_dis_Estuary.round({'longitude': 5, 'latitude': 5}) #rounding up the edges

        # River mouths of the rivers corrected in correct_largest_russian_rivers(): the mouth with the largest
        # upstream area within 2 degrees of the approximate position. Found on the first day (the upstream
        # area is not loaded in lazy mode) and reused.
        if lazy_mode==False:
            self.large_river_mouths = {}
            for river, (lat, lon) in self.LARGE_RIVERS.items():
                df = self.df_dis_Estuary
                near = np.hypot(df['latitude'] - lat, df['longitude'] - lon) < 2
                mouth = df.loc[df.loc[near, 'uparea'].idxmax()]
                self.large_river_mouths[river] = (mouth['latitude'], mouth['longitude'])
                print('{river} river mouth: {lat:.3f}, {lon:.3f}, upstream area {up:.0f} km2'.format(
                    river=river, lat=mouth['latitude'], lon=mouth['longitude'], up=mouth['uparea'] / 1e6))
        else:
            self.large_river_mouths = input_grids.dis_grid.large_river_mouths

    def add_geus_gates(self,geus,depth_grid,distance_threshold=300,lazy_mode=False,input_grids=None):
        """
        Adds the GEUS Greenland ice discharge gates (see geus.py) to self.df_dis_Estuary as extra river
        mouths, at the nearest ocean cell of the model grid, with the discharge of the current month.
        They are merged and spread like the rivers. Called after correct_largest_russian_rivers() and
        before merge_estuaries().
        """
        from scipy.spatial import cKDTree
        from geus import geus_discharge

        if lazy_mode==False:
            ocean=~np.asarray(depth_grid.depth.mask).ravel()
            ocean_idx=np.flatnonzero(ocean)
            lon=np.asarray(depth_grid.lon).ravel()[ocean]
            lat=np.asarray(depth_grid.lat).ravel()[ocean]
            tree=cKDTree(transform_coordinates(np.column_stack([lat, lon])))
            dist,k=tree.query(transform_coordinates(np.column_stack([geus['lat'], geus['lon']])))  # km
            keep=dist<=distance_threshold
            self.geus_gates=pd.DataFrame({'gate':np.flatnonzero(keep),
                                          'latitude':geus['lat'][keep],'longitude':geus['lon'][keep],
                                          'idx_glofas_on_topaz':ocean_idx[k[keep]],
                                          'latitude_glofas_on_topaz':lat[k[keep]],
                                          'longitude_glofas_on_topaz':lon[k[keep]],
                                          'distance_to_topaz':dist[keep]})
            print('GEUS: {n} of {N} gates within {d} km of the model grid'.format(
                n=keep.sum(),N=keep.size,d=distance_threshold))
        else:
            self.geus_gates=input_grids.dis_grid.geus_gates

        date=self.df_dis_Estuary['valid_time'].iloc[0]
        q=geus_discharge(geus,date)[self.geus_gates['gate'].values]
        gates=self.geus_gates.drop(columns='gate').assign(
            dis24=q.astype(np.float32),dis_Estuary=q,valid_time=date,
            surface=self.df_dis_Estuary['surface'].iloc[0],land=1.0,ldd=5.0,uparea=np.nan,elevation=np.nan,
            longitude_grid=self.geus_gates['longitude'],latitude_grid=self.geus_gates['latitude'])
        self.df_dis_Estuary=pd.concat([self.df_dis_Estuary,gates],ignore_index=True)

    def merge_estuaries(self,lazy_mode=False):
        """
        Combines the river mouths that are moved to the same model cell (same idx_glofas_on_topaz) into
        one estuary, summing their discharge. Called after correct_largest_russian_rivers().
        """
        #Careful here with the field names !! (some may be missing, change the aggregate functions if needed)
        aggregation_functions = {'dis24':'sum', 'valid_time':'first', 'surface':'first', 'land':'first', 'ldd':'first', 'uparea':'sum', 'elevation':'first',
       'longitude_grid':'first', 'latitude_grid':'first', 'latitude_glofas_on_topaz':'first',
       'longitude_glofas_on_topaz':'first', 'distance_to_topaz':'first',
       'dis_Estuary':'sum','latitude':'first','longitude':'first'}

        # In lazy mode, we dont have LDD, uparea, etc so we change agg function slightly
        lazy_aggregation_functions={'dis24':'sum', 'valid_time':'first', 'surface':'first', 'land':'first',
       'longitude_grid':'first', 'latitude_grid':'first', 'latitude_glofas_on_topaz':'first',
       'longitude_glofas_on_topaz':'first', 'distance_to_topaz':'first',
       'dis_Estuary':'sum','latitude':'first','longitude':'first'}

        if lazy_mode==False:
            self.df_dis_Estuary=self.df_dis_Estuary.groupby('idx_glofas_on_topaz',as_index=False).aggregate(aggregation_functions)
        else:
            self.df_dis_Estuary=self.df_dis_Estuary.groupby('idx_glofas_on_topaz',as_index=False).aggregate(lazy_aggregation_functions)



    def write_GloFAS_clim(self, output_file_name='GloFAS_clim.dat') :
        """
        Writes a txt files in the form of AHYPE climatology with the GloFAS estuaries
        :param output_file_name: Name of the output txt file [str]
        """
        df=self.df_dis_Estuary
        with open(output_file_name, 'w') as clim_file:
            for i in df.index.values:
                dis_km3_per_month=df.loc[[i]].dis_Estuary.values[0]/1e9*3600*24*30
                clim_file.write(f"T  NoName    GloFAS daily climatology\n")
                clim_file.write(f"{dis_km3_per_month:12.5f}") #Average

                #TODO: Monthly weights, altough it doesn't make a lot of sens here
                for j in range(12): clim_file.write(f"{1:8.5f}")

                clim_file.write("\n")
                clim_file.write(f"{df.loc[[i]].latitude.values[0]:9.2f}")
                clim_file.write(f"{df.loc[[i]].longitude.values[0]:9.2f}\n")

    def get_glofas_into_topaz(self,depth_grid,distance_threshold=300,lazy_mode=False,input_grids=None):
        """
        Interpolates the GloFAS point (originally on the GloFAS grid) onto the TOPAZ depth grid using a nearest neighbour method
        :param depth_grid: TOPAZ depth grid [Agrid object]
        :param distance_threshold: Maximum distance (in km) between a GloFAS estuary and the closest TOPAZ ocean point; estuaries further away are dropped in find_estuaries() [int]
        """

        if lazy_mode==False:
            bathy_mask=~depth_grid.depth.mask

            ds_depth=xr.Dataset()
            ds_depth['longitude_bathy']=depth_grid.lon[bathy_mask]
            ds_depth['latitude_bathy'] = depth_grid.lat[bathy_mask]
            ds_depth['depth'] = depth_grid.depth[bathy_mask]
            ds_depth['TOPAZ_index']=np.indices(bathy_mask.shape).flatten()[bathy_mask]

            grid_lon,grid_lat=np.meshgrid(self.ds.longitude,self.ds.latitude)

            self.ds['longitude_grid'] = (['latitude', 'longitude'], grid_lon)
            self.ds['latitude_grid'] = (['latitude', 'longitude'], grid_lat)


            distance,indice=refine_grid(ds_depth.longitude_bathy,ds_depth.latitude_bathy,self.ds['longitude_grid'],self.ds['latitude_grid'],ds_depth.longitude_bathy.shape,self.ds.longitude_grid.shape)

            self.ds['latitude_glofas_on_topaz']=(['latitude', 'longitude'],ds_depth['latitude_bathy'].values[indice[0]])
            self.ds['longitude_glofas_on_topaz'] = (['latitude', 'longitude'], ds_depth['longitude_bathy'].values[indice[0]])
            # self.ds['idx_glofas_on_topaz']= (['latitude', 'longitude'], indice[0])
            self.ds['idx_glofas_on_topaz'] = (['latitude', 'longitude'], ds_depth['TOPAZ_index'].values[indice[0]])
            self.ds['distance_to_topaz']=(['latitude', 'longitude'], distance)

            # Masking cells  that are too far from Topaz grid for given distance threhsold.If distance above threshold, it then becomes Nan (and will be dropped when finding estuaries

            self.ds['distance_to_topaz']=xr.where(self.ds['distance_to_topaz']>distance_threshold,np.nan,self.ds['distance_to_topaz'])

        elif lazy_mode==True: #Lazy mode, we use input grid and copy interpolation
            self.ds['longitude_grid'] = (['latitude', 'longitude'], input_grids.dis_grid.ds['longitude_grid'].data)
            self.ds['latitude_grid'] = (['latitude', 'longitude'], input_grids.dis_grid.ds['latitude_grid'].data)

            self.ds['latitude_glofas_on_topaz'] = (
            ['latitude', 'longitude'],input_grids.dis_grid.ds['latitude_glofas_on_topaz'].data.astype(float))
            self.ds['longitude_glofas_on_topaz'] = (
            ['latitude', 'longitude'], input_grids.dis_grid.ds['longitude_glofas_on_topaz'].data.astype(float))
            self.ds['idx_glofas_on_topaz'] = (['latitude', 'longitude'], input_grids.dis_grid.ds['idx_glofas_on_topaz'].data)
            self.ds['distance_to_topaz'] = (['latitude', 'longitude'], input_grids.dis_grid.ds['distance_to_topaz'].data)

            # Masking cells  that are too far from Topaz grid for given distance threhsold.If distance above threshold, it then becomes Nan (and will be dropped when finding estuaries

            self.ds['distance_to_topaz'] = xr.where(self.ds['distance_to_topaz'] > distance_threshold, np.nan,
                                                    self.ds['distance_to_topaz'])

    def correct_largest_russian_rivers(self,correction_file,correct_GloFAS=True,climatology=False):
        """
        Corrects the GloFAs data for the 3 largest russian rivers: Lena, Yenisei and Ob using a linear transformation specific to each day of the year and each river. This correction is based on ArcticGRO data.
        :param correction_file: NetCDF file containing the correction parameters [str]
        :param correct_GloFAS: Whether to correct the data (True) or not (False) [bool]
        """
        if correct_GloFAS:
            correction = xr.open_dataset(correction_file)

            if climatology==False:
                day_of_year=self.df_dis_Estuary.valid_time[0].dayofyear

                if self.df_dis_Estuary.valid_time[0].year%4==0: #Leap year:
                    if self.df_dis_Estuary.valid_time[0].dayofyear>=60: #After 29 february
                        correction = correction.sel({'day of the year': day_of_year - 2}) #Count 28 Feb twice, so we shift by 1
                    else:
                        correction = correction.sel({'day of the year': day_of_year-1})

                else:
                    correction = correction.sel({'day of the year': day_of_year-1})


            else:
                First_day_of_month=[1,32,60,91,121,	152,182,213,244,274,305,335,365]

                days_of_month=np.arange(First_day_of_month[self.df_dis_Estuary.valid_time[0]-1],
                                        First_day_of_month[self.df_dis_Estuary.valid_time[0]])
                # day_of_year=First_day_of_month[self.df_dis_Estuary.valid_time[0]-1]+15 #to reach approx middle of the month
                correction = correction.sel({'day of the year': days_of_month-1})
                correction = correction.mean(dim='day of the year')


            # The rows of the river mouths themselves (found in find_estuaries), before the river mouths
            # are merged per model cell, so that neighbouring rivers are not corrected
            for river in correction.river.data:
                lat, lon = self.large_river_mouths[river]
                idx = self.df_dis_Estuary.index[(self.df_dis_Estuary['latitude'] == lat) &
                                                (self.df_dis_Estuary['longitude'] == lon)][0]

                slope=correction.sel({'river':river})['Slope'].data #Select slope and intercept for day and river
                intercept=correction.sel({'river':river})['Intercept'].data

                corrected=self.df_dis_Estuary.at[idx,'dis24']*slope+intercept #Correcting the discharge value
                if corrected<0: #The linear correction can give negative discharge (e.g. Lena in late May); set it to zero
                    print('Corrected discharge of %s is negative (%.1f m3/s), set to 0' % (river, corrected))
                    corrected=0.
                self.df_dis_Estuary.at[idx,'dis24']=np.float32(corrected)

            # #To delete the main rivers from total discharge:
            #     print("!!!!!!!!!!!!We are deleting main rivers!!!!!!!!!!!!!!")
            #     self.df_dis_Estuary.at[idx, 'dis24']=0

            print('GloFAS data correction done')

        else:
            print('No correction applied to GloFAS data')


    def propagate_river_plume(self,depth_grid,rradius = 2.5*80000,alongshoreradius = 2.5*200000,river_index_df=None,re_propagate=False,
                              Propagation_cleaning_step=True,sigma=3,pointsource=False):
        """
        Spreads the discharge of each estuary over the nearby ocean cells of the model grid, with
        spreading_weights() of bin/river_spreading.py (shared with the AHYPE v4b river forcing).
        The weights are stored in self.ocean_weights and reused for the following days in lazy mode.

        :param depth_grid: TOPAZ depth grid [A-grid object]
        :param rradius: maximum distance to the coast [m]
        :param alongshoreradius: maximum distance to the estuary [m]
        :param Propagation_cleaning_step: remove weight from cells separated from the estuary by land
        :param sigma: width of the Gaussian smoothing of the weights [grid cells], 0 for none
        :param pointsource: put each estuary's discharge into its own model cell, without spreading
        """
        idm, jdm = depth_grid.idm, depth_grid.jdm #i and j are reversed in fortran (j for lines, i for columns). It's confusing...

        start_time=time.time()

        if re_propagate==False:
            Estuaries_index_list=self.df_dis_Estuary.index
        else:
            Estuaries_index_list=[river_index_df] #We only propagate the flux for the modified rivers

        land = np.asarray(depth_grid.depth.mask).reshape((jdm, idm))
        lon = np.asarray(depth_grid.lon).reshape((jdm, idm))
        lat = np.asarray(depth_grid.lat).reshape((jdm, idm))
        cells = self.df_dis_Estuary['idx_glofas_on_topaz'].loc[Estuaries_index_list].to_numpy().astype(int)
        if pointsource:
            new_weights = spreading_weights(lon, lat, land, cells, strategy='pointsource')
        else:
            new_weights = spreading_weights(lon, lat, land, cells, strategy='diffuse', cleaning=Propagation_cleaning_step,
                                            rradius=rradius, alongshore=alongshoreradius, sigma=sigma)

        # Columns of the estuaries just computed (all of them, or the one being re-propagated)
        positions = [self.df_dis_Estuary.index.get_loc(i) for i in Estuaries_index_list]
        if re_propagate:
            pos = positions[0]
            W = self.ocean_weights
            self.ocean_weights = sp.hstack([W[:, :pos], new_weights, W[:, pos + 1:]], format='csc')
        else:
            self.ocean_weights = new_weights  # reused for the following days in lazy mode

        # Total river flux on the model grid (m3/s per cell), same as in lazy mode
        dis_array = self.df_dis_Estuary['dis24'].to_numpy()
        self.river_flux_TOPAZ = self.ocean_weights @ dis_array

        # Check per-estuary flux conservation
        column_sums = np.asarray(self.ocean_weights.sum(axis=0, dtype=np.float64)).ravel()
        for pos, i in zip(positions, Estuaries_index_list):
            if abs(dis_array[pos] * column_sums[pos] - dis_array[pos]) > 1:
                warnings.warn('Fluxes do not add up for river at position:'+
                      str(self.df_dis_Estuary['longitude_glofas_on_topaz'].loc[i])+','+str(self.df_dis_Estuary['latitude_glofas_on_topaz'].loc[i]))
                print(dis_array[pos] * column_sums[pos],'!=',dis_array[pos])
        print('Plume weights: {n} estuaries, {nnz} nonzero values ({mb:.1f} MB)'.format(
            n=self.ocean_weights.shape[1], nnz=self.ocean_weights.nnz,
            mb=(self.ocean_weights.data.nbytes + self.ocean_weights.indices.nbytes + self.ocean_weights.indptr.nbytes) / 1e6))
        print('Total propagation time:', time.time()-start_time)


    def propagate_river_plume_lazily(self,depth_grid,input_grids):
        """
        Propagates the river discharge from an estuary to the neighbouring points on the TOPAZ grid.
        This process is done by computing a weight for ocean points within a certain radius from the estuary.
        This wieght depends on proximity to the estuary and to the coast.

        The ocean weights aren't computed here but are given directly as an input (lazy version)

        :param depth_grid: Topaz depth grid [Agrid]
        :param input_grids: Grids object containing (among other things) the ocean weights for each estuary [Grids]
        """
        self.ds_topaz_river_flux = xr.Dataset(
            coords={"longitude": depth_grid.lon.data, "latitude": depth_grid.lat.data})

        ocean_weights=input_grids.dis_grid.ocean_weights #Sparse (jdm*idm, n_estuaries) weights from the first day
        estuary_dis=self.df_dis_Estuary['dis24'].to_numpy()

        self.river_flux_TOPAZ=ocean_weights@estuary_dis  # (jdm*idm,) m3/s per cell

    def check_total_fluxes(self):
        """
        Computes the total flux from propagated river flux on the TOPAZ grid and from the Estuaries in the estuaries dataframe.
        And then compares them. These two values should be almost identical, if the different between them is bigger than 1m3/s,
        a warning message will appear.
        """
        Total_propagated_flux=self.river_flux_TOPAZ.sum()

        Total_estuaries_flux=self.df_dis_Estuary['dis24'].sum()

        print('For day {date} : \n'
              'Total propagated flux={Tot_prop_flux:.3f}m3/s \n'
              'Total GloFAs estuaries flux={Tot_estuary_flux:.3f}m3/s'.format(
            date=str(self.ds.valid_time.data)[:10],
            Tot_prop_flux=Total_propagated_flux, Tot_estuary_flux=Total_estuaries_flux))

        if abs(Total_estuaries_flux-Total_propagated_flux)>1:
            warnings.warn("!!!!!! Difference between fluxes sum is unusually high !!!!!!!!!!!!")


    def load_modif_file(self, input_edit_file,do_edit=False):
        """
        The modifications written in the input csv modification file are added in the Dataframe
        :param input_edit_file: input CSV file containing the modifications to be made to estuaries locations [str]
        :param do_edit: Boolean indicating whether or not to apply the edits [bool]
        :return:
        """

        if input_edit_file==None or do_edit==False: #Exit function if no input edit file or do_edit is False
            return None

        df_input_edit=pd.read_csv(input_edit_file)

        for i in range(len(df_input_edit)):

            # For each line in the file, we look for corresponding point in th estuary dataframe
            GloFAS_lon, GloFAS_lat = df_input_edit.loc[i].GloFAS_longitude, df_input_edit.loc[
                i].GloFAS_latitude
            TOPAZ_longitude, TOPAZ_latitude = df_input_edit.loc[i].TOPAZ_longitude, \
                                              df_input_edit.loc[i].TOPAZ_latitude

            try:  # If this point is in the estuary dataframe, we make the replacements (changing the point of estuary on topaz grid)
                self.df_dis_Estuary.loc[(self.df_dis_Estuary.longitude == GloFAS_lon) &
                                                       (self.df_dis_Estuary.latitude == GloFAS_lat),
                                                       ['idx_glofas_on_topaz', 'longitude_glofas_on_topaz',
                                                        'latitude_glofas_on_topaz']] = \
                    [df_input_edit.loc[i].TOPAZ_idx,
                     TOPAZ_longitude, TOPAZ_latitude]

                print("Point at:" + str(GloFAS_lon) + "," + str(GloFAS_lat) + " moved")

            except:  # But the point may not be in the estuary dataframe
                print("Point at:" + str(GloFAS_lon) + "," + str(GloFAS_lat) + " is not in the list of estuaries")



    def isolate_river(self,river_name):
        "Can be used to focus on only one river (and study discharge time series for example)"

        if river_name==None:
            return None

        #Coords of rivers (lon,lat)
        RIVERS={"Glomma":(10.6,56.2)}

        River=RIVERS[river_name]

        distance_to_river = np.sqrt((self.df_dis_Estuary['latitude'] - River[0]) ** 2 + (self.df_dis_Estuary['longitude'] -River[1]) ** 2)


        #We'll find the indices for the 3 rivers, we know what they should be but it may change with GloFAS version
        idx_River=distance_to_river.argmin()

        # self.df_dis_Estuary.at[idx_River,'dis24'] 
