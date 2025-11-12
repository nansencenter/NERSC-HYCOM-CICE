import copy
import time
import xarray as xr
from utils.grid_operations import refine_grid
from utils.grid_operations import transform_coordinates
from utils.grid_operations import apply_ocean_weight
from utils.grid_operations import cleaning_discharge
import numpy as np
import pandas as pd
# import imageio.v3 as iio
import skimage as ski
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.use('TkAgg')
import warnings
from scipy.ndimage import gaussian_filter #Import gausian_filter function

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

    def find_estuaries(self,uparea_threshold=1e9,elevation_threshold=100,lazy_mode=False,input_grids=None):
        """
        Finds the cells on the grid that correspond to estuaries by using condition on LDD, upstream area and elevation
        :param uparea_threshold: minimum upstream area (in m^2) for a cell to be a potential estuary [int]
        :param elevation_threshold: maximum elevation (in m) for a cell to be a potential estuary [int]
        """

        if lazy_mode==False:
            self.__uparea_threshold=uparea_threshold
            self.__elevation_threshold = elevation_threshold


            mask_estuaries=np.where(self.ds['ldd']== 5 , 1,np.nan) #Mask based on LDD
            mask_estuaries=mask_estuaries*np.where(self.ds['uparea']> uparea_threshold , 1,np.nan) #Mask based on upstream area: lowX=180000km2 (approx. size of the river Rhine)
            mask_estuaries=mask_estuaries*np.where(self.ds['elevation']< elevation_threshold , 1,np.nan) #Mask based on elevation

            self.mask_estuaries=mask_estuaries #We put the mask in GloGAS grid object so that it can be re used

        elif lazy_mode==True: #In lazy mode, we don't compute mask
            mask_estuaries=input_grids.dis_grid.mask_estuaries

        # And now we apply mask to Discharge data:
        # ds_dis*=mask_estuaries
        self.ds['dis_Estuary'] = self.ds[self._GloFAS_grid__field_name] * mask_estuaries



        # We transform xarray dataset into a panda Dataframe because otherwise
        # It's impossible to drop NaN along multiple dimensions

        self.df_dis_Estuary = self.ds.to_dataframe().dropna(subset='dis_Estuary') # Only keep the points that are estuaries

        self.df_dis_Estuary = self.df_dis_Estuary.reset_index()  # Just a trick to make life easier when locating an index

        self.df_dis_Estuary.drop(columns=['longitude_grid','latitude_grid']) #That's useless

        self.df_dis_Estuary = self.df_dis_Estuary.round({'longitude': 5, 'latitude': 5}) #rounding up the edges

        #We will aslo find the duplicated rows: if two estuary have the same idx on topaz grid (same point) then they are combined and fluxes are sumed

        #Careful here with the field names !! (some may be missing, change the aggregate functions if needed)
        aggregation_functions = {'dis24':'sum', 'valid_time':'first', 'surface':'first', 'land':'first', 'ldd':'first', 'uparea':'sum', 'elevation':'first',
       'longitude_grid':'first', 'latitude_grid':'first', 'latitude_glofas_on_topaz':'first',
       'longitude_glofas_on_topaz':'first', 'distance_to_topaz':'first',
       'dis_Estuary':'sum','latitude':'first','longitude':'first'}
       #  aggregation_functions = {'dis24': 'sum', 'valid_time': 'first', 'land': 'first',
       #                           'ldd': 'first', 'uparea': 'sum', 'elevation': 'first',
       #                           'longitude_grid': 'first', 'latitude_grid': 'first',
       #                           'latitude_glofas_on_topaz': 'first',
       #                           'longitude_glofas_on_topaz': 'first', 'distance_to_topaz': 'first',
       #                           'dis_Estuary': 'sum', 'latitude': 'first', 'longitude': 'first'}

        # In lazy mode, we dont have LDD, uparea, etc so we change agg function slightly
        lazy_aggregation_functions={'dis24':'sum', 'valid_time':'first', 'surface':'first', 'land':'first',
       'longitude_grid':'first', 'latitude_grid':'first', 'latitude_glofas_on_topaz':'first',
       'longitude_glofas_on_topaz':'first', 'distance_to_topaz':'first',
       'dis_Estuary':'sum','latitude':'first','longitude':'first'}
       #  lazy_aggregation_functions = {'dis24': 'sum', 'valid_time': 'first', 'land': 'first',
       #                                'longitude_grid': 'first', 'latitude_grid': 'first',
       #                                'latitude_glofas_on_topaz': 'first',
       #                                'longitude_glofas_on_topaz': 'first', 'distance_to_topaz': 'first',
       #                                'dis_Estuary': 'sum', 'latitude': 'first', 'longitude': 'first'}

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

    def get_glofas_into_topaz(self,depth_grid,distance_threshold=100,lazy_mode=False,input_grids=None):
        """
        Interpolates the GloFAS point (originally on the GloFAS grid) onto the TOPAZ depth grid using a nearest neighbour method
        :param depth_grid: TOPAZ depth grid [Agrid object]
        :param distance_threshold: Maximum distance (in m) between a GloFAS estuary and the closest TOPAZ grid depth point to take the estuary into account [int]
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


            Lena=(72,129) #LAT, LON, idx 354
            Yenisei=(70.90,83.30) #idx 268
            Ob=(66.5,69.8) #idx 234

            distance_to_lena = np.sqrt((self.df_dis_Estuary['latitude'] - Lena[0]) ** 2 + (self.df_dis_Estuary['longitude'] - Lena[1]) ** 2)

            distance_to_yenisei = np.sqrt((self.df_dis_Estuary['latitude'] - Yenisei[0]) ** 2 + (self.df_dis_Estuary['longitude'] - Yenisei[1]) ** 2)

            distance_to_ob = np.sqrt((self.df_dis_Estuary['latitude'] - Ob[0]) ** 2 + (self.df_dis_Estuary['longitude'] - Ob[1]) ** 2)

            #We'll find the indices for the 3 rivers, we know what they should be but it may change with GloFAS version
            idx_Lena=distance_to_lena.argmin()
            idx_Yenisei=distance_to_yenisei.argmin()
            idx_Ob=distance_to_ob.argmin()

            IDX={'Lena':idx_Lena,'Yenisei':idx_Yenisei,'Ob':idx_Ob}

            for river in correction.river.data:
                idx=IDX[river]

                slope=correction.sel({'river':river})['Slope'].data #Select slope and intercept for day and river
                intercept=correction.sel({'river':river})['Intercept'].data

                self.df_dis_Estuary.at[idx,'dis24']=(self.df_dis_Estuary.at[idx,'dis24']*slope+intercept).astype(np.float32)#Correcting the discharge value

            # #To delete the main rivers from total discharge:
            #     print("!!!!!!!!!!!!We are deleting main rivers!!!!!!!!!!!!!!")
            #     self.df_dis_Estuary.at[idx, 'dis24']=0

            print('GloFAS data correction done')

        else:
            print('No correction applied to GloFAS data')


    def propagate_river_plume(self,depth_grid,rradius = 2.5*80000,alongshoreradius = 2.5*200000,river_index_df=None,re_propagate=False,
                              Propagation_cleaning_step=True):
        """
        Propagates the river discharge from an estuary to the neighbouring points on the TOPAZ grid.
        This process is done by computing a weight for ocean points within a certain radius from the estuary.
        This weight depends on proximity to the estuary and to the coast.

        A cleaning process is done to ensure that the flux isn't propagated to points that are separated from the estuary by land.

        :param depth_grid: TOPAZ depth grid [A-grid object]
        :param rradius: river radius to define distance max to land [km]
        :param alongshoreradius: radius to define max distance to river origin (estuary) [km]
        """
        idm, jdm = depth_grid.idm, depth_grid.jdm #i and j are reversed in fortran (j for lines, i for columns). It's confusing...

        start_time=time.time()

        if re_propagate==False:

            # We create a dataset. Each dataArray in the dataset containns fluxes associated with a specific estuary on every point of the TOPAZ grid (this way we can easily modify fluxes for one estuary)
            self.ds_topaz_river_flux=xr.Dataset(coords={"longitude":depth_grid.lon.data,"latitude":depth_grid.lat.data})
            Estuaries_index_list=self.df_dis_Estuary.index

            self.ds_topaz_ocean_weights=xr.Dataset(coords={"longitude":np.arange(0,jdm,1),"latitude":np.arange(0,idm,1), "Estuary_index":self.df_dis_Estuary.index})
            self.ds_topaz_ocean_weights['ocean_weights']=(['longitude','latitude','Estuary_index'],np.zeros((jdm,idm,self.df_dis_Estuary.index.shape[0])))

        else:
            Estuaries_index_list=[river_index_df] #We only propagate the flux for the modified rivers

        #extracting depth from TOPAZ
        depth=depth_grid.grid.read_field('depth').data
        depth[depth>1e4]=0

        # distance to closest land cell
        TOPAZ_geodetic=np.vstack((depth_grid.lon.data,depth_grid.lat.data)).T
        TOPAZ_cartesian=transform_coordinates(TOPAZ_geodetic) #go from lon/lat to ECEF to comptue distances
        X,Y,Z=TOPAZ_cartesian[:,0],TOPAZ_cartesian[:,1],TOPAZ_cartesian[:,2]

        dist_to_land = np.zeros_like(depth_grid.depth.data) * np.nan

        # project ocean points on land to find closest landpoint with k-d tree
        distances, indices = refine_grid(depth_grid.lon[depth_grid.depth.mask], depth_grid.lat[depth_grid.depth.mask],
                                         depth_grid.lon, depth_grid.lat,
                                         depth_grid.lon[depth_grid.depth.mask].shape,depth_grid.lon.shape)

        # --Connected components method to find isolated lands
        depth_as_image = ((depth_grid.depth.mask).astype(int) * 200).reshape((jdm, idm))
        image = np.stack([depth_as_image] * 3, axis=2)
        image = image.astype(np.float64)
        gray_image = ski.color.rgb2gray(image)
        # denoise the image with a Gaussian filter
        blurred_image = ski.filters.gaussian(gray_image, sigma=0.1)
        # mask the image according to threshold
        binary_mask = blurred_image > 0.5
        labeled_image, count = ski.measure.label(binary_mask, connectivity=2, return_num=True)

        #For each ocean point,find label of closest land
        # closest_land=labeled_image.flatten()*np.nan

        labeled_image_onland=labeled_image.flatten()[depth_grid.depth.mask]
        closest_land=labeled_image_onland[indices[0]] #Label of closest isolated land area for each TOPAZ point
        closest_land = closest_land.reshape((jdm, idm))
        # colored_label_image = ski.color.label2rgb(labeled_image, bg_label=0)

        dist_to_land[~depth_grid.depth.mask] = distances[
            ~depth_grid.depth.mask]*1000  # For ocean points we apply distance to closest land point and get in into [m]

        outside_rradius=np.where(dist_to_land>rradius)

        for i in Estuaries_index_list:

            river_flux = depth_grid.depth.data * 0
            mini,minj=np.unravel_index(self.df_dis_Estuary['idx_glofas_on_topaz'].loc[i],(jdm,idm))

            pos_estuary=transform_coordinates((self.df_dis_Estuary['longitude_glofas_on_topaz'].loc[i],
                                               self.df_dis_Estuary['latitude_glofas_on_topaz'].loc[i]))


            x_estuary,y_estuary,z_estuary=pos_estuary[:,0],pos_estuary[:,1],pos_estuary[:,2]

            dist_to_estuary=np.sqrt((x_estuary-X)**2+(y_estuary-Y)**2+(z_estuary-Z)**2)*1000

            dist_to_estuary=np.ma.masked_array(dist_to_estuary,mask=depth_grid.depth.mask)


            #assign weight to ocean points based on distance to estuary
            ocean_weight=np.zeros_like(dist_to_estuary.data)


            outside_alongshore=np.where(dist_to_estuary>alongshoreradius)
            ocean_weight[~dist_to_estuary.mask]=np.exp(-dist_to_estuary[~dist_to_estuary.mask]/alongshoreradius)
            ocean_weight[outside_alongshore]*=0

            #Assign extra weight based on distance to shore
            ocean_weight[~dist_to_estuary.mask]*=np.exp(-2*dist_to_land[~dist_to_estuary.mask]/rradius)
            ocean_weight[outside_rradius] *= 0
            #
            # Correct shore distance weight
            land_area = closest_land[mini, minj] #Find label of land area closest to the estuary
            other_lands=np.where(closest_land.flatten()!=land_area) #Find points that have different closest land area
            ocean_weight[other_lands]*=0 #If points have different closest land area than estuary, then they have no weight



            #Remove isolated points
            #For every ocean point with weight, we try to see if there are land points between this point and the estuary
            weight_mask=np.where(ocean_weight!=0)
            #print("Start cleaning process for river ",i)

            #We do cleaning process using Numba Just In Time method to accelerate it
            if Propagation_cleaning_step:
                ocean_weight_copy=copy.deepcopy(ocean_weight)
                ocean_weight=cleaning_discharge(ocean_weight_copy,weight_mask,mini,minj,idm,jdm,
                                                depth_grid.depth.mask.reshape((jdm, idm)))
            else:
                print('skip river plume cleaning')

            #Add smoothing step to avoid discontinuities:
            ocean_weight_grid=np.reshape(ocean_weight,(jdm,idm))
            smoothed_ocean_weight = gaussian_filter(ocean_weight_grid, sigma=3)

            smoothed_ocean_weight[depth_grid.depth.mask.reshape((jdm, idm))]=0 #delete weights on land points

            ocean_weight=smoothed_ocean_weight.flatten() #And back to flat array
            
            #Normalize the weights:
            sum_weight=np.nansum(ocean_weight)
            ocean_weight=ocean_weight/sum_weight

            estuary_flux=self.df_dis_Estuary['dis24'].loc[i]*ocean_weight
            river_flux+=estuary_flux #We add it to the total river flux

            # River flux grid for this estuary
            self.ds_topaz_river_flux['River_flux_estuary_'+str(i)]=(["longitude"],river_flux)

            #We also save the ocean weight grid for this particular estuary:
            self.ds_topaz_ocean_weights['ocean_weights'].loc[dict(Estuary_index=i)]=ocean_weight.reshape((jdm,idm))

            # print("5:", time.time())

            #We check that the sum of propagated flux is equal to estuary discahrge
            if np.sqrt((np.nansum(estuary_flux)-self.df_dis_Estuary['dis24'].loc[i])**2)>1:
                warnings.warn('Fluxes do not add up for river at position:'+
                      str(self.df_dis_Estuary['longitude_glofas_on_topaz'].loc[i])+','+str(self.df_dis_Estuary['latitude_glofas_on_topaz'].loc[i]))
                print(np.nansum(estuary_flux),'!=',self.df_dis_Estuary['dis24'].loc[i])

            else:
#                print("river {river_number},  {Topaz_lon:.3f}, {Topaz_lat:.3f} : flux {river_flux:.3f} \n Sum of propagated flux={sum_prop_flux:.3f} with max local={max_prop_flux:.3f}".format(river_number=str(i), river_flux=self.df_dis_Estuary['dis24'].loc[i], sum_prop_flux=np.nansum(estuary_flux),max_prop_flux=np.nanmax(estuary_flux),Topaz_lon=self.df_dis_Estuary['longitude_glofas_on_topaz'].loc[i], Topaz_lat=self.df_dis_Estuary['latitude_glofas_on_topaz'].loc[i]))
                print("river {river_number},  {Topaz_lon:.3f}, {Topaz_lat:.3f} : flux {river_flux:.3f} ".format(river_number=str(i), river_flux=self.df_dis_Estuary['dis24'].loc[i], Topaz_lon=self.df_dis_Estuary['longitude_glofas_on_topaz'].loc[i], Topaz_lat=self.df_dis_Estuary['latitude_glofas_on_topaz'].loc[i]))


        #We know sum the fluxes of each estuary to obtain final grid with flux on every cell
        self.river_flux_TOPAZ=depth_grid.depth.data*0
        for data_vars in self.ds_topaz_river_flux.data_vars:
            self.river_flux_TOPAZ+=self.ds_topaz_river_flux[data_vars].data

        print('Total propagation time:', time.time()-start_time)
        #35s with no Numba for ocean_weight computation and Numba for cleaning process


    def propagate_river_plume_lazily(self,depth_grid,input_grids):
        """
        Propagates the river discharge from an estuary to the neighbouring points on the TOPAZ grid.
        This process is done by computing a weight for ocean points within a certain radius from the estuary.
        This wieght depends on proximity to the estuary and to the coast.

        The ocean weights aren't computed here but are given directly as an input (lazy version)

        :param depth_grid: Topaz depth grid [Agrid]
        :param input_grids: Grids object containing (among other things) the ocean weights for each estuary [Grids]
        """
        # We create a dataset. Each dataArray in the dataset containns fluxes associated with a specific estuary on every point of the TOPAZ grid (this way we can easily modify fluxes for one estuary)
        self.ds_topaz_river_flux = xr.Dataset(
            coords={"longitude": depth_grid.lon.data, "latitude": depth_grid.lat.data})

        ds_topaz_ocean_weights=input_grids.dis_grid.ds_topaz_ocean_weights #Getting the ocean weights from Grids

        ocean_weight_matrix=ds_topaz_ocean_weights['ocean_weights'].data
        estuary_dis_matrix=self.df_dis_Estuary['dis24'].to_numpy()

        river_flux=ocean_weight_matrix@estuary_dis_matrix #Faster that way

        self.river_flux_TOPAZ=river_flux[:,:]

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
