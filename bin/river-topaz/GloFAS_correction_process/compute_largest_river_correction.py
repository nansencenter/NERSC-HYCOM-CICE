import xarray as xr
import pandas as pd
from functions import *
import os
from datetime import date, timedelta

#----------------INFOS--------------------------

#ArcticGRO info

NAME=['Lena', 'Ob','Yenisei']

lat_gro=[70.68, 66.63,67.43]

lon_gro=[127.39, 66.60,86.48]

path_GRO_Dataframes='/cluster/projects/nn9481k/GloFAS_data/RiverCorrection/GRO_Dataframes/' #Path to ArcticGRO daily dataframes (dataframes computed using file compute_ArcticGRO_Dataframes.py)
path_GloFAS_NC_files='/cluster/work/users/adracc/GloFAS4.0/data/' #Path to daily GloFAS data as NetCDF files

start_year=2001
end_year=2020

def compute_largest_river_correction(path_GRO_Dataframes,path_GloFAS_NC_files,start_year=2010,end_year=2020):
    """
    Computes the corrections for 3 largest Russian rivers: Yenisei, Lena and Ob. This correction is done using ArtciGRO data.
    For each river, a linear regression is computed for each day of the year in the form: ArcticGRO=a*GloFAS+b. The a and b params for each river and each day are saved in a NetDCF file and will be used to Correct daily GloFAS discharge for the rivers in question
    :param path_GRO_Dataframes: path to pandas dataframes containg daily ArcticGRO data. These dataframes are computed using a different script [str]
    :param path_GloFAS_NC_files: path to folder containing yearly data (NetCDF files) for GloFAS [str]
    :param start_year: start year of period used to compute correction (one year=one point in regression) [int]
    :param end_year: end year of period used to compute correction [int]
    """
    YEAR=[i for i in range(start_year,end_year)] #list of years

    file=path_GloFAS_NC_files+'data_'+str(start_year)+'.nc' #select file for given year, it will only be used to extract GloFAS lon/lat grid. Any year will do.

    try:
        ds_find_cell=xr.open_dataset(file) #just used to find closest cell to a river
    except:
        raise FileNotFoundError
    
    LON,LAT=np.meshgrid(ds_find_cell.longitude.data,ds_find_cell.latitude.data)

    df_slope=pd.DataFrame(columns=NAME,index=np.arange(365))
    df_intercept=pd.DataFrame(columns=NAME,index=np.arange(365))

    for i in range(len(NAME)):
        print(NAME[i])

        river=NAME[i]
        lon_river_gro,lat_river_gro=lon_gro[i],lat_gro[i]

        ii_min=find_closest_cells_to_river(ds_find_cell,lon_river_gro,lat_river_gro) #We find closest cell for this river

        #Then we'll find the corresponding ArtcicGRO Dataframe
        path = path_GRO_Dataframes
        river_GRO_file = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and \
                    river in i and 'daily' in i][0]
        
        try:
            df_GRO=pd.read_csv(path+river_GRO_file)
        except:
            raise FileNotFoundError

        # day_cntr=0 #Restet the day counter for each river

        GloFAS_dis=np.zeros((365,len(YEAR)))
        ArcticGRO_dis=np.zeros((365,len(YEAR))) #will contain discharge for day of year j and for every year

        for j in range(len(YEAR)):
            print(YEAR[j])
            year=YEAR[j]


            file=path_GloFAS_NC_files+'data_'+str(year)+'.nc' #select file for given year

            try:
                ds=xr.open_dataset(file)
            except: raise FileNotFoundError


            ds_closest=ds.sel(longitude=LON[ii_min[:,0],ii_min[:,1]],latitude=LAT[ii_min[:,0],ii_min[:,1]],) #returns only cell clostest to river
        

            daily_dis_GloFAS = (ds_closest.dis24.to_numpy())  #daily discharge in m3/s
            daily_dis_GloFAS = np.max(daily_dis_GloFAS, axis=(1,2))

            try:
                daily_dis_GRO=df_GRO.loc[df_GRO['year']==year]['discharge'].values
            except:
                daily_dis_GRO=np.nan*np.ones_like(daily_dis_GloFAS)
                print('Cannot find '+str(year)+' for river '+river+' in ArcticGRO data')

            if daily_dis_GRO.shape!=daily_dis_GloFAS.shape:
                print("shape problem with river "+river+' in year '+str(year))
                daily_dis_GRO=np.nan*np.ones_like(daily_dis_GloFAS)

            if daily_dis_GloFAS.shape[0]==366: #in case of leap year, we delete february 29th
                daily_dis_GloFAS=np.delete(daily_dis_GloFAS,60,0)
                daily_dis_GRO=np.delete(daily_dis_GRO,60,0)


            GloFAS_dis[:,j]=daily_dis_GloFAS #adding daily discharge for year j into arrays
            ArcticGRO_dis[:,j]=daily_dis_GRO

        print('starting lin reg')
        for k in range(365): #Compute lin reg for every day of the year
            A=GloFAS_dis[k,:].reshape((len(YEAR),1))
            A=np.hstack((A,np.ones((len(YEAR),1))))

            B=ArcticGRO_dis[k,:].reshape((len(YEAR),1))

            m, c = np.linalg.lstsq(A, B)[0]

            df_slope.at[k,river]=m
            df_intercept.at[k,river]=c

    df_slope.to_csv('Correction_slope')
    df_intercept.to_csv('Correction_intercept')

    ds_corr=xr.Dataset(coords={'day of year':np.arange(1,366,1),'river':NAME})
    ds_corr['Slope']=(['day of the year','river'], df_slope.values)
    ds_corr['Intercept']=(['day of the year','river'], df_intercept.values)

    # print(ds_corr['Slope'].sel({'river':'Lena'}))

    ds_corr.to_netcdf('largest_river_correction_{start}_{end}.nc'.format(start=start_year,end=end_year))

compute_largest_river_correction(path_GRO_Dataframes=path_GRO_Dataframes,path_GloFAS_NC_files=path_GloFAS_NC_files,start_year=start_year,end_year=end_year)
