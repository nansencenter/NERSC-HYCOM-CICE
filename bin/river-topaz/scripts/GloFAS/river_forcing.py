import sys
import numpy as np
import xarray as xr
import pandas as pd
from compute_discharge_NoUI import compute_discharge_noUI
import abfile

#-----------------------------------PARAMETERS----------------------------------------------
# data path definition
GloFASdata_path='/cluster/projects/nn9481k/GloFAS_data/TOPAZrunoff_data'
if GloFASdata_path not in sys.path:
   sys.path.append(GloFASdata_path)

input_dataset=GloFASdata_path+'/../data_v40/data_all_year.nc' #path to netcdf file containing raw GloFAS data

output_path='/cluster/work/users/xiejp/work2025/Hriver/' #Where you want the forcing files to be 

start_date='2020-01-01T00:00:00.000000000' # Start date of forcing file with format: YYYY-MM-DDT00:00:00.000000000
end_date='2020-01-02T00:00:00.000000000'   # End date of forcing file


Apply_GloFAS_correction=True #True to Apply GloFAS correction starting from 2001

Edit_estuaries=False #True to take into account estuary edit file and modify the positions for certain estuaries

Propagation_cleaning_step=True #True to prevent river plume to extend beyond land (can create discontinuities)

rradius = 2*80000  #Radius for river plume propagation based on distance to shore [m]
alongshoreradius = 2*160000 #Radius for river plume propagation based on distance to estuary [m]

#---------------------------------------------------------------------------------------------



def date_as_Julian_day(date):
    return (pd.Timestamp(date)-pd.Timestamp(year=1900, month=12, day=31)).days

def create_river_forcing(input_dataset,start_date,end_date,dt=.25):
    """
    Creates river forcing [ab] files with dt period (6 hourly if dt=0.25) from start_date to end_date (included)
    :param input_dataset: Name of the dataset to be used (NetCDF file)
    :param start_date: Start date: yyyy-mm-ddT00:00:00.000000000 [str]
    :param end_date: End date: yyyy-mm-ddT00:00:00.000000000 [str]
    :param dt: Time step for output data AS FRACTION OF DAY--> dt=0.25 => 6 hours period , dt=1 => 1 day period
    :return: NetCDF file with river forcing on TOPAZ grid for every date
    """
    data=xr.open_dataset(input_dataset)

    #check that 'valid_time' field exists
    if 'valid_time' not in data:
        if 'time' not in data:
            print('There is not time information in the dataset')
            return False
        data['valid_time']=data['time'] #'valid_time' field is essential for processing,
        #But it sometimes disappears when merging GloFAS data from different years with CDO command

    data = data.sel(time=slice(start_date, end_date))

    #Create list of times in days since 01/01/1900 with dt as timestep
    start_of_time = pd.Timestamp(year=1900, month=12, day=31) #01/01/1901
    start_date_as_day,end_date_as_day=pd.Timestamp(start_date)-start_of_time,pd.Timestamp(end_date)-start_of_time #days since 01/01/1900
    date_list=np.arange(start_date_as_day.days,end_date_as_day.days+1+dt,dt)

    lazy_mode=False
    input_grids=None
    day_cntr=0


    for value in data.time.values:

        daily_river_discharge=data.sel(time=value)

        correct_GloFAS=(value.astype('datetime64[Y]').astype(int)+1970>=2001) & Apply_GloFAS_correction #correction starts in 2001
        

        grids = compute_discharge_noUI(river_data=daily_river_discharge, lazy_mode=lazy_mode,input_grids=input_grids,correct_GloFAS=correct_GloFAS,Edit_estuaries=Edit_estuaries,Propagation_cleaning_step=Propagation_cleaning_step,
                                       rradius = rradius ,alongshoreradius = alongshoreradius)

        day_cntr+=1

        if lazy_mode==False : #On the first day, we create Dataset

            idm,jdm=grids.topaz_depth_grid.idm,grids.topaz_depth_grid.jdm

            river_forcing = xr.Dataset(
                coords={'longitude': np.arange(0,jdm,1), 'latitude': np.arange(0,idm,1), 'dtime1': date_list})
            lazy_mode = True #After 1st day, go in lazy mode

            river_forcing['river_flux']=(['longitude','latitude','dtime1'],np.zeros((jdm,idm,date_list.shape[0])))

            # Create the output HYCOM file
            print("Creating forcing file")
            if dt<1.0:
               foutfile=output_path+"riverh_{start}_{end}_{dtt}h".format(start=start_date[:10],end=end_date[:10],dtt=dt*24)
            else:
               foutfile=output_path+"riverh_{start}_{end}_{dtt}d".format(start=start_date[:10],end=end_date[:10],dtt=dt)
            ffile = abfile.ABFileRiverForcing(foutfile, "w",idm=idm, jdm=jdm,
                cline1="GloFAS v4.0",
                cline2="river forcing (m/s)")

            input_grids=grids #The day one grids will be used for the other days

        #conversion from m3/s to m/s
        river_flux_topaz = grids.dis_grid.river_flux_TOPAZ.reshape((jdm, idm))/(grids.topaz_grid.scpx*grids.topaz_grid.scpy)
        # river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value))]=river_flux_topaz

 
        N=int(1//dt) #new dim size
        river_flux_topaz_extended=np.stack([river_flux_topaz] * N, axis=2)


        # Convert dtime1 slice to indices
        time_idx = river_forcing.indexes["dtime1"]
        # print(time_idx)
        i0 = time_idx.get_loc(date_as_Julian_day(value))
        # print('i0=',i0)
        i1 = time_idx.get_loc(date_as_Julian_day(value)+1-dt) + 1  # inclusive slice
        # print('i1=',i1)
        

        # print(river_forcing.river_flux.isel(dtime1=slice(i0, i1))[:])

        river_forcing.river_flux.isel(dtime1=slice(i0, i1))[:] =river_flux_topaz_extended


        #Linear interpolation for previous day:
        if date_as_Julian_day(value)!=start_date_as_day.days :
            flux_jm1,flux_j=river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value)-1)],river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value))]
            total_daily=flux_jm1*(N/2+0.5)+(N-1)/2*flux_j
            bias_correction=flux_jm1-total_daily/N
            bias_correction*=0

            bias_correction_extended = np.stack([bias_correction] * N, axis=2)

            weight_0,weight_1=np.arange(1,0,-dt),np.arange(0,1,dt)

            river_interpolated=(river_forcing['river_flux'].loc[dict(dtime1=slice(date_as_Julian_day(value)-1, date_as_Julian_day(value)-dt))]*weight_0).data+ \
                               (river_forcing['river_flux'].loc[dict(dtime1=slice(date_as_Julian_day(value), date_as_Julian_day(value)+1-dt))]*weight_1).data

            river_forcing['river_flux'].loc[
                dict(dtime1=slice(date_as_Julian_day(value) - 1, date_as_Julian_day(value) - dt))]=river_interpolated+bias_correction_extended


            #We also write it in HYCOM file
            period_in_hours=24*dt #For conversion
            # conv = 3600 * period_in_hours
            conv = 1

            for dtime in np.arange(date_as_Julian_day(value)-1,date_as_Julian_day(value),dt):
                ffile.write_field(river_forcing['river_flux'].loc[dict(dtime1=dtime)].data*conv,
                                river_forcing['river_flux'].loc[dict(dtime1=dtime)].data*conv,
                                'rivers', dtime, dt)

        #finally we write last day in the file
        if date_as_Julian_day(value)==end_date_as_day.days:
            period_in_hours = 24 * dt  # For conversion
            ffile.write_field(river_forcing['river_flux'].loc[dict(dtime1=end_date_as_day.days)].data * conv,
                                             river_forcing['river_flux'].loc[dict(dtime1=end_date_as_day.days)].data * conv,
                                            'rivers', end_date_as_day.days, dt)

    ffile.close()

    return river_forcing


river_forcing=create_river_forcing(input_dataset=input_dataset,start_date=start_date,end_date=end_date)
#river_netcdf=river_forcing['river_flux'].transpose("dtime1","latitude","longitude")
#river_netcdf.to_netcdf('{wrkdrt}/riverh_{start}_{end}.nc'.format(wrkdrt=output_path,start=start_date[:10],end=end_date[:10]),unlimited_dims='dtime1')
