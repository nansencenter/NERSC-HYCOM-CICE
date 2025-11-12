import numpy as np
import xarray as xr
import pandas as pd
from compute_discharge_NoUI import compute_discharge_noUI
import abfile


def date_as_Julian_day(date):
    return (pd.Timestamp(date)-pd.Timestamp(year=1900, month=1, day=1)).days

def create_GloFAS_clim_river_forcing(input_dataset):
    """
    Creates river forcing [ab] files as monthly climatology. In the output files, there are 12 grids, on for each month.
    No inter-annual variability
    :param input_dataset: Name of the dataset to be used (NetCDF file)
    """
    data=xr.open_dataset(input_dataset)


    #Create list of times in days since 01/01/1900 with dt as timestep
    date_list=np.arange(1,13,1)


    lazy_mode=False
    input_grids=None
    cntr=0
    for month in date_list:
        grids = compute_discharge_noUI(river_data=data, lazy_mode=lazy_mode,input_grids=input_grids,
                                       climatology=True,month=month,correct_GloFAS=False)

        if lazy_mode==False : #On the first day, we create Dataset


            idm,jdm=grids.topaz_depth_grid.idm,grids.topaz_depth_grid.jdm

            river_forcing = xr.Dataset(
                coords={'longitude': np.arange(0,jdm,1), 'latitude': np.arange(0,idm,1), 'month': date_list})
            lazy_mode = True #After 1st month, go in lazy mode

            river_forcing['river_flux']=(['longitude','latitude','month'],np.zeros((jdm,idm,date_list.shape[0])))

            # Create the output HYCOM file
            ffile = abfile.ABFileRiver(
                "../../data/river_forcing/GloFAS_climatology/GloFAS_clim_forcing_monthly_no_correction_LargeRadius_LowUparea", "w",
                idm=idm, jdm=jdm,
                cline1="GloFAS climatology",
                cline2="river forcing [m/s]")

            input_grids=grids #The day one grids will be used for the other days


        river_flux_topaz = grids.dis_grid.river_flux_TOPAZ.reshape((jdm, idm))
        # river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value))]=river_flux_topaz

        #conversion from m3/s to m/s
        river_forcing['river_flux'].loc[dict(month=month)] = river_flux_topaz/(grids.topaz_grid.scpx*grids.topaz_grid.scpy)
        # river_forcing['river_flux'].loc[dict(month=month)] = river_flux_topaz


        #We also write it in HYCOM file
        Nb_day_month=[31,28,31,30,31,30,31,31,30,31,30,31]
        conv=1

        ffile.write_field(river_forcing['river_flux'].loc[dict(month=month)].data*conv,
                        river_forcing['river_flux'].loc[dict(month=month)].data*conv,
                        'rivers', month)

        cntr+=1

    ffile.close()

    return river_forcing

def create_annual_GloFAS_clim_river_forcing(input_folder,start_year=1979,end_year=2015):
    """
    Creates river forcing [ab] files as monthly climatology with inter-annual variability.
    In the output files, we have 12 grids for each year (one for each month), so in total Nb_of_years*12 grids..
    This river forcing file can be used with HYCOM 2.3 but not HYCOM 2.2
    :param input_folder: Name of the folder containing the datasets to be used (NetCDF file)
    """

    #Create list of times in days since 01/01/1900 with dt as timestep
    MONTHS=np.arange(1,13,1)


    lazy_mode=False
    input_grids=None
    cntr=0
    YEARS=np.arange(start_year,end_year+1,1)

    for year in YEARS:
        for month in MONTHS:
            data=xr.open_dataset('{folder_path}monthly_clim_{year}.nc'.format(folder_path=input_folder,year=year))
            grids = compute_discharge_noUI(river_data=data, lazy_mode=lazy_mode,input_grids=input_grids,
                                           climatology=True,month=month,correct_GloFAS=False)

            if lazy_mode==False : #On the first day, we create Dataset


                idm,jdm=grids.topaz_depth_grid.idm,grids.topaz_depth_grid.jdm

                river_forcing = xr.Dataset(
                    coords={'longitude': np.arange(0,jdm,1), 'latitude': np.arange(0,idm,1), 'month': MONTHS})
                lazy_mode = True #After 1st month, go in lazy mode

                river_forcing['river_flux']=(['longitude','latitude','month'],np.zeros((jdm,idm,MONTHS.shape[0])))

                # Create the output HYCOM file
                ffile = abfile.ABFileRiver(
                    "../../data/river_forcing/GloFAS_climatology/GloFAS_annual_clim_forcing_{start_year}_{end_year}".format(start_year=start_year,end_year=end_year),
                    "w", idm=idm, jdm=jdm,
                    cline1="GloFAS climatology {start_year}-{end_year}".format(start_year=start_year,end_year=end_year),
                    cline2="river forcing [m/s]")

                input_grids=grids #The day one grids will be used for the other days


            river_flux_topaz = grids.dis_grid.river_flux_TOPAZ.reshape((jdm, idm))
            # river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value))]=river_flux_topaz

            #conversion from m3/s to m/s
            river_forcing['river_flux'].loc[dict(month=month)] = river_flux_topaz/(grids.topaz_grid.scpx*grids.topaz_grid.scpy)
            # river_forcing['river_flux'].loc[dict(month=month)] = river_flux_topaz


            #We also write it in HYCOM file
            Nb_day_month=[31,28,31,30,31,30,31,31,30,31,30,31]
            conv=1

            ffile.write_field(river_forcing['river_flux'].loc[dict(month=month)].data*conv,
                            river_forcing['river_flux'].loc[dict(month=month)].data*conv,
                            'rivers', month)

            cntr+=1

    ffile.close()

    return river_forcing


def create_GloFAS_yearly_discharge(start_year=1979,end_year=2015):
    """
    Creates yearly flux from GloFAS data using monthly average.

    For each year we get the montly average file,
    this file contains average discharge for the 12 month of a specif year.
    The montly discharge is processed to get total estuary discharge from that month on TOPAZ grid.
    And we sum monthly discharge to get average yearly discharge
    :param input_dataset: Name of the dataset to be used (NetCDF file)
    :param start_date: Start date: yyyy-mm-ddT00:00:00.000000000 [str]
    :param end_date: End date: yyyy-mm-ddT00:00:00.000000000 [str]
    :return: NetCDF file with river forcing on TOPAZ grid for every date
    """

    YEAR=np.arange(start_year,end_year+1,1)

    river_forcing = xr.Dataset(
        coords={'year': YEAR})

    river_forcing['yearly_flux']=('year',np.zeros(len(YEAR)))

    Nb_day_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    for yr in range(start_year,end_year+1,1):
        print(yr)
        if yr>=2002: correctec_GloFAS=True
        else: correctec_GloFAS=False
        input_dataset='../../data/river/GloFAS_monthly/monthly_clim/monthly_clim_{year}.nc'.format(year=yr)
        data=xr.open_dataset(input_dataset)


        #Create list of times in days since 01/01/1900 with dt as timestep
        date_list=np.arange(1,13,1)


        lazy_mode=False
        input_grids=None
        cntr=0

        yearly_river_flux=0
        for month in date_list:
            grids = compute_discharge_noUI(river_data=data, lazy_mode=lazy_mode,input_grids=input_grids,
                                           climatology=True,month=month,correct_GloFAS=correctec_GloFAS)

            if lazy_mode==False : #On the first day, we create Dataset

                lazy_mode = True #After 1st month, go in lazy mode

                input_grids=grids #The day one grids will be used for the other days


            yearly_river_flux += grids.dis_grid.river_flux_TOPAZ.sum()*Nb_day_month[cntr]*3600*24*1e-9
            # river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value))]=river_flux_topaz
            cntr+=1
        river_forcing['yearly_flux'].loc[dict(year=yr)]=yearly_river_flux

    return river_forcing

def create_GloFAS_yearly_discharge_single_river(start_year=1979,end_year=2015):
    """
    Creates yearly flux from GloFAS data using monthly average.

    For each year we get the montly average file,
    this file contains average discharge for the 12 month of a specif year.
    The montly discharge is processed to get total estuary discharge from that month on TOPAZ grid.
    And we sum monthly discharge to get average yearly discharge
    :param input_dataset: Name of the dataset to be used (NetCDF file)
    :param start_date: Start date: yyyy-mm-ddT00:00:00.000000000 [str]
    :param end_date: End date: yyyy-mm-ddT00:00:00.000000000 [str]
    :return: NetCDF file with river forcing on TOPAZ grid for every date
    """

    YEAR = np.arange(start_year, end_year + 1, 1)

    river_forcing = xr.Dataset(
        coords={'year': YEAR})

    river_forcing['yearly_flux'] = ('year', np.zeros(len(YEAR)))

    Nb_day_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    for yr in range(start_year, end_year + 1, 1):
        print(yr)
        if yr >= 2002:
            correctec_GloFAS = True
        else:
            correctec_GloFAS = False
        input_dataset = '../../data/river/GloFAS_monthly/monthly_clim/monthly_clim_{year}.nc'.format(year=yr)
        data = xr.open_dataset(input_dataset)

        # Create list of times in days since 01/01/1900 with dt as timestep
        date_list = np.arange(1, 13, 1)

        lazy_mode = False
        input_grids = None
        cntr = 0

        yearly_river_flux = 0
        for month in date_list:
            grids = compute_discharge_noUI(river_data=data, lazy_mode=lazy_mode, input_grids=input_grids,
                                           climatology=True, month=month, correct_GloFAS=correctec_GloFAS)

            if lazy_mode == False:  # On the first day, we create Dataset

                lazy_mode = True  # After 1st month, go in lazy mode

                input_grids = grids  # The day one grids will be used for the other days

            yearly_river_flux += grids.dis_grid.river_flux_TOPAZ.sum() * Nb_day_month[cntr] * 3600 * 24 * 1e-9
            # river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value))]=river_flux_topaz
            cntr += 1
        river_forcing['yearly_flux'].loc[dict(year=yr)] = yearly_river_flux

    return river_forcing


input_dataset='../../data/river/GloFAS_monthly/GloFAS_global_clim_1979_2015.nc'
input_folder='../../data/river/GloFAS_monthly/monthly_clim/'
river_forcing=create_annual_GloFAS_clim_river_forcing(input_folder=input_folder,start_year=1979,end_year=1981)
# river_forcing=create_GloFAS_yearly_discharge()

# river_forcing.to_netcdf('../../data/river_forcing/GloFAS_climatology/GloFAS_clim_forcing.nc')