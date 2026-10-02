import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('TkAgg')

def correct_GloFAS_monthly_average(input_dataset):
    """
    This function corrects discharge from rivers Lena, Ob and Yenisei for a monthly discharge .nc GloFAS file
    corresponding to specific year.
    Both input and output are on GloFAS grid, we just correct the cells corresponding to the mentionned rivers
    :param input_dataset: Name of file contaiing monthly average on GloFAS grid.
    """
    GloFAS_grid=xr.open_dataset(input_dataset)


    #First we have to find the cells corresponding to estuaries for our rivers

    LDD_grid=xr.open_dataset('../data/river/auxiliary_files__glofas_v4_0/ldd_glofas_v4_0.nc')
    UPAREA_grid=xr.open_dataset('../data/river/auxiliary_files__glofas_v4_0/uparea_glofas_v4_0.nc')

    LON,LAT=np.meshgrid(GloFAS_grid.longitude.data,GloFAS_grid.latitude.data)
    m, n = LON.shape

    LDD=LDD_grid.ldd.data
    UPAREA=UPAREA_grid.uparea.data

    estuaries = np.where((LDD == 5) & (UPAREA>2e9))


    Lena = (72, 129)  # LAT, LON, idx 354
    Yenisei = (70.90, 83.30)  # idx 268
    Ob = (66.5, 69.8)  # idx 234

    dist_to_lena=np.sqrt((LAT-Lena[0])**2+(LON-Lena[1])**2)
    arg_lena=estuaries[0][np.argmin(dist_to_lena[estuaries])],estuaries[1][np.argmin(dist_to_lena[estuaries])]

    dist_to_yenisei=np.sqrt((LAT-Yenisei[0])**2+(LON-Yenisei[1])**2)
    arg_yenisei=estuaries[0][np.argmin(dist_to_yenisei[estuaries])],estuaries[1][np.argmin(dist_to_yenisei[estuaries])]

    dist_to_ob=np.sqrt((LAT-Ob[0])**2+(LON-Ob[1])**2)
    arg_ob=estuaries[0][np.argmin(dist_to_ob[estuaries])],estuaries[1][np.argmin(dist_to_ob[estuaries])]

    Args_rivers={'Lena':arg_lena,'Ob':arg_ob,'Yenisei':arg_yenisei}

    #Now we will do the correction
    correction = xr.open_dataset('../data/river/correction_files/largest_river_correction_2001_2019.nc')

    for month in GloFAS_grid.valid_time.data:
        First_day_of_month = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365]

        days_of_month = np.arange(First_day_of_month[month - 1],
                                  First_day_of_month[month])

        correction_month = correction.sel({'day of the year': days_of_month})
        correction_month = correction_month.mean(dim='day of the year')

        for river in correction_month.river.data:

            arg=Args_rivers[river]

            slope = correction_month.sel({'river': river})['Slope'].data  # Select slope and intercept for day and river
            intercept = correction_month.sel({'river': river})['Intercept'].data

            GloFAS_grid.dis24.values[month-1,arg[0],arg[1]]=\
                GloFAS_grid.dis24.values[month-1,arg[0],arg[1]]*slope + intercept

    return GloFAS_grid




for year in range(2001,2016):
    print(year)
    file='../data/river/GloFAS_monthly/monthly_clim/monthly_clim_{year}.nc'.format(year=year)
    modified_Grid=correct_GloFAS_monthly_average(file)
    modified_Grid.to_netcdf('../data/river/GloFAS_monthly/monthly_clim_corrected/monthly_mean_{year}_corrected.nc'.format(year=year))