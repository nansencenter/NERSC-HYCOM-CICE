#!/usr/bin/env python
import argparse
import os
import sys
import numpy as np
import xarray as xr
import pandas as pd
from compute_discharge_NoUI import compute_discharge_noUI
import abfile
import logging
import datetime

#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False # Dont propagate to parent in hierarchy (determined by "." in __name__)

#---------------------------------------------------------------------------------------------
def date_as_Julian_day(date):
    return (pd.Timestamp(date)-pd.Timestamp(year=1900, month=12, day=31)).days

MAX_GAP = 30  # [days] how far to look before/after a missing value for valid data


def fill_missing(full: xr.Dataset, valid: np.ndarray, iday: int, today: np.ndarray) -> np.ndarray:
    """Return today's discharge with missing values linearly interpolated in time.

    A value is missing where it is NaN today but the cell has data on the first or last day
    of the GloFAS file (`valid`, the river network). GloFAS has a few such gaps (e.g. the
    Mackenzie estuary on 2021-05-07), which would otherwise change the estuary list and crash
    the lazy mode. Each missing value is interpolated linearly between the nearest days with
    data before and after (at most MAX_GAP days away), using only the original file values,
    so the result does not depend on the processed period. At the start or end of the file
    the nearest value is used.

    :param full: GloFAS dataset (whole file) with dis24(time, latitude, longitude)
    :param valid: True where the cell has data on the first or last day of the file
    :param iday: index of today in full.time
    :param today: dis24 of today (latitude, longitude)
    :return: dis24 of today with the gaps filled (today itself if nothing is missing)
    """
    missing = np.flatnonzero(np.isnan(today) & valid)
    if missing.size == 0:
        return today

    def nearest(step: int) -> tuple[np.ndarray, np.ndarray]:
        """Nearest valid values and their distance in days, going back (-1) or forward (+1)."""
        val, dist = np.full(missing.size, np.nan), np.zeros(missing.size)
        for k in range(1, MAX_GAP + 1):
            j = iday + step * k
            todo = np.flatnonzero(np.isnan(val))
            if j < 0 or j >= full.sizes["time"] or todo.size == 0:
                break
            v = full.dis24.isel(time=j).values.reshape(-1)[missing[todo]]
            ok = np.isfinite(v)
            val[todo[ok]], dist[todo[ok]] = v[ok], k
        return val, dist

    v0, k0 = nearest(-1)
    v1, k1 = nearest(+1)
    both = np.isfinite(v0) & np.isfinite(v1)
    fill = np.where(both, v0 + (v1 - v0) * k0 / np.where(both, k0 + k1, 1.), np.where(np.isfinite(v0), v0, v1))
    if np.isnan(fill).any():
        sys.exit("Error: %d GloFAS values on %s are missing for more than %d days"
                 % (np.isnan(fill).sum(), str(full.time.values[iday])[:10], MAX_GAP))
    filled = today.copy()
    filled.reshape(-1)[missing] = fill
    jj, ii = np.unravel_index(missing, today.shape)
    for j, i, f in zip(jj, ii, fill):
        logger.warning("%s: missing GloFAS discharge at %.3fN %.3fE, filled by time interpolation: %.1f m3/s"
                       % (str(full.time.values[iday])[:10], full.latitude.values[j], full.longitude.values[i], f))
    return filled


def create_river_forcing(input_dataset,start_date,end_date,topaz_grid_file,topaz_depth_file,dt=.25):
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

    # Stop if the GloFAS data do not cover the whole period: sel() would silently return
    # fewer days, and the last requested day would then not be written at all.
    t0, t1 = data.time.values[0], data.time.values[-1]
    if np.datetime64(start_date) < t0 or np.datetime64(end_date) > t1:
        sys.exit("Error: requested period %s to %s, but %s only covers %s to %s"
                 % (start_date[:10], end_date[:10], input_dataset, str(t0)[:10], str(t1)[:10]))
    full = data
    valid = np.isfinite(full.dis24.isel(time=0).values) | np.isfinite(full.dis24.isel(time=-1).values)
    data = data.sel(time=slice(start_date, end_date))

    #Create list of times in days since 01/01/1900 with dt as timestep
    start_of_time = pd.Timestamp(year=1900, month=12, day=31) #01/01/1901
    start_date_as_day,end_date_as_day=pd.Timestamp(start_date)-start_of_time,pd.Timestamp(end_date)-start_of_time #days since 01/01/1900
    date_list=np.arange(start_date_as_day.days,end_date_as_day.days+1+dt,dt)

    lazy_mode=False
    input_grids=None
    day_cntr=0

    import time as _time
    t_total_start = _time.time()

    for value in data.time.values:

        t_day_start = _time.time()
        date_str = str(value)[:10]

        daily_river_discharge=data.sel(time=value).load()
        today = daily_river_discharge.dis24.values
        filled = fill_missing(full, valid, full.indexes["time"].get_loc(value), today)
        if filled is not today:
            daily_river_discharge["dis24"] = daily_river_discharge.dis24.copy(data=filled)

        correct_GloFAS=(value.astype('datetime64[Y]').astype(int)+1970>=2001) & Apply_GloFAS_correction #correction starts in 2001


        grids = compute_discharge_noUI(river_data=daily_river_discharge,topaz_grid_file=topaz_grid_file,topaz_depth_file=topaz_depth_file,
                                       lazy_mode=lazy_mode,input_grids=input_grids,correct_GloFAS=correct_GloFAS,Edit_estuaries=Edit_estuaries,Propagation_cleaning_step=Propagation_cleaning_step,
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
            ffile = abfile.ABFileRiverForcing(foutfile, "w", mask=True, idm=idm, jdm=jdm,
                cline1="GloFAS v4.0",
                cline2="river forcing (m/s)")

            land_mask = grids.topaz_depth_grid.depth.mask.reshape((jdm, idm))

            input_grids=grids #The day one grids will be used for the other days

        #conversion from m3/s to m/s
        river_flux_topaz = grids.dis_grid.river_flux_TOPAZ.reshape((jdm, idm))/(grids.topaz_grid.scpx*grids.topaz_grid.scpy)
        # river_forcing['river_flux'].loc[dict(dtime1=date_as_Julian_day(value))]=river_flux_topaz

 
        N=int(1//dt) #new dim size
        river_flux_topaz_extended=np.repeat(river_flux_topaz[:, :, np.newaxis], N, axis=2)


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
                                land_mask,
                                'rivers', dtime, dt)

        #finally we write last day in the file
        if date_as_Julian_day(value)==end_date_as_day.days:
            period_in_hours = 24 * dt  # For conversion
            ffile.write_field(river_forcing['river_flux'].loc[dict(dtime1=end_date_as_day.days)].data * conv,
                                             land_mask,
                                            'rivers', end_date_as_day.days, dt)

        print("Day %s done in %.1f s  (total elapsed: %.1f s)" % (
              date_str, _time.time()-t_day_start, _time.time()-t_total_start))

    ffile.close()

    n_days = day_cntr
    t_total = _time.time() - t_total_start
    print("\nFinished: %d days processed in %.1f s  (%.1f s/day)" % (n_days, t_total, t_total/max(n_days,1)))

    return river_forcing




if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          #tmp = datetime.datetime.strptime( values, "%Y-%m-%d")
          setattr(args, self.dest, tmp)

   parser = argparse.ArgumentParser(description='Prepare HYCOM forcing files from a set of input files')
   parser.add_argument('start_time',    action=DateTimeParseAction, help='Start time in UTC zone. Format = YYYY-mm-dd<THH:MM:SS>')
   parser.add_argument('end_time',      action=DateTimeParseAction, help='Stop  time in UTC zone. Format = YYYY-mm-dd<THH:MM:SS>')
   parser.add_argument('Outfile_drt',   type=str, help='directory to save the output files')
   parser.add_argument('depth_file',    type=str, help='path to the HYCOM depth .a file for this grid configuration')
   parser.add_argument('--dt',          type=float, default=0.25,
                       help='Output record interval as fraction of a day (default 0.25 = 6-hourly, '
                            'required by HYCOM 2.2.98 where rivers must match the atmospheric forcing times; '
                            'use 1.0 = daily for HYCOM 2.3, priver=2, which interpolates between any records)')

   args = parser.parse_args()
   if args.dt <= 0. or args.dt > 1. or abs(1./args.dt - round(1./args.dt)) > 1e-6:
      parser.error('--dt must divide one day evenly (e.g. 0.25, 0.5, 1.0)')

   # Start date of forcing file with format: YYYY-MM-DDT00:00:00.000000000
   start_date=str(args.start_time)[0:10] + 'T00:00:00.000000000' 
   print(start_date)
   # End date of forcing file
   end_date=str(args.end_time)[0:10] + 'T00:00:00.000000000' 
   #end_date=args.end_time
   print(end_date)
   #Where you want the forcing files to be 
   if args.Outfile_drt !="":
      output_path=args.Outfile_drt 
   else:
      output_path='./'

#-----------------------------------PARAMETERS----------------------------------------------

   GloFASdata_path='/cluster/projects/nn9481k/GloFAS_data/TOPAZrunoff_data'
   if GloFASdata_path not in sys.path:
       sys.path.append(GloFASdata_path)
   input_dataset=GloFASdata_path+'/../data_v40/data_all_year.nc' # path to netcdf file
                                                                 #containing raw GloFAS data

   topaz_depth_file=args.depth_file
   topaz_grid_file =os.path.join(os.path.dirname(topaz_depth_file), 'regional.grid.a')

   Apply_GloFAS_correction=True #True to Apply GloFAS correction starting from 2001

   Edit_estuaries=False         #True to take into account estuary edit file and modify
                                # the positions for certain estuaries

   Propagation_cleaning_step=True #True to prevent river plume to extend beyond land 
                                  #(can create discontinuities)

   rradius = 2*90000              #Radius for river plume propagation based on distance
                                  # to shore [m]
   alongshoreradius = 2*180000    #Radius for river plume propagation based on distance 
                                  # to estuary [m]
#---------------------------------------------------------------------------------------------

   river_forcing=create_river_forcing(input_dataset=input_dataset,
                    start_date=start_date,end_date=end_date,
                    topaz_grid_file=topaz_grid_file,topaz_depth_file=topaz_depth_file,
                    dt=args.dt)
   river_netcdf=river_forcing['river_flux'].transpose("dtime1","longitude","latitude")
   river_netcdf.to_netcdf('{wrkdrt}/riverh_{start}_{end}.nc'.format(wrkdrt=output_path,
                    start=start_date[:10],end=end_date[:10]),unlimited_dims='dtime1')
