#! /bin/bash

# Files from the historical run can be found here on NIRD:
# /projects/NS9034K/CMIP6/CMIP/NCC/NorESM2-MM/historical/r1i1p1f1/

module load NCO/5.1.3-iimpi-2022a
module load CDO/1.9.10-iimpi-2022a

dirname_h=Forcing_NorESM        # Historical
dirname_p=Forcing_NorESM_RCP8.5 # Projection
#dirname_p=New_2079 # Projection

decade=$1

startyear=${decade}0
endyear=${decade}0
fleapyear=${decade}$2

# Check leap year:
if [ $(($fleapyear % 4)) -eq 0 ]; then
   echo "Leap year OK"
else
   echo "Leap year not OK"
   exit
fi

# remove the variable "time_bnds" since the atmospheric forcing "klikkker i vinkel" when it tries to read it and set caledar to standard time
rmtbound=true
if [ $rmtbound = true ]; then
   for variable in FLDS FSDS PS QREFHT TREFHT U V; do
      echo ${variable}
      for ((year=${startyear}; year<=${endyear}; year+=1)); do
         if [ $year -le 2014 ]; then
            ncks -C -O -x -v time_bnds ${dirname_h}/${variable}_y${year}.nc tmp.nc
         else
            ncks -C -O -x -v time_bnds ${dirname_p}/${variable}_y${year}.nc tmp.nc
         fi
         cdo setcalendar,standard tmp.nc ${variable}_nt_y${year}.nc
         rm tmp.nc
      done                          
   done
fi

# Add leap day to all varibles except precipitation
addleap=false
if [ $addleap = true ]; then
   for variable in FLDS FSDS PS QREFHT TREFHT U V; do
      for ((year=${fleapyear}; year<=${endyear}; year+=4)); do
         cdo seldate,${year}-02-28T00:00:00 ${variable}_nt_y${year}.nc tmp.nc
         cdo setdate,${year}-02-29T00:00:00 tmp.nc tmp0000.nc
         cdo seldate,${year}-02-28T03:00:00 ${variable}_nt_y${year}.nc tmp.nc
         cdo setdate,${year}-02-29T03:00:00 tmp.nc tmp0300.nc
         cdo seldate,${year}-02-28T06:00:00 ${variable}_nt_y${year}.nc tmp.nc
         cdo setdate,${year}-02-29T06:00:00 tmp.nc tmp0600.nc
         cdo seldate,${year}-02-28T09:00:00 ${variable}_nt_y${year}.nc tmp.nc   
         cdo setdate,${year}-02-29T09:00:00 tmp.nc tmp0900.nc
         cdo seldate,${year}-02-28T12:00:00 ${variable}_nt_y${year}.nc tmp.nc
         cdo setdate,${year}-02-29T12:00:00 tmp.nc tmp1200.nc
         cdo seldate,${year}-02-28T15:00:00 ${variable}_nt_y${year}.nc tmp.nc
         cdo setdate,${year}-02-29T15:00:00 tmp.nc tmp1500.nc
         cdo seldate,${year}-02-28T18:00:00 ${variable}_nt_y${year}.nc tmp.nc
         cdo setdate,${year}-02-29T18:00:00 tmp.nc tmp1800.nc
         cdo seldate,${year}-02-28T21:00:00 ${variable}_nt_y${year}.nc tmp.nc
         cdo setdate,${year}-02-29T21:00:00 tmp.nc tmp2100.nc
         cdo -O mergetime ${variable}_nt_y${year}.nc tmp*00.nc  ${variable}_ntl_y${year}.nc
         rm tmp.nc tmp*00.nc   
         cp ${variable}_ntl_y${year}.nc ${variable}_nt_y${year}.nc
      done 
   done
fi

# Interpolate precipitaion down to 3 hours.
tinterprec=true
fpath=/cluster/work/users/annettes/NORESM_forcing/Daily_NORESM/
fileh=${fpath}/pr_day_nt_NorESM2-MM_historical_r1i1p1f1_gn  # Historical
filep=${fpath}/pr_day_nt_NorESM2-MM_ssp585_r1i1p1f1_gn      # Projection

if [ $tinterprec = true ]; then
   for variable in TPRECIP; do
      for ((year=${startyear}; year<=${endyear}; year+=1)); do
         prev_year=$((year-1))
         next_year=$((year+1))
         # Transition between historical and projction
         if [ $year -lt 2014 ]; then
            cdo mergetime ${fileh}_${prev_year}.nc ${fileh}_${year}.nc ${fileh}_${next_year}.nc tmp.nc
         elif [ $year -eq 2014 ]; then
            cdo mergetime ${fileh}_${prev_year}.nc ${fileh}_${year}.nc ${filep}_${next_year}.nc tmp.nc
         elif  [ $year -eq 2015 ]; then
            cdo mergetime ${fileh}_${prev_year}.nc ${filep}_${year}.nc ${filep}_${next_year}.nc tmp.nc
         elif  [ $year -eq 2100 ]; then
            cdo mergetime ${filep}_${prev_year}.nc ${filep}_${year}.nc tmp.nc
         elif  [ $year -gt 2015 ]; then
            cdo mergetime ${filep}_${prev_year}.nc ${filep}_${year}.nc ${filep}_${next_year}.nc tmp.nc
         fi
         cdo inttime,$year-01-01,00:00:00,3hours tmp.nc tmp1.nc
         cdo seldate,$year-01-01T00:00:00,$next_year-01-01T00:00:00 tmp1.nc tmp2.nc
         cdo mulc,0.001 tmp2.nc tmp3.nc
         cdo chname,pr,TPRECIP tmp3.nc tmp4.nc
         cdo setattribute,TPRECIP@units="m/s" tmp4.nc ${variable}_nt_y${year}.nc
         rm tmp*.nc
      done
   done
fi

# copy the last field from the previous year to the next year to facilitate the change of year in hycom_atmfor.py
startyear=${decade}0
copylastfield=false  # Only necessary for historical it seems.
if [ $copylastfield = true ]; then
   for variable in FLDS FSDS PS QREFHT TREFHT U V TPRECIP; do
#   for variable in TPRECIP; do
      for ((year=${startyear}; year<=${endyear}; year+=1)); do
         prev_year=$((year-1))
         cdo seldate,${year}-01-01T00:00:00 ${variable}_nt_y${prev_year}.nc tmp.nc
         cdo mergetime tmp.nc ${variable}_nt_y${year}.nc tmp1.nc
         mv tmp1.nc ${variable}_nt_y${year}.nc
      done
   done 
fi


