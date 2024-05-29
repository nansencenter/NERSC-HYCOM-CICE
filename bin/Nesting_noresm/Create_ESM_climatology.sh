#! /bin/bash
module load CDO/1.9.10-iimpi-2022a

syear=1981
eyear=2010

cstr="_Omon_NorCPM1_dcppA-assim_r20i2p1f1_gr_" #196001-196912.nc

# Extract individual years from 10-year files
for ((year=$syear; year<=$eyear; year+=1)); do
   dec=`echo -n $year | head -3c`
   echo $year $dec 
   for vari in so thetao; do 
      cdo selyear,${year} ${vari}${cstr}${dec}001-${dec}912.nc ${vari}${cstr}${year}.nc
      for ((mon=1; mon<=12;mon+=1)); do
         smon=`echo -n 0$mon | tail -2c`
	 cdo selmon,${mon} ${vari}${cstr}${year}.nc ${vari}${cstr}${year}_${smon}.nc
      done
   done
done

# merge to monhtly files
for ((mon=1; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   for vari in so thetao; do 
      echo merging for month $smon
      cdo mergetime ${vari}${cstr}*_${smon}.nc ${vari}${cstr}${smon}_all.nc
      cdo timmean ${vari}${cstr}${smon}_all.nc ${vari}${cstr}${smon}_clim.nc
   done
done

for vari in so thetao; do
   cdo mergetime ${vari}${cstr}*_clim.nc ${vari}${cstr}clim.nc
done

