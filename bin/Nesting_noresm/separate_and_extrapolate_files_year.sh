#! /bin/bash
module load CDO/1.9.10-iimpi-2022a

# Variables on regular grid
lname=NorESM2-MM_historical_r1i1p1f1_gr
#lname=NorESM2-MM_ssp585_r1i1p1f1_gr

year=$1
vari=$2
varie=${vari}_Omon
BIASCORR=false

if [ $year -le 2014 ]; then
   d1=${year:0:3}
elif [ $(($year % 10)) -eq 0 ]; then
   dyear=$((year-1))
   d1=${dyear:0:3}
else
   d1=${year:0:3}
fi

if [ "$vari" != "zos" ]; then
#### Extract to monthly files for 3D varibales
   cdo selyear,${year} ${varie}_${lname}_${d1}001-${d1}912.nc \
                      ${varie}_${lname}_${year}_grid.nc
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      cdo selmon,${mon} ${varie}_${lname}_${year}_grid.nc \
                        ${varie}_${lname}_${year}${smon}_grid.nc
   done

#### Bias correct temperature and salinity
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      if [ "$vari" == "thetao" ]; then BIASCORR=true; fi
      if [ "$vari" == "so" ]; then BIASCORR=true; fi
      if [ $BIASCORR = true ]; then
         cdo sub ${vari}_Omon_${lname}_${year}${smon}_grid.nc \
                 Climatology/bias_${vari}_decal_${smon}.nc \
                 ${vari}_Omon_${lname}_${year}${smon}_gridbc.nc
      fi
   done

#### Extrapolate missing data
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      if  [ $BIASCORR = true ]; then
         cdo setmisstonn ${varie}_${lname}_${year}${smon}_gridbc.nc \
                         ${varie}_${lname}_${year}${smon}_extrap.nc
      else
         cdo setmisstonn ${varie}_${lname}_${year}${smon}_grid.nc \
                         ${varie}_${lname}_${year}${smon}_extrap.nc
      fi
   done
else
# Variables on native grid
#   lname=NorESM2-MM_historical_r1i1p1f1_gn
   lname=NorESM2-MM_ssp585_r1i1p1f1_gn

#### Extract to monthly files
   cdo selyear,${year} ${varie}_${lname}_${d1}101-${year:0:1}10012.nc \
                       ${varie}_${lname}_${year}_grid.nc
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      cdo selmon,${mon} ${varie}_${lname}_${year}_grid.nc \
                        ${varie}_${lname}_${year}${smon}_grid.nc
   done

#### Extrapolate missing data
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      cdo setmisstonn ${varie}_${lname}_${year}${smon}_grid.nc \
                      ${varie}_${lname}_${year}${smon}_extrap.nc
   done
fi
