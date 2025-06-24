#! /bin/bash
module load CDO/1.9.10-iimpi-2022a

source_nesting_experiment_path=$PWD
source ../REGION.src
#
echo ''
echo 'Using '${ESM_ID}' scenario'
echo 'The script will look for decadal files under '${Nesting_Files_PATH}'. Modify REGION.src if necessary'
echo 'The script assumes you are running in an experiment folder located in the Nesting Region (e.g. ESMa1.00/expt_01.0)'
echo ''
#
# Ensure Nesting_Files_PATH is set
if [ -z "$Nesting_Files_PATH" ]; then
    echo "Error: Nesting_Files_PATH is not set."
    exit 1
else
    echo $Nesting_Files_PATH' found'
fi

sleep 5 # so the user can read whats above
syear=$1
eyear=$2
cstr=$3 # e.g. _Omon_NorESM2-MM_historical_r1i1p1f1_gr_

cd $Nesting_Files_PATH

# Extract individual years from 10-year files
for ((year=$syear; year<=$eyear; year+=1)); do
   dec=`echo -n $year | head -3c`
   echo $year $dec 
   for vari in so thetao uo vo no3 o2 po4 si; do 
      cdo selyear,${year} ${vari}${cstr}${dec}001-${dec}[0123456789]12.nc ${vari}${cstr}${year}.nc
      for ((mon=1; mon<=12;mon+=1)); do
         smon=`echo -n 0$mon | tail -2c`
	 cdo selmon,${mon} ${vari}${cstr}${year}.nc ${vari}${cstr}${year}_${smon}.nc
      done
   done
done

# merge to monhtly files
for ((mon=1; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   for vari in so thetao uo vo no3 o2 po4 si; do 
      echo merging for month $smon
      cdo mergetime ${vari}${cstr}*_${smon}.nc ${vari}${cstr}${smon}_all.nc
      cdo timmean ${vari}${cstr}${smon}_all.nc ${vari}${cstr}${smon}_clim.nc
   done
done

for vari in so thetao uo vo no3 o2 po4 si; do
   cdo mergetime ${vari}${cstr}*_clim.nc ${vari}${cstr}clim.nc
done

