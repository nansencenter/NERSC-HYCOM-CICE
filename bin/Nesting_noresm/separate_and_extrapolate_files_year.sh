#! /bin/bash
module load CDO/1.9.10-iimpi-2022a

#
# Ensure Nesting_Files_PATH is set
if [ -z "$Nesting_Files_PATH" ]; then
    echo "Error: Nesting_Files_PATH is not set."
    echo "Check if you sourced REGION.src"
    exit 1
fi

sleep 5 # so the user can read whats above

# Variables on regular grid
lname=${ESM_Scenario}_gr

year=$1
vari=$2
varie=${vari}_Omon
BIASCORR=false

# if [ $year -le 2014 ]; then
#    d1=${year:0:3}
# elif [ $(($year % 10)) -eq 0 ]; then
#    dyear=$((year-1))
#    d1=${dyear:0:3}
#    d2=${year:0:3}
# else
#    d1=${year:0:3}
#    dyear=$((year+10))
#    d2=${dyear:0:3}
# fi

# Calculate the start of the decade
d1=$(( (year / 10) * 10 ))
# Calculate the end of the decade
d2=$(( d1 + 9 ))

file_found=false
if [ "$vari" != "zos" ]; then
   # Loop through the years of the decade
   for (( i=d1; i<=d2; i++ )); do
    filename=${Nesting_Files_PATH}${varie}_${lname}_${d1}01-${i}12.nc
    if [ -e "$filename" ]; then
        file_found=true
        d2=$i
        break
    fi
   done
   # If no file was found
   if [ "$file_found" = false ]; then
    echo $filename "not found"
    exit 1
   fi


#### Extract to monthly files for 3D varibales
   cdo selyear,${year} ${filename} \
                      ${Nesting_Files_PATH}/${varie}_${lname}_${year}_grid.nc
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      cdo selmon,${mon} ${Nesting_Files_PATH}/${varie}_${lname}_${year}_grid.nc \
                        ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_grid.nc
   done

#### Bias correct temperature and salinity
#### Bias correct nutrient and oxygen, but make sure thye are not negative.
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      if [ "$vari" == "thetao" ]; then BIASCORR=true; fi
      if [ "$vari" == "so" ]; then BIASCORR=true; fi
      if [ "$vari" == "no3" -o "$vari" == "si" ]; then BIASCORR=true; fi
      if [ "$vari" == "po4" -o "$vari" == "o2" ]; then BIASCORR=true; fi
      if [ $BIASCORR = true ]; then
         # cdo sub ${vari}_Omon_${lname}_${year}${smon}_grid.nc \
         #         /cluster/projects/nn9481k/NORESM_bias/bias_${vari}_decal_${smon}.nc \
         #         ${vari}_Omon_${lname}_${year}${smon}_gridbc.nc
         cdo sub ${Nesting_Files_PATH}/${vari}_Omon_${lname}_${year}${smon}_grid.nc \
                 ${Nesting_Files_PATH}/bias_${vari}_decal_${smon}.nc \
                 ${Nesting_Files_PATH}/${vari}_Omon_${lname}_${year}${smon}_gridbc.nc
      fi
   done

   #### For nutrients and oxygen, set negative values to a small value.
   if [ $BIASCORR = true ]; then
      for ((mon=1; mon<=12; mon+=1)); do
         smon=`echo -n 0$mon | tail -2c`
         if [ "$vari" == "no3" -o "$vari" == "po4" -o "$vari" == "o2" -o "$vari" == "si" ]; then
            cdo setrtoc,-inf,0.0001,0.0001 ${Nesting_Files_PATH}/${vari}_Omon_${lname}_${year}${smon}_gridbc.nc \
               ${Nesting_Files_PATH}/${vari}_Omon_${lname}_${year}${smon}_gridbc2.nc
         else
            cp ${Nesting_Files_PATH}/${vari}_Omon_${lname}_${year}${smon}_gridbc.nc \
               ${Nesting_Files_PATH}/${vari}_Omon_${lname}_${year}${smon}_gridbc2.nc
	 fi
     done
   fi

#### Extrapolate missing data
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      if  [ $BIASCORR = true ]; then
         cdo setmisstonn ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_gridbc2.nc \
                         ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_extrap.nc
      else
         cdo setmisstonn ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_grid.nc \
                         ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_extrap.nc
      fi
   done
else
# Variables on native grid
   lname=${ESM_Scenario}_gn
   # Calculate the start of the decade
   d1=$(( (year / 10) * 10 ))
   # Calculate the end of the decade
   d2=$(( d1 + 9 ))

   file_found=false
   # Loop through the years of the decade
   for (( i=d1; i<=d2; i++ )); do
    filename=${Nesting_Files_PATH}${varie}_${lname}_${d1}01-${i}12.nc
    if [ -e "$filename" ]; then
        file_found=true
        d2=$i
        break
    fi
   done
   # If no file was found
   if [ "$file_found" = false ]; then
    echo $filename "not found"
    exit 1
   fi


#### Extract to monthly files
#   cdo selyear,${year} ${Nesting_Files_PATH}/${varie}_${lname}_${d1}01-${year:0:1}10012.nc \
   cdo selyear,${year} ${Nesting_Files_PATH}/${varie}_${lname}_${d1}01-${d2}12.nc \
                       ${Nesting_Files_PATH}/${varie}_${lname}_${year}_grid.nc
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      cdo selmon,${mon} ${Nesting_Files_PATH}/${varie}_${lname}_${year}_grid.nc \
                        ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_grid.nc
   done

#### Extrapolate missing data
   for ((mon=1; mon<=12; mon+=1)); do
      smon=`echo -n 0$mon | tail -2c`
      cdo setmisstonn ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_grid.nc \
                      ${Nesting_Files_PATH}/${varie}_${lname}_${year}${smon}_extrap.nc
   done
fi
