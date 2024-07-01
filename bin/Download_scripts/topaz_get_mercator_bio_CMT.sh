#!/bin/bash
##--------------------------------User-specification
#SBATCH --job-name=TPZ_get_mercator_bio
#SBATCH --time=01:00:00
#SBATCH -n 1
#SBATCH -o  ${HOME}/TOPAZ5/NEMO/GETMERC_bio_t.out 

exec 1> ${HOME}/TOPAZ5/NEMO/get_mercator_bio.log 2>&1

DATADIR=${HOME}/TOPAZ5/NEMO
start_time=`date +%s`

# Define the start and end dates
start_date="2024-01-01"
end_date="2024-01-10"
mamba activate copernicusmarine
username="UUUUUU"
password="XXXXXX"
#--------------------------------
#
#
date
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& now download the BIO file for the ame day &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
cmt_service='GLOBAL_ANALYSISFORECAST_BGC_001_028-TDS'
cmt_dataset1='cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m'
cmt_region1=' -x -180 -X 179.75 -y 30 -Y 90'
cmt_depth=' -z 0.493 -Z 5727.918'
cmt_variables1=' -v no3 -v si -v po4'

# Convert the dates to seconds since the epoch
start_sec=$(date -d "$start_date" +%s)
end_sec=$(date -d "$end_date" +%s)
# Loop through each day in the date range
for (( current_sec=$start_sec; current_sec<=$end_sec; current_sec+=86400 )); do
   #
   current_date=$(date -d "@$current_sec" +%Y-%m-%d)
   echo $current_date
   #dateT1="2024-05-20"
   dateT1="$current_date"
   date1=$(echo $(echo $dateT1|cut -c1-4)$(echo $dateT1|cut -c6-7)$(echo $dateT1|cut -c9-10))
   echo "Downloading GLOBAL Mercator data produced on date: $dateT1 >>>>> $date1 "
   #Assume file does not exist
   #exit_bio=1
   exit_bio=1
   if [ -s ${DATADIR}/global_analysis_forecast_phy_${date1}_bio.nc ] ; then
      echo "Already exist, skip downloading GLOBAL data for $dateT1 12:00:00"
      exit_bio=0
   else
      echo "try to download from MDS server"
      for attempts in 1 2 3 # Try three times (because of frequent connections timeouts)
      do
       #echo "copernicusmarine subset --force-download --overwrite --no-metadata-cache --dataset-id ${cmt_dataset1}\
       #${cmt_variables1} -t "${dateT1}T00:00:00" -T "${dateT1}T00:00:00" ${cmt_depth}\
       #--username $username --password $password  -o ${DATADIR} -f global_analysis_forecast_phy_${date1}_bio.nc"
       copernicusmarine subset --force-download --overwrite --no-metadata-cache --dataset-id ${cmt_dataset1}\
       ${cmt_variables1} -t "${dateT1}T00:00:00" -T "${dateT1}T00:00:00" ${cmt_region1} ${cmt_depth}\
       --username $username --password $password  -o ${DATADIR} -f global_analysis_forecast_phy_${date1}_bio.nc
       #  
       if [ -s ${DATADIR}/global_analysis_forecast_phy_${date1}_bio.nc ] ; then
          echo "SUCCESS downloading GLOBAL bi data for $dateT1 12:00:00"
          exit_bio=0
          break
       else
          echo "FAILURE downloading GLOBAL data for $dateT1 12:00:00"
          exit_bio=1
       	if [[ "${attempts}" < "3" ]] ; then
                echo "Attempt ${attempts}/3 failed. Waiting 30 s until next attempt."
       		sleep 30
       	fi
       fi
      done
      # End try to download from Motu server:----------------
   fi
done
#--------
#--------

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& end download the bio file for the ame day &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#
end_time=`date +%s`
elapsed_time=$((end_time-start_time))
echo "--------->>>>>>>>>>>>>  elapsed_time=$(( $elapsed_time)) sec"
### echo "get_mercator job FINISHED" > ${DATADIR}/get_mercator_9.end_msg
