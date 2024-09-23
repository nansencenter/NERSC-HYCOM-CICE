#!/bin/bash
##--------------------------------User-Specification
#SBATCH --job-name=TPZ_get_mercator
#SBATCH --time=01:00:00
#SBATCH -n 1
#SBATCH -o  ${HOME}/TOPAZ5/NEMO/GETMERC_t.out 

exec 1> ${HOME}/TOPAZ5/NEMO/get_mercator_phys_9.log 2>&1

DATADIR=${HOME}/TOPAZ5/NEMO
start_time=`date +%s`

# Define the start and end dates
start_date="2024-01-01"
end_date="2024-01-10"
mamba activate copernicusmarine
module load CDO/2.3.0-eccodes-aec-cmor-hpc1-intel-2023a-eb
username="UUUUUU"
password="XXXXXX"
#--------------------------------
#
#
date
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& now download the PHYS file for the ame day &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
cmt_service='GLOBAL_ANALYSISFORECAST_PHY_001_024'
#cmt_dataset='cmems_mod_glo_phy_anfc_0.083deg_P1D-m'
cmt_region=' -x -180 -X 179.91668701171875 -y 30 -Y 90'
cmt_depth=' -z 0.493 -Z 5727.918'
echo "get mercator job STARTED" > ${DATADIR}/get_mercator_9.start_msg

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
   exit1_st=1
   exit2_st=1
   if [ -s ${DATADIR}/global_analysis_forecast_phy_${date1}.nc ] ; then
      echo "Already exist, skip downloading GLOBAL data for $dateT1 12:00:00"
      exit1_st=0
      exit2_st=0
   else
      echo "Try to download from CMEMS catalogue"
      count_f=0
      # Variables are found in different datasets
      for foo in {1..4} 
      do
         if [ ${foo} -eq 1 ] ; then
            cmt_variables=' -v zos'
            cmt_dataset='cmems_mod_glo_phy_anfc_0.083deg_P1D-m'
         elif [ ${foo} -eq 2 ] ; then
            cmt_variables=' -v thetao'
            cmt_dataset='cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m'
         elif [ ${foo} -eq 3 ] ; then
            cmt_variables=' -v so'
            cmt_dataset='cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m'
         else
            cmt_variables=' -v uo -v vo'
            cmt_dataset='cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m'
         fi
         #
         echo " "
         for attempts in 1 2 3  # Try five times (because of frequent connections timeouts)
         do
           date
           #echo "copernicusmarine subset --force-download --overwrite --no-metadata-cache --dataset-id ${cmt_dataset}\
           #${cmt_variables}  -t "${dateT1}T00:00:00" -T "${dateT1}T00:00:00" ${cmt_region} ${cmt_depth}\
           #    --username $username --password $password  -o ${DATADIR} -f global_analysis_forecast_phy${foo}_${date1}.nc"
           copernicusmarine subset --force-download --overwrite --no-metadata-cache --dataset-id ${cmt_dataset}\
              ${cmt_variables}  -t "${dateT1}T00:00:00" -T "${dateT1}T00:00:00" ${cmt_region} ${cmt_depth}\
               --username $username --password $password  -o ${DATADIR} -f global_analysis_forecast_phy${foo}_${date1}.nc  
           #---
           if [ -s ${DATADIR}/global_analysis_forecast_phy${foo}_${date1}.nc ] ; then
               echo "SUCCESS downloading GLOBAL data for global_analysis_forecast_phy${foo}_${date1}.nc"
               count_f=$(( $count_f + 1 ))
               echo "count=$count_f"
               break
            else
              echo "FAILURE downloading GLOBAL data for global_analysis_forecast_phy${foo}_${date1}.nc"
              if [[ "${attempts}" < "3" ]] ; then
                 echo "Attempt ${attempts}/5 failed. Waiting 60 s until next attempt."
                 sleep 10
              fi
            fi
         done
      done
      if [ ${count_f} -eq 4 ] ; then
          exit1_st=0
          exit2_st=0
      else
          echo "Error---count=${count_f}"
      fi
   fi
   #---
   if [ ${exit1_st} -eq 0 ] && [ ${exit2_st} -eq 0 ] ; then
      echo "get mercator finished"
      # if exist merge files with CDO
      if [ -s ${DATADIR}/global_analysis_forecast_phy1_${date1}.nc ] ; then
         pushd ${DATADIR} 
         echo "cdo -z zip_1 merge global_analysis_forecast_phy{1..4}_${date1}.nc global_analysis_forecast_phy_${date1}.nc"
         cdo -z zip_1 merge global_analysis_forecast_phy{1..4}_${date1}.nc global_analysis_forecast_phy_${date1}.nc
         # remove partial files
         rm -f global_analysis_forecast_phy[1234]_${date1}.nc
         popd
      fi
      echo "get_mercator job FINISHED" > ${DATADIR}/get_mercator_9.end_msg
   else
      echo "error! get_mercator job unsuccessfully FINISHED file"
   fi
   #
done
#--------
#--------

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& end download the bio file for the ame day &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#
end_time=`date +%s`
elapsed_time=$((end_time-start_time))
echo "--------->>>>>>>>>>>>>  elapsed_time=$(( $elapsed_time)) sec"
### echo "get_mercator job FINISHED" > ${DATADIR}/get_mercator_9.end_msg
