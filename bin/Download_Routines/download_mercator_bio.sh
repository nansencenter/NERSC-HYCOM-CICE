#!/bin/bash
# usage:
# ./download_mercator_bio.sh jday1 jday2 year
# e.g. 
# ./download_mercator_bio.sh 1 10 1994 

start_time=`date +%s`
p_indx=1

#-------------------Start modify these according your case----------------------------------
# load some modules (if necessary) - may need to be updated depending on machin
# module load Python/3.8.2-GCCcore-9.3.0

# add path to local motuclient
MOTUCLIENT="/cluster/home/annettes/.local/lib/python3.8/site-packages/motuclient/motuclient.py" # Instruction: put link to motuclient

#
DATADIR=""#Instruction: write path to folder to download the data to.
bindir="~/HYCOM-CICE/NERSC-HYCOM-CICE/hycom/MSCPROGS/bin/" #Instruction: write path to execuables in MSCPROGS
username="" #Instruction: write your Copernicus marine user name
password="" #Instruction: write your Copernicus marine password
julday1=1  #set a julian day  
#-------------------End moodify these according your case----------------------------------

for ((julday1=$1; julday1<=$2; julday1+=1)); do
source ~/.bashrc
date1=$(echo $(${bindir}/jultodate $julday1 $3 1 1))
dateT1=$(echo $(echo $date1|cut -c1-4)-$(echo $date1|cut -c5-6)-$(echo $date1|cut -c7-8))
echo "Target day is on julian day-----------------------------------> $julday1 = $dateT1"
echo "Downloading GLOBAL Mercator data produced on julian day $julday1 = $dateT1"
echo " "
echo "Output>>: ${DATADIR}/global_analysis_forecast_bio_${date1}.nc" 
echo " "
motu_site='http://my.cmems-du.eu/motu-web/Motu'
motu_service='GLOBAL_MULTIYEAR_BGC_001_029-TDS'
motu_dataset='cmems_mod_glo_bgc_my_0.25_P1D-m'
motu_region=' -x -180 -X 179.91667175293 -y 30 -Y 90' # for TOPAZ5
motu_depth=' -z 0.493 -Z 5727.918'                    # for TOPAZ5
#motu_region=' -x -10 -X 70 -y -60 -Y 10'               # for Agulhas
#motu_depth=' -z 0 -Z 6500'                             # for Agulhas
#motu_region=' -x -10 -X 70 -y 55 -Y 80'                # for Lofoten
#motu_depth=' -z 0 -Z 6000'                             # for Lofoten
motu_variables1=' -v no3 -v si'
motu_variables2=' -v po4 -v o2'

#Assume file does not exist
exit1_st=1
exit2_st=1

echo "get mercator job STARTED" > ${DATADIR}/get_mercator_1.start_msg

if [ -s ${DATADIR}/global_analysis_forecast_bio_${date1}.nc ] ; then
   echo "Already exist, skip downloading GLOBAL data for $dateT1 12:00:00"
   exit1_st=0
   exit2_st=0
else
   # Try to download from Motu server if ftp fails:----------------
   #if [ ! -f ${DATADIR}/global_analysis_forecast_bio_${date1}.nc ] ; then
   if [ ${exit1_st} -eq 1 ] && [ ${exit2_st} -eq 1 ] ; then
      echo "try to download from motu-server server"
      for attempts in 1 2 3 4 5 # Try five times (because of frequent connections timeouts)
      do
         echo $MOTUCLIENT -u $username -p $password -m ${motu_site} -s ${motu_service} -d ${motu_dataset} ${motu_region} -t "${dateT1} 12:00:00" -T "${dateT1} 12:00:00" ${motu_depth} ${motu_variables1} -o ${DATADIR} -f global_analysis_forecast_bio1_${date1}.nc
         $MOTUCLIENT -u $username -p $password -m ${motu_site} -s ${motu_service} -d ${motu_dataset} ${motu_region} -t "${dateT1} 12:00:00" -T "${dateT1} 12:00:00" ${motu_depth} ${motu_variables1} -o ${DATADIR} -f global_analysis_forecast_bio1_${date1}.nc
         if [ -s ${DATADIR}/global_analysis_forecast_bio1_${date1}.nc ] ; then
            echo "SUCCESS downloading GLOBAL data for $dateT1 12:00:00"
            exit1_st=0
            sleep 5
            break
         else
            echo "FAILURE downloading GLOBAL data for $dateT1 12:00:00"
            exit1_st=1
      		if [[ "${attempts}" < "5" ]] ; then
                  echo "Attempt ${attempts}/5 failed. Waiting 60 s until next attempt."
      			sleep 60
      		fi
         fi
      done
      # Downloading second file
      for attempts in 1 2 3 4 5
      do
         $MOTUCLIENT -u $username -p $password -m ${motu_site} -s ${motu_service} -d ${motu_dataset} ${motu_region} -t "${dateT1} 12:00:00" -T "${dateT1} 12:00:00" ${motu_depth} ${motu_variables2} -o ${DATADIR} -f global_analysis_forecast_bio2_${date1}.nc
         if [ -s ${DATADIR}/global_analysis_forecast_bio2_${date1}.nc ] ; then
            echo "SUCCESS downloading GLOBAL data for $dateT1 12:00:00"
            exit2_st=0
            break
         else
            echo "FAILURE downloading GLOBAL data for $dateT1 12:00:00"
            exit2_st=1
      		if [[ "${attempts}" < "5" ]] ; then
                  echo "Attempt ${attempts}/5 failed. Waiting 60 s until next attempt."
      			sleep 60
      		fi
         fi
      done
      # End try to download from Motu server:----------------
   fi
fi
#---
#---

# load CDO to use CDO
module load CDO/1.9.10-intel-2021b
# check that each portion downloaded
if [ ${exit1_st} -eq 0 ] && [ ${exit2_st} -eq 0 ] ; then
   echo "get mercator finished"
   # if exist merge files with CDO
   if [ -s ${DATADIR}/global_analysis_forecast_bio1_${date1}.nc ] ; then
      cdo merge ${DATADIR}/global_analysis_forecast_bio1_${date1}.nc \
                ${DATADIR}/global_analysis_forecast_bio2_${date1}.nc \
                ${DATADIR}/global_analysis_forecast_bio_${date1}.nc
      # remove partial files
      rm -f ${DATADIR}/global_analysis_forecast_bio[12]_${date1}.nc
   fi
   echo "get_mercator job FINISHED" 
  echo "get_mercator job FINISHED" > ${DATADIR}/get_mercator_1.end_msg
else
   echo "error! get_mercator job unsuccessfully FINISHED" > ${DATADIR}/get_mercator_1.err_msg
fi


end_time=`date +%s`
elapsed_time=$((end_time-start_time))
echo "--------->>>>>>>>>>>>> motu-client elapsed_time=$(( $elapsed_time)) sec "
done
### [ $exit_st == 0 ] && echo "get mercator job finished" \
### || echo "error! get_mercator job unsuccessfully FINISHED" > ${DATADIR}/get_mercator_1.err_msg
### 
### echo "get_mercator job FINISHED" > ${DATADIR}/get_mercator_1.end_msg
