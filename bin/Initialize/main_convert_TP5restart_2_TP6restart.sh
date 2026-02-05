#!/bin/bash
#Main script to convert TP5 restart file to TP6 restart file
# It calls:
# 1) topaz_restart2archv_TP5_CML.sh     to extract TP5 archv file from TP5 restart file
# 2) topaz_TP5archv_2_TP6archv_CML.sh   to interpolate TP5 archv file to TP6 archv file
# 3) topaz_archv2restart_TP6_CML.sh     to create TP6 restart file from TP6 archv file 
# 4) convert_cice_restart.sh           to interpolate CICE TP5 restart file to TP6 grid
# New restart files are created under the initialize directory of TP6a0.03 experiment
BIND="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/bin/Initialize/"
#O="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize"
SCRTH="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/initialize/SCRATCH"
New_RST="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/initialize/New_Restart"
echo "Convert TP5 restart  file to TP6  restart file"
echo "-------------------------------------------------------"
if [ "$#" -ne 1 ]; then
   echo "$0: Usage: $0 TP5_restart_file "
   echo " "
   echo "Example ../main_convert_TP5restart_2_TP6restart.sh TP5restart.2024_184_00_0000_mem000.a"
   echo " "
   # --- input  TP5  restart  file:
   # --- output TP6  restart  file
   # ../main_convert_TP5restart_2_TP6restart.sh TP5restart.2024_184_00_0000.a
   exit 1
fi
mkdir -p "${SCRTH}"
pushd ${SCRTH}
# All operations are done in the SCRATCH directory
echo "-------------------------------------------------------"
echo "All operations are done in the SCRATCH directory : ${SCRTH} "
echo "-------------------------------------------------------"

TP5_restart_file_in=$1
yyyy=$(echo $TP5_restart_file_in | cut -c12-15)
ddd=$(echo $TP5_restart_file_in | cut -c17-19)
icedate=`date -d "${yyyy}-01-01 + $((10#$ddd-1)) days" "+%F"`
TP5_cice_restart="iced.${icedate}-00000_mem000.nc"
echo "cp ${HOME}/work/sea/TOPAZ5/Forecast14/Forecast14_000/${TP5_restart_file_in}   ${SCRTH}/"
echo "cp ${HOME}/work/sea/TOPAZ5/Forecast14/Forecast14_000/cice/${TP5_cice_restart} ${SCRTH}/"
cp ${HOME}/work/sea/TOPAZ5/Forecast14/Forecast14_000/${TP5_restart_file_in}   ${SCRTH}/
cp ${HOME}/work/sea/TOPAZ5/Forecast14/Forecast14_000/cice/${TP5_cice_restart} ${SCRTH}/
#echo "copying restart ${yyyy}_${strt_day} from ${src} to --------> ${dest} "
#echo "CICE date_f= `date -d "${yyyy}-01-01 + $((10#$strt_day-1)) days" "+%F"` "


echo " %%%%%%%% 1) Extract TP5archv from TP5restarr file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 1) Extract TP5archv from TP5restarr file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 1) Extract TP5archv from TP5restarr file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 1) Extract TP5archv from TP5restarr file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 1) Extract TP5archv from TP5restarr file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
TP5_archv_file_out=$(echo ${TP5_restart_file_in} | sed -e 's/restart/archv/g' -e 's/_0000//g' )
echo "input file ${TP5_restart_file_in}  and output_archv_file= ${TP5_archv_file_out}"
echo ""
echo "exec ${BIND}/topaz_restart2archv_TP5_CML.sh ${TP5_restart_file_in} > topaz_restart2archv_TP5_CML_${TP5_restart_file_in::-2}.log"
${BIND}/topaz_restart2archv_TP5_CML.sh ${TP5_restart_file_in} > topaz_restart2archv_TP5_CML_${TP5_restart_file_in::-2}.log
echo ""
if [ $? -ne 0 ]; then
    echo "Command failed, exiting."
    echo "Check log file: topaz_restart2archv_TP5_CML_${TP5_restart_file_in::-2}.log"
    exit
else
    echo "Successfully created TP5 archv file: ${TP5_archv_file_out}" 
fi
echo ""

echo " %%%%%%%% 2) Interpolate TP5archv file to create TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 2) Interpolate TP5archv file to create TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 2) Interpolate TP5archv file to create TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 2) Interpolate TP5archv file to create TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 2) Interpolate TP5archv file to create TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
TP6_archv_file_out=$(echo ${TP5_archv_file_out} | sed -e 's/TP5/TP6/g')
TP6_restart_file_out=$(echo ${TP5_restart_file_in} | sed -e 's/TP5/TP6/g')
remove_mem000=$(echo ${TP6_restart_file_out} | sed -e 's/_mem000//g' )
TP6_restart_file_out=${remove_mem000}   
echo "convert ${TP5_archv_file_out} to ${TP6_archv_file_out}"
echo ""
echo "exec ${BIND}/topaz_TP5archv_2_TP6archv_CML.sh ${TP5_archv_file_out} > topaz_TP5archv_2_TP6archv_${TP5_archv_file_out::-2}.log"
${BIND}/topaz_TP5archv_2_TP6archv_CML.sh ${TP5_archv_file_out} > topaz_TP5archv_2_TP6archv_${TP5_archv_file_out::-2}.log
echo ""
if [ $? -ne 0 ]; then
    echo "Command failed, exiting."
    echo "Check log file: topaz_TP5archv_2_TP6archv_${TP5_archv_file_out::-2}.log"
    exit
else    
   echo "Successfully created TP6 archv file: ${TP6_archv_file_out}" 
fi
echo ""
echo " %%%%%%%% 3) Create TP6restart file from TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 3) Create TP6restart file from TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 3) Create TP6restart file from TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 3) Create TP6restart file from TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 3) Create TP6restart file from TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo " %%%%%%%% 3) Create TP6restart file from TP6archv file  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------"
echo ""
echo "exec ${BIND}/topaz_archv2restart_TP6_CML.sh ${TP6_archv_file_out} TP6restart.2024_012_00_0000.a > topaz_archv2restart_TP6_CML_${TP6_archv_file_out::-2}.log"
${BIND}/topaz_archv2restart_TP6_CML.sh ${TP6_archv_file_out} TP6restart.2024_012_00_0000.a > topaz_archv2restart_TP6_CML_${TP6_archv_file_out::-2}.log
if [ $? -ne 0 ]; then
    echo "Command failed, exiting."
    echo "Check log file: topaz_archv2restart_TP6_CML_${TP6_archv_file_out::-2}.log"
    exit
else    
   echo "Successfully created restart file ${TP6_restart_file_out} "
fi
echo ""
echo " ----  ---- "
echo " ----  ---- "
echo " %%%%%%%% 4) Now interpolate cice TP5 restart file on TP6 grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ---- "
echo " %%%%%%%% 4) Now interpolate cice TP5 restart file on TP6 grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ---- "
echo " %%%%%%%% 4) Now interpolate cice TP5 restart file on TP6 grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ---- "
echo " %%%%%%%% 4) Now interpolate cice TP5 restart file on TP6 grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ---- "
echo " %%%%%%%% 4) Now interpolate cice TP5 restart file on TP6 grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ---- "
echo " %%%%%%%% 4) Now interpolate cice TP5 restart file on TP6 grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ---- "
echo ""
echo "exec ${BIND}/convert_cice_restart.sh ${TP5_cice_restart} > convert_cice_restart_${TP5_cice_restart::-2}.log"
${BIND}/convert_cice_restart.sh ${TP5_cice_restart} > convert_cice_restart_${TP5_cice_restart::-2}.log   
echo " "
echo " ----  ---- "
if [ $? -ne 0 ]; then
    echo "Command failed, exiting."
    echo "Check log file: convert_cice_restart_${TP5_cice_restart::-2}.log"
    exit
else    
   echo "Successfully created CICE restart file TP6_iced.${icedate}-00000.nc"
fi
echo ""
echo " ----  ---- "
echo " checking ${SCRTH}/${TP6_restart_file_out} and moving new restart files to ${New_RST}/"
if [ -e ${SCRTH}/${TP6_restart_file_out} ]; then
   echo "mv ${SCRTH}/${TP6_restart_file_out} ${New_RST}/"
   mv ${SCRTH}/${TP6_restart_file_out} ${New_RST}/
   mv ${SCRTH}/${TP6_restart_file_out%.a}.b ${New_RST}/
   # move CICE restart file
   echo "mv ${SCRTH}/TP6_iced.${icedate}-00000.nc ${New_RST}/"
   mv ${SCRTH}/TP6_iced.${icedate}-00000.nc ${New_RST}/
   echo "Moved new restart files to ${New_RST}/"
   echo "Successfully created restart file ${SCRTH}/${TP6_restart_file_out::-2}.a"
else
   echo "Error: something went wrong"
fi
echo ""
echo "-------------------------------------------------------"
echo "Exiting  ${SCRTH}"
popd
echo "-------------------------------------------------------"
