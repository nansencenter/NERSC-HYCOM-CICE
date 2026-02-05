#!/bin/bash
#Script to convert CICE TP5 restart file to CICE TP6 restart file
# 
RID="TP6a0.03"
O="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize"
SCRTH="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize/SCRATCH"
New_RST="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize/New_Restart"
echo "Convert TP5 cice restart file to TP6  restat file  file"
#example ./interp_isub_TP5_2_TP6.sh TP5archm.2024_184_12.a
echo "-------------------------------------------------------"
if [ "$#" -ne 1 ]; then
   echo "$0: Usage: $0 TP5 input_TP5_cice_restart_file "
   echo " "
   echo "Example ./convert_cice_restart.sh iced.2024-10-22-00000_mem000.nc"
   echo " "
   # --- input  TP5  cice  file:
   # --- output TP6  cice  file
   exit 1
fi

TP5_cice_input_file=$1
echo "input file $TP5_cice_input_file"


source  ${HOME}/sea/TOPAZ5/Realtime_exp/bin/mod_load_py3



echo ""
echo "------ converting TP5 file to TP6 file------"
mkdir -p "${SCRTH}"
#
pushd ${SCRTH}
echo "Entering  --------${SCRTH}--------"
echo ""
#clean up
echo "PWD=$PWD"
echo ""
echo "Transferring TP5 cice restart file to SCRATCH dir"
echo "cp ${HOME}/work/sea/TOPAZ5/Forecast14/Forecast14_000/cice/${TP5_cice_input_file} ./"
echo "ln -sf ${HOME}/sea/TOPAZ5/NERSC-HYCOM-CICE/TP5a0.06/topo/cice_grid.nc cice_grid_coarse.nc"
echo "ln -sf ${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/topo/cice_grid_merged_kmt.nc cice_grid_fine.nc"
cp ${HOME}/work/sea/TOPAZ5/Forecast14/Forecast14_000/cice/${TP5_cice_input_file} .
ln -sf ${HOME}/sea/TOPAZ5/NERSC-HYCOM-CICE/TP5a0.06/topo/cice_grid.nc cice_grid_coarse.nc
ln -sf ${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/topo/cice_grid_merged_kmt.nc cice_grid_fine.nc
#
echo "python ${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/bin/cice_convert_restart.py ${TP5_cice_input_file} cice_grid_coarse.nc cice_grid_fine.nc"
python ${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/bin/cice_convert_restart.py ${TP5_cice_input_file} cice_grid_coarse.nc cice_grid_fine.nc

#

out_cicef=TP6_$(echo "$TP5_cice_input_file" | sed -e 's/_mem000//g')
echo mv "interpolated_${TP5_cice_input_file}" "$out_cicef"
mv "interpolated_${TP5_cice_input_file}" "$out_cicef"
#out_cicef="interpolated_${TP5_cice_input_file}"
echo "output file $out_cicef"
echo " ----  ---- "
echo " ----  ---- "

if [ -e ${SCRTH}/${out_cicef} ]; then
   #mv ${SCRTH}/${TP6_output_file::-2}.a ${O}/
   #mv ${SCRTH}/${TP6_output_file::-2}.b ${O}/
   echo "Successfully created restart file ${SCRTH}/${out_cicef}"
   #echo " ${SCRTH}/${out_cicef} ${New_RST}/"
   #mv  ${SCRTH}/${out_cicef} ${New_RST}/
else
   echo "Error: something went wrong"
   exit 1
fi
echo "Exiting  ${SCRTH}"
popd
echo "-------------------------------------------------------"
