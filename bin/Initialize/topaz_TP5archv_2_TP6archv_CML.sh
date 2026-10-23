#!/bin/bash
#
RID="TP6a0.03"
O="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize"
SCR="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize/SCRATCH"
TP6_subregion="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/ConvertTP5"
echo "Convert TP5  archv  file to TP6  archv  file"
#example ./interp_isub_TP5_2_TP6.sh TP5archm.2024_184_12.a
echo "-------------------------------------------------------"
if [ "$#" -ne 1 ]; then
   echo "$0: Usage: $0 TP5 input_archive_file "
   echo " "
   echo "Example ./interp_isub_TP5_2_TP6.sh TP5archv.2024_184_00.a"
   echo " "
   # --- input  TP5  archv  file:
   # --- output TP6  archv  file
   exit 1
fi

TP5_input_file=$1
echo "input file $TP5_input_file"
#TP6_output_file=TP6$(echo $TP5_input_file | rev | cut -c1-19 | rev)
TP6_output_file=$(echo ${TP5_input_file} | sed -e 's/TP5/TP6/g')
echo "TP6_output_file= ${TP6_output_file::-2}.a"
echo ""
echo "------ converting TP5 file to TP6 file------"
mkdir -p "${SCR}"
#
pushd ${SCR}
echo "Entering  --------${SCR}--------"
echo ""
#clean up
/bin/rm -f TP6*.[ab]
/bin/rm -f regional.depth.* regional.grid.*

logfile=Interpolate_tp5.log
echo "ln -sf ${TP6_subregion}/TP6_regional.grid.a TP6_regional.grid.a"
ln -sf ${TP6_subregion}/TP6_regional.grid.a TP6_regional.grid.a
ln -sf ${TP6_subregion}/TP6_regional.grid.b TP6_regional.grid.b
echo "ln -sf ${TP6_subregion}/TP6_regional.depth.a TP6_regional.depth.a"
ln -sf ${TP6_subregion}/TP6_regional.depth.a TP6_regional.depth.a
ln -sf ${TP6_subregion}/TP6_regional.depth.b TP6_regional.depth.b
echo "ln -sf ${TP6_subregion}/TP6_gmap.a TP6_gmap.a"
ln -sf ${TP6_subregion}/TP6_gmap.a TP6_gmap.a
ln -sf ${TP6_subregion}/TP6_gmap.b TP6_gmap.b
#
#
echo "ln -sf ${TP6_subregion}/TP5_regional.depth.a TP5_regional.depth.a"
ln -sf ${TP6_subregion}/TP5_topo/regional.depth.a  regional.depth.a
ln -sf ${TP6_subregion}/TP5_topo/regional.depth.b  regional.depth.b
ln -sf ${TP6_subregion}/TP5_topo/regional.grid.a   regional.grid.a
ln -sf ${TP6_subregion}/TP5_topo/regional.grid.b   regional.grid.b
echo " "
#./isubaregion<<EOF > $logfile
${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL/subregion/src_2.2.12/isubaregion<<EOF 
${SCR}/TP6_regional.grid.a
${SCR}/TP6_gmap.a
${SCR}/TP6_regional.depth.a
${SCR}/regional.depth.a
${SCR}/${TP5_input_file::-2}.a
${SCR}/${TP6_output_file::-2}.a
Interpolating from region TP5 to region TP6
1600   'idm   '
1520   'jdm   '
0	    'iceflg' = ice in output archive flag (0=none,1=energy loan model)
0	    'smooth' = smooth interface depths    (0=F,1=T)
15     'iscan'
EOF
#
echo " ----  ---- "
echo " ----  ---- "
if [ -e ${SCR}/${TP6_output_file::-2}.a ]; then
   #mv ${SCR}/${TP6_output_file::-2}.a ${O}/
   #mv ${SCR}/${TP6_output_file::-2}.b ${O}/
   echo "Successfully created restart file ${SCR}/${TP6_output_file::-2}.a"
else
   echo "Error: something went wrong"
   exit 1
fi
#echo "Exiting  ${SCR}"
#popd
echo "-------------------------------------------------------"
