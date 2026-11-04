#!/bin/bash
#set -x  # Enables debugging

# --- Form a HYCOM restart file from a HYCOM archive file.
RID="TP5a0.06"
R="${HOME}/sea/TOPAZ5/NERSC-HYCOM-CICE/${RID}"
D="${HOME}/work/sea/TOPAZ5/Forecast14/Forecast14_000/"
#process under TP6
O="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/initialize"
SCR="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/initialize/SCRATCH"
#ArcSRC="/${HOME}/sm_alfal/sea/TOPAZ6/NERSC-HYCOM-CICE/TP6a0.03/ConvertTP5"
mkdir -p "${SCR}"
#
# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
   echo "$0: Usage: $0 input_archive_file template_restart_file"
   echo " "
   echo "Example: bash ./topaz_restart2archv_TP5_CML.sh input_restart_file"
   echo " "
   # --- input restart file:
   # --- output archive file
   exit 1
fi

#input_archive_file_in=$1
input_restart_file=$1
#echo "input_archvf= ${input_archive_file_in} " 
echo "input_restatf= ${input_restart_file}"
#

#input_archive_file=$(echo "$input_archive_file_in" | sed -e 's/TP6//g' -e 's/_12/_00/g' )
#input_archive_file=${input_archive_file/archm/archv}
#input_archive_file=${input_archive_file/_12/_00}
#in_restart_file=${input_restart_file/TP6/}
output_archv_file=$(echo ${input_restart_file} | sed -e 's/restart/archv/g' -e 's/_0000//g' )
#output_restart_file=$(echo "$input_archive_file_in" | sed -e 's/archm/restart/g' -e 's/12/00_0000/g'  )
echo "output_archv_file= ${output_archv_file}"
#
#echo "ln -sf ${template_restart_file_in} ${template_restart_file}"
#ln -sf ${template_restart_file_in} ${template_restart_file}
#ln -sf ${template_restart_file_in%.a}.b ${template_restart_file%.a}.b
#

# Create output directory if it does not exist
pushd ${SCR} 
echo "Entering  -----${SCR}------"
#clean dir
/bin/rm -f TP*.[ab]
/bin/rm -f regional.depth.* regional.grid.*

echo "cp ${D}/${input_restart_file}      ${input_restart_file}"
echo "cp ${D}/${input_restart_file%.a}.b ${input_restart_file%.a}.b"
cp ${D}/${input_restart_file}       ${input_restart_file}
cp ${D}/${input_restart_file%.a}.b  ${input_restart_file%.a}.b
#exit
# Handling regional.depth.a and regional.depth.b
touch regional.depth.a regional.depth.b
if [ ! -s regional.depth.a ]; then
  /bin/rm -f regional.depth.a
  /bin/ln -s "${R}/topo/depth_${RID}_05.a" regional.depth.a
fi
if [ ! -s regional.depth.b ]; then
  /bin/rm -f regional.depth.b
  /bin/ln -s "${R}/topo/depth_${RID}_05.b" regional.depth.b
fi
# Handling regional.grid.a and regional.grid.b
touch regional.grid.a regional.grid.b

if [ ! -s regional.grid.a ]; then
  /bin/rm -f regional.grid.a
  /bin/ln -s "${R}/topo/regional.grid.a" .
fi

if [ ! -s regional.grid.b ]; then
  /bin/rm -f regional.grid.b
  /bin/ln -s "${R}/topo/regional.grid.b" .
fi
####
#/bin/rm -f ${D}/051_archvR2.2011_274_00.a ${D}/051_archvR2.2011_274_00.b
#/bin/rm -f ${D}/${in_restart_file} ${D}/${in_restart_file%.a}.b
/bin/rm -f ./${output_archv_file} ./${output_archv_file%.a}.b


${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL/archive/src/restart2archv <<E-o-D
${SCR}/${input_restart_file}
${SCR}/${output_archv_file}
 020   'iexpt ' = experiment number x10 (000=from archive file)
   3    'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
   0	  'sshflg' = diagnostic SSH flag (0=SSH,1=SSH&stericSSH)
   0	  'iceflg' = ice model flag (0=none,1=energy loan model)
 800    'idm   ' = longitudinal array size
 760    'jdm   ' = latitudinal  array size
   0	  'kapref' = thermobaric ref. state (-1=input,0=none,1,2,3=constant)
  50	  'kdm   ' = number of layers
   2	  'n     ' = extract restart time slot number (1 or 2)
  25.0   'thbase' = reference density (sigma units)
   0.1    'sigma ' = layer  1  density (sigma units)
   0.2    'sigma ' = layer  2  density (sigma units)
   0.3    'sigma ' = layer  3  density (sigma units)
   0.4    'sigma ' = layer  4  density (sigma units)
   0.5    'sigma ' = layer  5  density (sigma units)
   0.6    'sigma ' = layer  6  density (sigma units)
   0.7    'sigma ' = layer  7  density (sigma units)
   0.8    'sigma ' = layer  8  density (sigma units)
   0.9    'sigma ' = layer  9  density (sigma units)
   1.0    'sigma ' = layer 10  density (sigma units)
  24.40   'sigma ' = layer 11  density (sigma units)
  25.72   'sigma ' = layer 12  density (sigma units)
  26.40   'sigma ' = layer 13  density (sigma units)
  26.90   'sigma ' = layer 14  density (sigma units)
  27.13   'sigma ' = layer 15  density (sigma units)
  27.19   'sigma ' = layer 16  density (sigma units)
  27.25   'sigma ' = layer 17  density (sigma units)
  27.29   'sigma ' = layer 18  density (sigma units)
  27.34   'sigma ' = layer 19  density (sigma units)
  27.38   'sigma ' = layer 20  density (sigma units)
  27.42   'sigma ' = layer 21  density (sigma units)
  27.50   'sigma ' = layer 22  density (sigma units)
  27.56   'sigma ' = layer 23  density (sigma units)
  27.63   'sigma ' = layer 24  density (sigma units)
  27.66   'sigma ' = layer 25  density (sigma units)
  27.69   'sigma ' = layer 26  density (sigma units)
  27.72   'sigma ' = layer 27  density (sigma units)
  27.74   'sigma ' = layer 28  density (sigma units)
  27.76   'sigma ' = layer 29  density (sigma units)
  27.78   'sigma ' = layer 30  density (sigma units)
  27.80   'sigma ' = layer 31  density (sigma units)
  27.82   'sigma ' = layer 32  density (sigma units)
  27.85   'sigma ' = layer 33  density (sigma units)
  27.87   'sigma ' = layer 34  density (sigma units)
  27.90   'sigma ' = layer 35  density (sigma units)
  27.93   'sigma ' = layer 36  density (sigma units)
  27.95   'sigma ' = layer 37  density (sigma units)
  27.97   'sigma ' = layer 38  density (sigma units)
  27.99   'sigma ' = layer 39  density (sigma units)
  28.01   'sigma ' = layer 40  density (sigma units)
  28.02   'sigma ' = layer 41  density (sigma units)
  28.03   'sigma ' = layer 42  density (sigma units)
  28.04   'sigma ' = layer 43  density (sigma units)
  28.05   'sigma ' = layer 44  density (sigma units)
  28.06   'sigma ' = layer 45  density (sigma units)
  28.07   'sigma ' = layer 46  density (sigma units)
  28.08   'sigma ' = layer 47  density (sigma units)
  28.09   'sigma ' = layer 48  density (sigma units)
  28.10   'sigma ' = layer 49  density (sigma units)
  28.13   'sigma ' = layer 50  density (sigma units)
E-o-D

echo " ----  ---- "
echo " ----  ---- "
echo " ----  ---- "
if [ -e ${SCR}/${output_archv_file} ]; then
   #mv ${SCR}/${output_archv_file}        ${O}/${output_archv_file}
   #mv ${SCR}/${output_archv_file%.a}.b   ${O}/${output_archv_file%.a}.b
   #echo "Successfully created archv file ${O}/${output_archv_file}"
   echo "Successfully created archv file ${SCR}/${output_archv_file}"
else
   echo "Error: something went wrong"
   exit 1
fi

