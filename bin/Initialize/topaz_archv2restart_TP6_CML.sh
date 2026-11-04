#!/bin/bash
#set -x  # Enables debugging

# --- Form a HYCOM restart file from a HYCOM archive file.
RID="TP6a0.03"
R="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}"
D="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/expt_07.0/data"
O="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize/"
SCR="${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/${RID}/initialize/SCRATCH"
mkdir -p "${SCR}"
#
# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
   echo "$0: Usage: $0 input_archive_file template_restart_file"
   echo " "
   echo "Example: bash ./topaz_archv2restart_TP6_CML.sh input_archive_file template_restart_file"
   echo "Example: bash ./topaz_archv2restart_TP6_CML.sh TP6archv.2024_274_00.a TP6restart.2024_190_00_0000.a"
   echo " "
   # --- input archive file
   # --- input restart file: (template) can be any restart file from any date
   # --- output restart file
   exit 1
fi

input_archive_file_in=$1
template_restart_file_in=$2
echo "input_archvf= ${input_archive_file_in} " 
echo "input_template_restatf= ${template_restart_file_in}"
#
#input_archive_file=$(echo "$input_archive_file_in" | sed -e 's/TP6//g' -e 's/_12/_00/g' )
#input_archive_file=${input_archive_file/archm/archv}
#input_archive_file=${input_archive_file/_12/_00}
#template_restart_file=${template_restart_file_in/TP6/}
# Make linkage after removing TP6
pushd ${SCR} 
echo "Entering  ${SCR}"
#/bin/rm -f TP*.[ab]
/bin/rm -f regional.depth.* regional.grid.*

#AAA    echo "ln -sf ${D}/${input_archive_file_in} ${input_archive_file_in}"
#AAA    ln -sf ${D}/${input_archive_file_in}      ${input_archive_file_in}
#AAA    ln -sf ${D}/${input_archive_file_in%.a}.b ${input_archive_file_in%.a}.b
#
#output_restart_file="${input_archive_file/archm/restart}"
output_restart_file=$(echo "$input_archive_file_in" | sed -e 's/archv/restart/g' -e 's/_00_/_00_0000/g'  )
#output_restart_file=$(echo "$output_restart_file" | sed -e 's/archv/restart/g'  )
#
echo "ln -sf ${O}/${template_restart_file_in} ${template_restart_file_in}"
ln -sf ${O}/${template_restart_file_in}       ${template_restart_file_in}
ln -sf ${O}/${template_restart_file_in%.a}.b  ${template_restart_file_in%.a}.b
#

# Create output directory if it does not exist

# Handling regional.depth.a and regional.depth.b
touch regional.depth.a regional.depth.b

if [ ! -s regional.depth.a ]; then
  /bin/rm -f regional.depth.a
  /bin/ln -s "${R}/topo/depth_${RID}_08.a" regional.depth.a
fi

if [ ! -s regional.depth.b ]; then
  /bin/rm -f regional.depth.b
  /bin/ln -s "${R}/topo/depth_${RID}_08.b" regional.depth.b
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

/bin/rm  -f ./${output_restart_file} ./${output_restart_file%.a}.b
# Call to the external executable (archv2restart)
# Assuming you want to execute a command like in the original script
${HOME}/sea/TOPAZ6/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL/archive/src/archv2restart <<E-o-D
${SCR}/${input_archive_file_in}
${SCR}/${template_restart_file_in}
${SCR}/${output_restart_file}
 070     'iexpt ' = experiment number x10 (000=from archive file)
  3     'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
1600     'idm   ' = longitudinal array size
1520     'jdm   ' = latitudinal  array size
  0     'kapref' = thermobaric ref. state (-1=input,0=none,1,2,3=constant)
 50     'kdm   ' = number of layers
 25.0   'thbase' = reference density (sigma units)
180.0   'baclin' = baroclinic time step (seconds), int. divisor of 86400
  0     'rmontg' = pbavg correction from relax.montg file  (0=F,1=T)
E-o-D
echo " ----  ---- "
echo " ----  ---- "
echo " ----  ---- "
if [ -e ${SCR}/${output_restart_file} ]; then
  echo "Successfully created restart file ${SCR}/${output_restart_file} "
  #rename output file to remove mem000 if exists
  newtp6_name=$(echo "${output_restart_file}" | sed -e 's/mem000//g' )
  echo "Renaming output file to ${newtp6_name} "
  echo "mv ${output_restart_file} ${newtp6_name}"
  mv ${output_restart_file}      ${newtp6_name}
  mv ${output_restart_file%.a}.b ${newtp6_name%.a}.b
  echo "Successfully created restart file ${SCR}/${newtp6_name}"
else
   echo "Error: something went wrong"
   exit 1
fi

#echo "Exiting  ${SCR}"
#popd

