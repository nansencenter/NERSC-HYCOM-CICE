#! /bin/bash -l
#
if [ -s env.drt ]; then
   source env.drt
else
   Inidir=$(pwd)
   Rundir=/cluster/work/users/xiejp/ERA5hourly
fi    

[ ! -s ${Rundir} ] && { echo "Cannot find ${Rundir}!";  exit 1; }
cd ${Rundir}

Downvars="10m_u_component_of_wind 10m_v_component_of_wind 2m_dewpoint_temperature 2m_temperature boundary_layer_height mean_sea_level_pressure total_cloud_cover sea_ice_cover skin_temperature surface_net_solar_radiation surface_net_thermal_radiation surface_solar_radiation_downwards surface_thermal_radiation_downwards runoff total_precipitation"

Y1=$1
Y2=$2
if [ $# -eq 3 ]; then
   MM=$3
else
   MM=1
fi

if [ $MM -gt 0 ]; then
   ${Inidir}/creat_mmdown.sh ${Y1} ${Y2} ${MM}
else
   ${Inidir}/creat_downpython.sh ${Y1} ${Y2}
fi

echo ${Rundir}
#[ ! -s ${Rundir} ] && exit ("Cannot find ${Rundir}!")
cd ${Rundir}
[ ! -r ./download ] && mkdir download 

for iyy in `seq ${Y1} ${Y2}`; do
   if [ $# -gt 3 ]; then
      # point out the variable inside of the 15 variabls
      ii=$4
   else
      ii=0
   fi
   for ivarnam in ${Downvars} ; do
       (( ii = ii + 1 ))
       case ${ii} in
       1) 
         varname='10U'  
          ;;
       2)
         varname='10V' 
          ;;
       3)
         varname='2D'  
          ;;
       4)
         varname='2T'  
          ;;
       5)
         varname='BLH'  
          ;;
       6)
         varname='MSL'  
          ;;
       7)
         varname='TCC'  
          ;;
       8)
         varname='CI'  
         ;;
       9)
         varname='SKT'  
          ;;
      10)
         varname='SSR'  
          ;;
      11)
         varname='STR'  
          ;;
      12)
         varname='SSRD'  
          ;;
      13)
         varname='STRD'  
          ;;
      14)
         varname='RO'  
          ;;
      15)
         varname='TP'  
          ;;
      esac 

      if [ -s down_${iyy}_${varname}.py -a ! -s ./download/hourly.${varname}_${iyy}.nc ]; then
         python down_${iyy}_${varname}.py  
         [ -s hourly.${varname}_${iyy}.nc ] && mv hourly.${varname}_${iyy}.nc ./download/.
      fi
done
done
echo "Download done!"
exit $?

