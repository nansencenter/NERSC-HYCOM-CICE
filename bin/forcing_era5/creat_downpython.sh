
# double-check the directories' existence
# Inidir : the relevant home directory for the script
# Rundir : the working direcotory to download and further processes 

Inidir=$(pwd)
Inidir=/cluster/home/xiejp/work2023/ERA5_down

Rundir=/cluster/work/users/xiejp/ERA5hourly
source ${Inidir}/env.drt


echo ${Rundir}
#[ ! -s ${Rundir} ] && exit ("Cannot find ${Rundir}!")
cd ${Rundir}

Downvars="10m_u_component_of_wind 10m_v_component_of_wind 2m_dewpoint_temperature 2m_temperature boundary_layer_height mean_sea_level_pressure total_cloud_cover sea_ice_cover skin_temperature surface_net_solar_radiation surface_net_thermal_radiation surface_solar_radiation_downwards surface_thermal_radiation_downwards runoff total_precipitation"

Y1=$1
Y2=$2

for iyy in `seq ${Y1} ${Y2}`; do

if [ $# -eq 3 ]; then
   ii=$3
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
  echo $ii ' + ' ${ivarnam} ' - ' ${varname} 
  cat ${Inidir}/Infile/hourly_era5yearly_nc.mal |\
     sed "s#DOWNVAR#${ivarnam}#" |\
     sed "s/VARNAME/${varname}/" |\
     sed "s/YYYY/${iyy}/g" > down_${iyy}_${varname}.py
  #if [ -s down_${iyy}_${varname}.py -a ! -s ./download/hourly.${varname}_${iyy}.nc ]; then
  #   python down_${iyy}_${varname}.py
  #fi

done

done
