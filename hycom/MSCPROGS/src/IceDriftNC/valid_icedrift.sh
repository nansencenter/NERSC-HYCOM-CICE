
#
# Assumes ice drift from OSISAF is accessible as bellow:
Obsdir='/cluster/projects/nn2993k/TP4b0.12/idrft_osisaf/'

# Model grid file and grid.info should be available:
#ln -sf ${Moddir}/regional.*.? .
#ln -sf ${Moddir}/grid.inof .

# Assumes ice drift routine is in same directory as this script
prog0=$(dirname $0)
prog=$prog0/icedriftnc_osisaf
iy=2019
im='02'
idd=10
#for idd in `seq 1 31`; do
   cdd=`echo 0${idd}|tail -3c`
   Fmod=iceh.${iy}-${im}-${cdd}.nc   
   if [ -s ${Fmod} ]; then
      ${prog} ${Fmod} ${Obsdir} 
      if [ -s Mod_drift000.nc -a -s Obs_drift.txt ]; then
         ml load matplotlib/3.2.1-intel-2020a-Python-3.8.2
         python Firstview_icespeed.py    
      fi
   fi
#done
