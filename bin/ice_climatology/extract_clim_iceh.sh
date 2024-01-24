###############################################################################
# Usage:
#       extract_clim_iceh.sh <modgridpath> 
#  or
#       extract_clim_iceh.sh 
#
# It will create an ice climatology file named ice_initial.nc 
# which should be under $model/scratch direcotry because the hycom_cice will access this file.
# (In fact, its location is predefined by ocn_data_dir in ice_in
# 
# The later will use the default run of the test in TP5
# Based on the sea ice climatology from the TP4b reanalysis in 1991-2020,
# interpolating the sea ice, SST, and SSS on the target model grids
# 
# Basic software requires:
# 1) python 3; cdo and nco modules; 
#    For example on Betzy:
#    ml load CDO/1.9.9-iompi-2020a
#    ml load NCO/4.9.7-iomkl-2020a
# 2) target model grid information and ice grid files;
# 3) link the ice_kmt.nc to the present work directory;
# 4) the iceh climatology file as the source file for interpolation.
# 
## Created by JX at 7th December 2022 on Betzy 
#
###############################################################################
# First check the relevant model grid files:
#
echo "It requires the model grid files to be available!"
if [ $# -eq 1 ]; then
   Targetmoddir=$1
else
   echo "Default run in a test case in TP5"
   Targetmoddir='/cluster/work/users/xiejp/TOPAZ/TP5a0.06/topo'
fi
files1="regional.grid regional.depth"
i0=0
for ifile in ${files1}; do
   file0=${Targetmoddir}/${ifile}
   if [ ! -s ${file0}.a -o ! -s ${file0}.b ]; then
      (( i0 = i0 + 1 ))
   fi
done

# sea ice climatology files:
Fini=TP4b_1991-2020_AssimSurf.nc

Fraw=cice_kmt.nc
files2="${Fini} ${Fraw} createmask.py"
j0=0
for ifile in ${files2}; do
   if [ ! -s ${ifile} ]; then
      (( j0 = j0 + 1 ))
   fi
done
if [ $i0 -gt 0 ]; then
   echo "Check the input model diretory: " ${Targetmoddir}
   exit 0
fi
if [ $j0 -gt 0 ]; then
   echo "Check the source files: ${files2} under "$(pwd) 
   exit 1
fi

# Second prepare the interpolation mask:
# python createmask.py /cluster/work/users/xiejp/TOPAZ/TP5a0.06/topo/ testmask
echo "Creating the mask for interpolation ... "
Ftmpmask=inter_mask
[ -s ${Ftmpmask} ] && rm ${Ftmpmask}
python createmask.py ${Targetmoddir} ${Ftmpmask}
if [ ! -s ${Ftmpmask} ]; then
   echo "The interpolation mask is not available "${Ftmpmask} 
   exit 1
fi


# Third interpolating
Vars="siconc,sithick,thetao,so"

Fnew=ice_initial.nc
[ -s ${Fnew} ] && rm ${Fnew}

Fout=ice_0.nc
Ftmp=ice_1.nc
Ftmp2=ice_2.nc

Isst=0   # Switch on sss/sst processing ofr Isst=1,otherwise switch off

if [ -s ${Fini} -a -s ${Ftmpmask} ]; then
   [ -s ${Fout} ] && rm ${Fout}
   echo "interpolating ... "
   module load CDO/1.9.9-iompi-2020a
   cdo remapbil,${Ftmpmask} ${Fini} ${Fout}

   ml load NCO/4.9.7-iomkl-2020a
   echo "Defaulting the values for different masks ... "
   ncks -v ${Vars} ${Fout} ${Ftmp}  
   ncrename -h -O -v siconc,aice_raw -v sithick,hi_raw ${Ftmp}  
   ncrename -h -O -v so,sss_raw -v thetao,sst_raw ${Ftmp}  

   if [ -s cice_kmt.nc ]; then
     ncks -4 cice_kmt.nc kmt4.nc
     ncks -A -v kmt kmt4.nc ${Ftmp} 
     [ -s kmt4.nc ] && rm kmt4.nc
     [ -s ${Ftmp2} ] && rm ${Ftmp2}
     if [ ${Isst} -eq 1 ]; then
        ncap2 -F -s "hi=kmt*0.0; aice=kmt*0.0; sst=kmt*0.0; sss=kmt*0.0" ${Ftmp} ${Ftmp2}
     else
        ncap2 -F -s "hi=kmt*0.0; aice=kmt*0.0" ${Ftmp} ${Ftmp2}
     fi
     [ -s ${Ftmp} ] && rm ${Ftmp}
     ncap2 -F -s "where(aice_raw>0) hi=hi+hi_raw; where(kmt==0) hi=hi_raw.get_miss()" ${Ftmp2} ${Ftmp}
     [ -s ${Ftmp2} ] && rm ${Ftmp2}

     if [ ${Isst} -eq 1 ]; then
        echo "sss/sst process needs to improve for the default missing values ... "
        ncap2 -F -s "where(kmt==1) sst=0.0; where(kmt==1&&sst_raw>-2.) sst=sst_raw; where(kmt==0) sst=sst_raw.get_miss()" ${Ftmp} ${Ftmp2}
        [ -s ${Ftmp} ] && rm ${Ftmp}
        ncap2 -F -s "where(kmt==1) sss=0.0; where(kmt==1&&sss_raw>0.0) sss=sss_raw; where(kmt==0) sss=sss_raw.get_miss()" ${Ftmp2} ${Ftmp}
        [ -s ${Ftmp2} ] && rm ${Ftmp2}
     fi

     if [ -s ${Ftmp}  ]; then
       ncap2 -F -s "where(aice_raw>0) aice=aice+aice_raw; where(kmt==0) aice=aice_raw.get_miss()" ${Ftmp} ${Ftmp2}
       # cutting the tiny sea ice dots
       ncap2 -F -s "where(aice_raw<.15) aice=0.0; where(aice_raw<.15) hi=0.0" ${Ftmp2} ${Fnew}
       rm ice_?.nc inter_mask
     fi
   fi 
fi
