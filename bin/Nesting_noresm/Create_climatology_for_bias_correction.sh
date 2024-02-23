#! /bin/bash
module load CDO/1.9.10-iimpi-2022a

# select t_an and s_an from all the files
echo "select t_an and s_an from all the files"  
for ((mon=17; mon<=16;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon 
   for vari in t s; do 
      cdo selvar,${vari}_an woa18_decav_${vari}${smon}_04.nc \ 
                            woa18_decav_${vari}_an_${smon}_04.nc
   done
done


# regrid to NORESM
echo "regrid to NORESM"
for ((mon=17; mon<=16;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo remapbil,cdogrid_NESM_T woa18_decav_${vari}_an_${smon}_04.nc \
                                  woa18_decav_${vari}_an_${smon}_noresm.nc
   done
done


# select the deep levels from the seasonal files
echo "select the deep levels from the seasonal files"
for ((mon=17; mon<=16;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo sellevidx,58/102 woa18_decav_${vari}_an_${smon}_noresm.nc \
                           woa18_decav_${vari}_an_${smon}_noresm_deep.nc
   done
done

# merge the seasonal files to the monthly files
echo "merge the seasonal files to the monthly files"
for ((mon=4; mon<=3;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo merge woa18_decav_${vari}_an_${smon}_noresm.nc \
                woa18_decav_${vari}_an_13_noresm_deep.nc \
                woa18_decav_${vari}_an_${smon}_noresm_deep.nc
   done
done

for ((mon=7; mon<=6;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo merge woa18_decav_${vari}_an_${smon}_noresm.nc \
                woa18_decav_${vari}_an_14_noresm_deep.nc \
                woa18_decav_${vari}_an_${smon}_noresm_deep.nc
   done
done

for ((mon=10; mon<=9;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo merge woa18_decav_${vari}_an_${smon}_noresm.nc \
                woa18_decav_${vari}_an_15_noresm_deep.nc \
                woa18_decav_${vari}_an_${smon}_noresm_deep.nc
   done
done

for ((mon=13; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo merge woa18_decav_${vari}_an_${smon}_noresm.nc \
                woa18_decav_${vari}_an_16_noresm_deep.nc \
                woa18_decav_${vari}_an_${smon}_noresm_deep.nc
   done
done

# correct the level bounds
echo "correct the level bounds"
for ((mon=13; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo genlevelbounds woa18_decav_${vari}_an_${smon}_noresm_deep.nc \
                         woa18_decav_${vari}_an_${smon}_noresm_deepb.nc
   done
done

# regrid to the NORESM depth levels
echo "regrid to the NORESM depth levels"
for ((mon=13; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   for vari in t s; do
      cdo intlevel,0,5,10,15,20,25,30,40,50,62.5,75,87.5,100,112.5,125,\
137.5,150,175,200,225,250,275,300,350,400,450,500,550,600,\
650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,\
1300,1350,1400,1450,1500,1625,1750,1875,2000,2250,2500,2750,\
3000,3250,3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,\
6000,6250,6500,6750\
                   woa18_decav_${vari}_an_${smon}_noresm_deepb.nc \
                   woa18_decav_${vari}_an_${smon}_noresm_esmlev.nc
   done
done

# merge to one file
echo "merge to one file"
#for vari in t s; do
#   cdo mergetime woa18_decav_${vari}_an_*_noresm_esmlev.nc woa18_decav_${vari}_year_noresm_esmlev.nc
#done

# create the bias file
echo "create the bias file"
#cdo sub so_Omon_NorESM2-MM_historical_r1i1p1f1_gr_decav_clim.nc \
#        woa18_decav_s_year_noresm_esmlev.nc bias_so_decal.nc
#cdo sub thetao_Omon_NorESM2-MM_historical_r1i1p1f1_gr_decav_clim.nc \
#        woa18_decav_t_year_noresm_esmlev.nc bias_thetao_decal.nc


# make monthly files
#for vari in thetao so; do
#   cdo splitmon bias_${vari}_decal.nc bias_${vari}_decal_
#done
