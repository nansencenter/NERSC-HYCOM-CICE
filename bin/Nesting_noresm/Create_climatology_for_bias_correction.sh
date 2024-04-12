#! /bin/bash
module load CDO/1.9.10-iimpi-2022a

#vari=t; long_name=temperature; version=decav; resol=04; esmvar=thetao;
#vari=s; long_name=salinity;    version=decav; resol=04; esmvar=so;
#vari=n; long_name=nitrate;     version=all;   resol=01; esmvar=no3;
#vari=p; long_name=phosphate;   version=all;   resol=01; esmvar=po4;
vari=i; long_name=silicate;    version=all;   resol=01; esmvar=si;
#vari=o; long_name=oxygen;      version=all;   resol=01; esmvar=o2;

# select analyzed mean from all the files
echo "select analyzed mean of" ${long_name} "from all the files"  
for ((mon=23; mon<=16;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon 
   cdo selvar,${vari}_an woa18_${version}_${vari}${smon}_${resol}.nc \
                         woa18_${version}_${vari}_an_${smon}_${resol}.nc
done


# regrid to NORESM
echo "regrid to NORESM"
for ((mon=23; mon<=16;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo remapbil,cdogrid_NESM_T woa18_${version}_${vari}_an_${smon}_${resol}.nc \
                               woa18_${version}_${vari}_an_${smon}_noresm.nc
done


# select the deep levels from the seasonal files for t and s and annual files for nutrients/oxygen
echo "select the deep levels from the seasonal files"
if [ "$vari" = "t" -o  "$vari" = "s" -o "$vari" = "o" ]; then
    mon1=23; mon2=16; numl=58; #temperature, salinity or oxygen
    winnum=13; sprnum=14; sumnum=15; autnum=16; #season files in the deep
else
    mon1=20; mon2=0;   numl=44; # nutrients
    winnum=00; sprnum=00; sumnum=00; autnum=00; #annual file in the deep    
fi

for ((mon=$mon1; mon<=$mon2;mon+=1)); do 
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo sellevidx,$numl/102 woa18_${version}_${vari}_an_${smon}_noresm.nc\
                        woa18_${version}_${vari}_an_${smon}_noresm_deep.nc
done

# merge the seasonal files to the monthly files
echo "merge the seasonal files to the monthly files"
for ((mon=13; mon<=3;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge woa18_${version}_${vari}_an_${smon}_noresm.nc \
             woa18_${version}_${vari}_an_${winnum}_noresm_deep.nc \
             woa18_${version}_${vari}_an_${smon}_noresm_deep.nc
done

for ((mon=46; mon<=6;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge woa18_${version}_${vari}_an_${smon}_noresm.nc \
             woa18_${version}_${vari}_an_${sprnum}_noresm_deep.nc \
             woa18_${version}_${vari}_an_${smon}_noresm_deep.nc
done

for ((mon=79; mon<=9;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge woa18_${version}_${vari}_an_${smon}_noresm.nc \
             woa18_${version}_${vari}_an_${sumnum}_noresm_deep.nc \
             woa18_${version}_${vari}_an_${smon}_noresm_deep.nc
done

for ((mon=100; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge woa18_${version}_${vari}_an_${smon}_noresm.nc \
             woa18_${version}_${vari}_an_${autnum}_noresm_deep.nc \
             woa18_${version}_${vari}_an_${smon}_noresm_deep.nc
done

# correct the level bounds
echo "correct the level bounds"
for ((mon=21; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo genlevelbounds woa18_${version}_${vari}_an_${smon}_noresm_deep.nc \
                      woa18_${version}_${vari}_an_${smon}_noresm_deepb.nc
done

# regrid to the NORESM depth levels
echo "regrid to the NORESM depth levels"
for ((mon=21; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo intlevel,0,5,10,15,20,25,30,40,50,62.5,75,87.5,100,112.5,125,\
137.5,150,175,200,225,250,275,300,350,400,450,500,550,600,\
650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,\
1300,1350,1400,1450,1500,1625,1750,1875,2000,2250,2500,2750,\
3000,3250,3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,\
6000,6250,6500,6750\
                woa18_${version}_${vari}_an_${smon}_noresm_deepb.nc \
                woa18_${version}_${vari}_an_${smon}_noresm_esmlev.nc
done

# merge to one file
echo "merge to one file"
#cdo mergetime woa18_${version}_${vari}_an_*_noresm_esmlev.nc woa18_${version}_${vari}_year_noresm_esmlev.nc

# convert climatology to same unit as the ESM
if [ "$vari" = "n" -o  "$vari" = "i" -o "$vari" = "p" -o "$vari" = "o" ]; then
    # from microm/kg to mole/m3
    # assming constant seawter density of 1028 kg/m3.
    cdo divc,1028.0 woa18_${version}_${vari}_year_noresm_esmlev.nc\
	woa18_${version}_${vari}_year_noresm_esmlev_unit.nc
else
    #already the same unit
    mv woa18_${version}_${vari}_year_noresm_esmlev.nc\
       woa18_${version}_${vari}_year_noresm_esmlev_unit.nc
fi
# create the bias file
echo "create the bias file"
cdo sub ${esmvar}_Omon_NorESM2-MM_historical_r1i1p1f1_gr_clim.nc \
        woa18_${version}_${vari}_year_noresm_esmlev_unit.nc bias_${esmvar}_decal.nc



# make monthly files
cdo splitmon bias_${esmvar}_decal.nc bias_${esmvar}_decal_

