#! /bin/bash
module load CDO/2.2.2-gompi-2023b

source_nesting_experiment_path=$PWD
source ../REGION.src
#
echo ''
echo 'The script will look for decadal files under '${Nesting_Files_PATH}'.'
echo 'Modify REGION.src if necessary'
echo 'The script assumes you are running in an experiment folder located in the Nesting Region (e.g. ESMa1.00/expt_01.0)'
echo ''
#
# Ensure Nesting_Files_PATH is set
if [ -z "$Nesting_Files_PATH" ]; then
    echo "Error: Nesting_Files_PATH is not set."
    echo "Check issues regarding REGION.src"
    exit 1
fi

long_name=$1

# Default WOA directory
directory="/cluster/projects/nn2993k/ModelInput/WOA2018"
# you can overwrite this by:
# e.g. ../bin/Nesting_noresm/Create_climatology_for_bias_correction.sh silicate woa=/path/to/other/directory

# Parse arguments
for arg in "$@"; do
    case $arg in
        woa=*)
            directory="${arg#*=}"
            shift
            ;;
        *)
            long_name="$arg"
            ;;
    esac
done

case "$long_name" in
    temperature)
        vari="t"
        version="decav"
        resol="04"
        esmvar="thetao"
        ;;
    salinity)
        vari="s"
        version="decav"
        resol="04"
        esmvar="so"
        ;;
    nitrate)
        vari="n"
        version="all"
        resol="01"
        esmvar="no3"
        ;;
    phosphate)
        vari="p"
        version="all"
        resol="01"
        esmvar="po4"
        ;;
    silicate)
        vari="i"
        version="all"
        resol="01"
        esmvar="si"
        ;;
    oxygen)
        vari="o"
        version="all"
        resol="01"
        esmvar="o2"
        ;;
    *)
        echo "Error: long_name '$long_name' is not recognized."
        exit 1
        ;;
esac

echo "vari=$vari, version=$version, resol=$resol, esmvar=$esmvar"
echo "Using WOA directory: $directory"
echo ''

sleep 5 # give the user chance to control

# select analyzed mean from all the files
echo "select analyzed mean of" ${long_name} "from all the files in ${directory}"  
for ((mon=00; mon<=16;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon 
   cdo selvar,${vari}_an "${directory}/woa18_${version}_${vari}${smon}_${resol}.nc" \
                         "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_${resol}.nc"
done


# regrid to NORESM
echo "regrid to NORESM"
for ((mon=00; mon<=16;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo remapbil,${Nesting_Files_PATH}/cdogrid_NESM_T "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_${resol}.nc" \
                               "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm.nc"
done


# select the deep levels from the seasonal files for t, s and oxygen (monthly WOA fields only
# go to 1500 m/57 levels, seasonal fields go to the full 5500 m/102 levels), and from the
# annual file for nitrate/phosphate/silicate (monthly fields only go to 800 m/43 levels)
echo "select the deep levels from the seasonal files"
if [ "$vari" = "t" -o  "$vari" = "s" -o "$vari" = "o" ]; then
    mon1=13; mon2=16; numl=58; #temperature, salinity or oxygen
    winnum=13; sprnum=14; sumnum=15; autnum=16; #season files in the deep
else
    mon1=0; mon2=0;   numl=44; # nutrients
    winnum=00; sprnum=00; sumnum=00; autnum=00; #annual file in the deep
fi

for ((mon=$mon1; mon<=$mon2;mon+=1)); do 
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo sellevidx,$numl/102 "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm.nc"\
                        "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deep.nc"
done

# merge the seasonal files to the monthly files
echo "merge the seasonal files to the monthly files"
for ((mon=1; mon<=3;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${winnum}_noresm_deep.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deep.nc"
done

for ((mon=4; mon<=6;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${sprnum}_noresm_deep.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deep.nc"
done

for ((mon=7; mon<=9;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${sumnum}_noresm_deep.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deep.nc"
done

for ((mon=10; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo merge "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${autnum}_noresm_deep.nc" \
             "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deep.nc"
done

# correct the level bounds
echo "correct the level bounds"
for ((mon=1; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo genlevelbounds "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deep.nc" \
                      "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deepb.nc"
done

# regrid to the NORESM depth levels
echo "regrid to the NORESM depth levels"
for ((mon=1; mon<=12;mon+=1)); do
   smon=`echo -n 0$mon | tail -2c`
   echo $smon
   cdo intlevel,0,5,10,15,20,25,30,40,50,62.5,75,87.5,100,112.5,125,\
137.5,150,175,200,225,250,275,300,350,400,450,500,550,600,\
650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,\
1300,1350,1400,1450,1500,1625,1750,1875,2000,2250,2500,2750,\
3000,3250,3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,\
6000,6250,6500,6750\
                "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_deepb.nc" \
                "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_${smon}_noresm_esmlev.nc"
done

# merge to one file
echo "merge to one file"
cdo mergetime "${Nesting_Files_PATH}/woa18_${version}_${vari}_an_*_noresm_esmlev.nc"  \
   "${Nesting_Files_PATH}/woa18_${version}_${vari}_year_noresm_esmlev.nc"

# convert climatology to same unit as the ESM
if [ "$vari" = "n" -o  "$vari" = "i" -o "$vari" = "p" -o "$vari" = "o" ]; then
    # from microm/kg to mole/m3
    # assming constant seawter density of 1028 kg/m3.
    cdo mulc,1.028e-3 "${Nesting_Files_PATH}/woa18_${version}_${vari}_year_noresm_esmlev.nc"\
	"${Nesting_Files_PATH}/woa18_${version}_${vari}_year_noresm_esmlev_unit.nc"
else
    #already the same unit
    mv "${Nesting_Files_PATH}/woa18_${version}_${vari}_year_noresm_esmlev.nc"\
       "${Nesting_Files_PATH}/woa18_${version}_${vari}_year_noresm_esmlev_unit.nc"
fi
# create the bias file
echo "create the bias file"
cdo sub "${Nesting_Files_PATH}/${esmvar}_Omon_NorESM2-MM_historical_r1i1p1f1_gr_clim.nc" \
        "${Nesting_Files_PATH}/woa18_${version}_${vari}_year_noresm_esmlev_unit.nc" \
        "${Nesting_Files_PATH}/bias_${esmvar}_decal.nc"



# make monthly files
cdo splitmon "${Nesting_Files_PATH}/bias_${esmvar}_decal.nc" "${Nesting_Files_PATH}/bias_${esmvar}_decal_"

