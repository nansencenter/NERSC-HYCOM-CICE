#! /bin/bash
module load CDO/2.2.2-gompi-2023b

source_nesting_experiment_path=$PWD
source ../REGION.src
#
echo ''
echo 'Processing GLODAP climatology for bias correction (dissic or talk).'
echo 'The script assumes you are running in an experiment folder located in the Nesting Region (e.g. ESMa1.00/expt_01.0)'
echo 'Example usage: ../bin/Nesting_noresm/Create_climatology_for_bias_correction_dissic_talk.sh dissic'
echo ''
#
# Ensure Nesting_Files_PATH is set
if [ -z "$Nesting_Files_PATH" ]; then
    echo "Error: Nesting_Files_PATH is not set."
    echo "Check issues regarding REGION.src"
    exit 1
fi

long_name=$1

# Default GLODAP directory
directory="/cluster/projects/nn2993k/ModelInput/GLODAPV2/"
# you can overwrite this by:
# e.g. ../bin/Nesting_noresm/Create_climatology_for_bias_correction_GLODAP.sh dissic glodap=/path/to/other/directory

# Parse arguments
for arg in "$@"; do
    case $arg in
        glodap=*)
            directory="${arg#*=}"
            shift
            ;;
        *)
            long_name="$arg"
            ;;
    esac
done

case "$long_name" in
    dissic)
        glodapvar="TCO2"
        glodapfile="GLODAPv2.2016b.TCO2.nc"
        esmvar="dissic"
        ;;
    talk)
        glodapvar="TAlk"
        glodapfile="GLODAPv2.2016b.TAlk.nc"
        esmvar="talk"
        ;;
    *)
        echo "Error: long_name '$long_name' is not recognized."
        echo "Valid options: dissic, talk"
        exit 1
        ;;
esac

echo "glodapvar=$glodapvar, esmvar=$esmvar"
echo "Using GLODAP directory: $directory"
echo "Input file: ${directory}/${glodapfile}"
echo ''

sleep 5 # give the user a chance to stop

# Extract the relevant variable from the GLODAP annual mean file
echo "Extracting ${glodapvar} from ${glodapfile}"
cdo selvar,${glodapvar} "${directory}/${glodapfile}" \
    "${Nesting_Files_PATH}/glodap_${esmvar}_annual.nc"


# GLODAP only provides an annual mean — no seasonal variability.
# Copy the annual mean to 12 monthly files, each stamped at the 15th of the month.
echo "Creating 12 monthly copies of the annual mean"
for ((mon=1; mon<=12; mon+=1)); do
    smon=$(printf "%02d" $mon)
    echo $smon
    cdo settaxis,2002,${mon},15 \
        "${Nesting_Files_PATH}/glodap_${esmvar}_annual.nc" \
        "${Nesting_Files_PATH}/glodap_${esmvar}_${smon}.nc"
done


# Regrid to NORESM horizontal grid
echo "Regridding to NORESM horizontal grid"
for ((mon=1; mon<=12; mon+=1)); do
    smon=$(printf "%02d" $mon)
    echo $smon
    cdo remapbil,${Nesting_Files_PATH}/cdogrid_NESM_T \
        "${Nesting_Files_PATH}/glodap_${esmvar}_${smon}.nc" \
        "${Nesting_Files_PATH}/glodap_${esmvar}_${smon}_noresm.nc"
done


# Correct the level bounds
# (GLODAP has 33 depth surfaces; genlevelbounds adds the bounds needed by intlevel)
echo "Correcting level bounds"
for ((mon=1; mon<=12; mon+=1)); do
    smon=$(printf "%02d" $mon)
    echo $smon
    cdo genlevelbounds \
        "${Nesting_Files_PATH}/glodap_${esmvar}_${smon}_noresm.nc" \
        "${Nesting_Files_PATH}/glodap_${esmvar}_${smon}_noresm_deepb.nc"
done


# Regrid to NORESM depth levels
# GLODAP reaches 5500 m; CDO will extrapolate with the deepest value for levels below that.
echo "Regridding to NORESM depth levels"
for ((mon=1; mon<=12; mon+=1)); do
    smon=$(printf "%02d" $mon)
    echo $smon
    cdo intlevel,0,5,10,15,20,25,30,40,50,62.5,75,87.5,100,112.5,125,\
137.5,150,175,200,225,250,275,300,350,400,450,500,550,600,\
650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,\
1300,1350,1400,1450,1500,1625,1750,1875,2000,2250,2500,2750,\
3000,3250,3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,\
6000,6250,6500,6750\
        "${Nesting_Files_PATH}/glodap_${esmvar}_${smon}_noresm_deepb.nc" \
        "${Nesting_Files_PATH}/glodap_${esmvar}_${smon}_noresm_esmlev.nc"
done


# Merge to one file
echo "Merging monthly files to one annual file"
cdo mergetime "${Nesting_Files_PATH}/glodap_${esmvar}_"??_noresm_esmlev.nc \
    "${Nesting_Files_PATH}/glodap_${esmvar}_year_noresm_esmlev.nc"


# Convert units from micro-mol/kg to mol/m3
# mol/m3 = micro-mol/kg * density(kg/m3) / 1e6 = * 1.028e-3 (assuming density = 1028 kg/m3)
echo "Converting units: micro-mol/kg -> mol/m3"
cdo mulc,1.028e-3 \
    "${Nesting_Files_PATH}/glodap_${esmvar}_year_noresm_esmlev.nc" \
    "${Nesting_Files_PATH}/glodap_${esmvar}_year_noresm_esmlev_unit.nc"


# Create the bias file
echo "Creating bias file"
cdo sub "${Nesting_Files_PATH}/${esmvar}_Omon_NorESM2-MM_historical_r1i1p1f1_gr_clim.nc" \
        "${Nesting_Files_PATH}/glodap_${esmvar}_year_noresm_esmlev_unit.nc" \
        "${Nesting_Files_PATH}/bias_${esmvar}_decal.nc"


# Split into monthly files
cdo splitmon "${Nesting_Files_PATH}/bias_${esmvar}_decal.nc" \
    "${Nesting_Files_PATH}/bias_${esmvar}_decal_"
