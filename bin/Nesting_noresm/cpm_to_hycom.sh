#!/bin/bash
# Usage: ../bin/Nesting_noresm/cpm_to_hycom.sh /path/to/nesting_expt/ /path/to/preprocessed_dir/*merged_*.nc

options=$(getopt -o m -- "$@")
[ $? -eq 0 ] || { echo "Error: Incorrect options"; exit 1; }
eval set -- "$options"
while true; do
    case "$1" in
        -m) shift ;; # Kept for backward compatibility if flags are used
        --) shift; break ;;
    esac
    shift
done

if [ $# -lt 2 ] ; then
    echo "Usage: $0 [nest_expt_path] [merged_nc_files...]"
    exit 1
fi

if [ -f EXPT.src ] ; then
    export BASEDIR=$(cd .. && pwd)
else
    echo "Error: Run this script inside your active experiment directory containing EXPT.src"
    exit 1
fi

source ${BASEDIR}/REGION.src || { echo "Could not source REGION.src"; exit 1; }
source ./EXPT.src || { echo "Could not source EXPT.src"; exit 1; }

export nest_expt=$1
shift
export input_files="$@"
export esm_gridfile=${ESM_NATIVE_MESH}

iexpt=01        
iversn=22
yrflag=3

function model_datetime() {
python - <<END
from netCDF4 import Dataset, num2date
with Dataset("$1", "r") as ncid:
    time = ncid.variables["time"][0]
    tunit = ncid.variables["time"].units                                                                                                                                                                                                            
    t_cal = ncid.variables["time"].calendar if "calendar" in ncid.variables["time"].ncattrs() else "standard"
date = num2date(time, units=tunit, calendar=t_cal)
print(date.strftime("archv.%Y_%j_00"))
END
}

echo "=== Processing Unified NorCPM Datasets ==="
for file in ${input_files}; do
    echo "Processing: ${file}"
    
    # 1. Run Python converter on the single unified file
    python ${BASEDIR}/bin/Nesting_noresm/cpm2archvz.py "${esm_gridfile}" "${file}" \
         --iexpt ${iexpt} --iversn ${iversn} --yrflag ${yrflag} 
          
    # 2. Get target date layout name
    archv_name=$(model_datetime "${file}")
        
    # 3. Call downstream HYCOM interpolator tool
    ${BASEDIR}/bin/archvz2hycom_biophys.sh ${nest_expt} ${archv_name} -n 70
done
echo "=== Processing Pipeline Finished ==="
