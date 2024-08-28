#!/bin/bash -l

#SBATCH --account=nn9481k
#SBATCH --job-name=gennest
#SBATCH -t 0:30:00
#SBATCH --qos=devel
#SBATCH --nodes=1   # number of nodes
#SBATCH  --mail-type=ALL

#SBATCH -o log/gnnjob.out
#SBATCH -e log/gnnjob.err

##SBATCH --array=1-10
##
##  Give the job a name
##
##         
##module restore system

#export SLURM_SUBMIT_DIR=$(pwd)

#cd $SLURM_SUBMIT_DIR

#year=$((2090+$SLURM_ARRAY_TASK_ID))

year=$1
target_region=$2 # target region (e.g. TP2a0.10) 
target_experiment=$3 # target experiment (e.g. 010)
experiment_full_name="expt_${target_experiment:0:2}.${target_experiment:2}"

source_nesting_experiment_path=$PWD
source ../REGION.src
#
echo ''
echo 'Using '${ESM_Scenario}' scenario'
echo 'The script will look for decadal files under '${Nesting_Files_PATH}'. Modify REGION.src if necessary'
echo 'The script assumes you are running in an experiment folder located in the Nesting Region (e.g. ESMa1.00/expt_01.0)'
echo ''
#
# Ensure Nesting_Files_PATH is set
if [ -z "$Nesting_Files_PATH" ]; then
    echo "Error: Nesting_Files_PATH is not set."
    exit 1
fi

sleep 5 # so the user can read whats above

genphynest=true
genbgcnest=true
#Generate the physical nesting files
if [ "$genphynest" = "true" ]; then

### Make sure all the files are prepared:
cd $Nesting_Files_PATH
for vari in thetao so uo vo zos; do
   num=`ls ${vari}_Omon_${ESM_Scenario}_g*_${year}*extrap* | wc -l`
   echo $year $vari $num
   if [ $num -ne 12 ]; then
     echo "Number of files incorrect for variable " $vari
     ${BINDIR}/Nesting_noresm/separate_and_extrapolate_files_year.sh $year $vari
   else
     echo "Extrapolation was already performed for "${vari}" before"
     echo ""
   fi
done

# Compute the nesting files
cd $source_nesting_experiment_path

${BINDIR}/Nesting_noresm/noresm_to_hycom.sh /cluster/work/users/$USER/$target_region/$experiment_full_name  \
   ${Nesting_Files_PATH}/thetao_Omon_${ESM_Scenario}_gr_${year}*_extrap.nc

# Compute the montgomery potential
#mkdir -p /cluster/work/users/${USER}/$target_region/nest
#mkdir -p /cluster/work/users/${USER}/$target_region/nest/$target_experiment
cd /cluster/work/users/${USER}/$target_region/nest/$target_experiment
#mkdir -p /cluster/work/users/${USER}/$target_region/nest/$target_experiment/Orig
for sday in 016 046 075 106 136 167 197 228 259 289 320 350 ; do
   mv archv.${year}_${sday}_00.* Orig/
done

cd /cluster/work/users/${USER}/$target_region/$experiment_full_name/
python ${BINDIR}/calc_montg1.py /cluster/work/users/${USER}/$target_region/nest/$target_experiment/Orig/archv.${year}_*_00.b \
       /cluster/work/users/${USER}/$target_region/$experiment_full_name/data/restart.2010_001_12_0000.b \
          /cluster/work/users/${USER}/$target_region/nest/$target_experiment/Montg/

cd $source_nesting_experiment_path #$SLURM_SUBMIT_DIR
fi

#Generate the biogeochemical nesting files
if [ "$genbgcnest" = "true" ]; then
### Make sure all the files are prepared:                                                                                  
cd ${Nesting_Files_PATH}
for vari in no3 po4 o2 si; do
   num=`ls ${vari}_Omon_${ESM_Scenario}_g*_${year}*extrap* | wc -l`
   echo $year $vari $num
   if [ $num -ne 12 ]; then
     echo "Number of files incorrect for variable " $vari
     srun -n1 -c2 --overlap ${BINDIR}/Nesting_noresm/separate_and_extrapolate_files_year.sh $year $vari
   fi
done
fi
exit
