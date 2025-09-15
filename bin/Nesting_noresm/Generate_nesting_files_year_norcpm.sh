#!/bin/bash -l

#SBATCH --account=nn9481k
#SBATCH --job-name=gennest
#SBATCH -t 0:30:00
#SBATCH --qos=devel
#SBATCH --nodes=1   # number of nodes
#SBATCH  --mail-type=ALL
#SBATCH --mail-user=

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

year=$1
# How to call (example)
# Generate_nesting_files_year_norcpm.sh 1980

# These variables have to be set prior to running the script
lname=NorESM2-MM_historical_r1i1p1f1
#lname=NorESM2-MM_ssp585_r1i1p1f1
genphynest=false
genbgcnest=true
ESM_datadir=/cluster/work/users/uname/NORESM/               # Name of directory containing ESM files
Nest_dir=/cluster/work/users/uname/ESMa1.00/expt_01.0/      # Name of directory for making nesting files
Model_dir=/cluster/work/users/uname/TZ4a0.10/expt_01.1/     # Name of directory for nested model
E=011

#Generate the physical nesting files
if [ "$genphynest" = "true" ]; then
### Make sure all the files are prepared:
cd ESM_datadir
for vari in thetao so uo vo zos; do
   num=`ls ${vari}_Omon_${lname}_g*_${year}*extrap* | wc -l`
   echo $year $vari $num
   if [ $num -ne 12 ]; then
     echo "Number of files incorrect for variable " $vari
     ./separate_and_extrapolate_files_year.sh $year $vari
   fi
done

# Compute the nesting files
cd $Nest_dir
../bin/Nesting_noresm/noresm_to_hycom.sh $Model_dir  \
$ESM_dir/thetao_Omon_${lname}_gr_${year}*_extrap.nc

# Compute the montgomery potential
cd $Model_dir
cd ../nest/$E
for sday in 016 046 075 106 136 167 197 228 259 289 320 350 ; do
   mv archv.${year}_${sday}_00.* Orig/
done

cd $Model_dir
python ./bin/calc_montg1.py ../nest/011/Orig/archv.${year}_*_00.b data/restart.2005_079_12_0000.b ../nest/$E/Montg/

cd $SLURM_SUBMIT_DIR
fi

#Generate the biogeochemical nesting files
if [ "$genbgcnest" = "true" ]; then
### Make sure all the files are prepared:                                                                                  
cd $ESM_dir
for vari in no3 po4 o2 si; do
   num=`ls ${vari}_Omon_${lname}_g*_${year}*extrap* | wc -l`
   echo $year $vari $num
   if [ $num -ne 12 ]; then
     echo "Number of files incorrect for variable " $vari
     srun -n1 -c2 --overlap ./separate_and_extrapolate_files_year.sh $year $vari
   fi
done
fi
exit
