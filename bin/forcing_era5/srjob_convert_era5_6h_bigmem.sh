#! /bin/bash -l
#
#  Make sure I use the correct shell.
#
#SBATCH -A nn9481k 
#SBATCH --partition=bigmem
#SBATCH --exclusive
#SBATCH --ntasks=4 --cpus-per-task=4
#SBATCH --mem-per-cpu=15G
#
#
#SBATCH -J ERA5conv
#SBATCH -N 1
#
#SBATCH -t 48:00:00
#SBATCH -e /cluster/work/users/xiejp/ERA5_hourly/Econ-%J.err
#SBATCH -o /cluster/work/users/xiejp/ERA5_hourly/Econ-%J.out
#

prgdrt=/cluster/home/xiejp/ERA5down_2025
prg=${prgdrt}/convert/do_6hourly.sh
#
[ ! -s ${prgdrt}/env.drt ] && { echo "warning: check ${prgdrt}env.drt to setup the directory"; exit 0; }
source ${prgdrt}/env.drt

if [ ! -r ${Rundir} ]; then
   echo "warning: check the directory ${Rundir} to start later."
   exit 0
else
   cd ${Rundir}
   echo $(pwd)
fi

# keep the fixed order:
Vars="10U 10V 2D 2T BLH MSL TCC CI SKT SSR STR SSRD STRD RO TP"


if [ $# -eq 2 ]; then
   Y1=$1
   Y2=$2
elif [ $# -eq 1 ]; then
   Y1=$1
   Y2=$1
else
   Y1=2009
   Y2=2009
fi
echo "${Y2}"
for iyy in `seq $Y1 $Y2`; do
    ii=0
    for ivar in ${Vars}; do
       (( ii = ii + 1 ))	   
       Fout=6h.${ivar}_${iyy}.nc
       echo "${Fout}"
       Fini=./hourly/hourly.${ivar}_${iyy}.nc
       if [ ! -s ./6h/${Fout} -a -s ${Fini} ]; then
          echo "${Fout}"
          echo "${prg} ${iyy} ${ii} "
          #srun -N1 -n1 ${prg} ${iyy} ${ii} 
          ${prg} ${iyy} ${ii} & 
       fi
    done
    wait
    sleep 10s
done
exit $?
