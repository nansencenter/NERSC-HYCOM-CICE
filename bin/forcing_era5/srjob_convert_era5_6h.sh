#! /bin/bash -l
#
#  Make sure I use the correct shell.
#
#SBATCH -A nn9481k 
#SBATCH --partition=preproc
#SBATCH --ntasks=5 --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#
#SBATCH -J ERA5C
#
#SBATCH -t 24:00:00
#SBATCH -e /cluster/work/users/xiejp/ERA5hourly/E0%J.err
#SBATCH -o /cluster/work/users/xiejp/ERA5hourly/E0%J.out
#

# double check the diretory of "prgdrt" is consitent
#  with the definition of "Inidir" in evn.drt:
prgdrt=/cluster/home/xiejp/Developing/ERA5
prg=${prgdrt}/convert/do_6hourly.sh
#
[ ! -s ${prgdrt}/env.drt ] && { echo "warning: check ${prgdrt}env.drt to setup the directory"; exit 0; }
source ${prgdrt}/env.drt

if [ ! -r ${Rundir} ]; then
   echo "warning: check the directory ${Rundir} to start later."
   exit 0
fi

# keep the fixed order:
Vars="10U 10V 2D 2T BLH MSL TCC CI SKT SSR STR SSRD STRD RO TP"

V0=0


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

cd ${Rundir}
echo $(pwd)
echo "${Y2}"
i0=0
for iyy in `seq $Y1 $Y2`; do
    if [ $V0 -gt 0 ]; then
       ii=$V0
    else
       ii=0
    fi
    for ivar in ${Vars}; do
       (( ii +=1 ))	   
       Fout=6h.${ivar}_${iyy}.nc
       Fini=./hourly/hourly.${ivar}_${iyy}.nc
       echo ${Fini}
       if [ ! -s ./6h/${Fout} -a -s ${Fini} ]; then
          echo "${Fout}"
          echo "${prg} ${iyy} ${ii} "
          ##srun -N1 -n1 ${prg} ${iyy} ${ii} 
          ${prg} ${iyy} ${ii} & 
	  (( i0++ ))
	  if [ $i0 -gt 4 ]; then
             wait
	     sleep 3s
	     i0=0
	  fi
       fi
    done
done
wait
sleep 2s
exit $?
