#! /bin/bash -l
#
#  Make sure I use the correct shell.
#
###SBATCH -A nn2993k 
#SBATCH -A nn9481k 
#SBATCH --partition=bigmem
#SBATCH --exclusive
#SBATCH --ntasks=16 --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#
#SBATCH -J ERA5prep
#SBATCH -N 1
#
#SBATCH -t 24:00:00
#SBATCH -e /cluster/work/users/xiejp/ERA5_hourly/Epre%J.err
#SBATCH -o /cluster/work/users/xiejp/ERA5_hourly/Epre%J.out
#

# Before convert the hourly NC into 6-hour files
# Modules of NCO and CDO are default loaded by the system.
# Otherwise the loading should be added here. 

prgdrt='/cluster/home/xiejp/ERA5down_2025/'
prg=${prgdrt}convert/prep_hourly24.sh

[ ! -s ${prgdrt}env.drt ] && { echo "warning: check ${prgdrt}env.drt to setup the directory"; exit 0; }
source env.drt

if [ ! -r ${Rundir} ]; then
   echo "warning: check the directory ${Rundir} to start later."
   exit 0
else
   cd ${Rundir}
   echo $(pwd)
if
Outdir=${Rundir}/hourly
[ ! -r ${Outdir} ] && mkdir hourly

if [ $# -lt 1 ]; then
   yr=2025
else
   yr=$1
fi

echo "$(ls ./download/hourly.*_${yr}.nc)"
i0=0
for ii in $(ls ./download/hourly.*_${yr}.nc); do
    echo ${ii}
    Fend=./hourly/hourly.${ii#*hourly.}
    if [ ! -s ${Fend} -a -s $ii ]; then
       ${prg} ${ii}  &
       (( i0 = i0 + 1 ))
       if [ $i0 -eq 8 ]; then
          wait 
          sleep 2s
          i0=0
       fi
    fi
done 
wait
sleep 2s

exit $?

