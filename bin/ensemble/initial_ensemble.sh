#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# To extract the initial restart file for an ensemble
# It requires to prefine the time period a sepific time period (one or some months)
# and used the relevant script of restart_hycomcice.sh
# 
# input or predefine parameters:
#
#  Example:
#
#   <> start_year end_year start_month end_month Oy Om Od <Sdir>
#   <> start_year end_year start_month end_month Oy Om Od <Sdir> <Odir>
#   <> start_year end_year start_month end_month Oy Om Od <Sdir> <Sdir_ice> <Odir>
#  
#  suggest:
# 
#  initial_ensemble.sh 1993 2012 7 9 2019 7 20 /cluster/work/users/xiejp/work_temp_data/data_081 /cluster/work/users/xiejp/TOPAZ/TP5a0.06/expt_04.0  
#  3Jan 2019 by Jiping Xie  
# Updated 31 March 2023 to add a requirement for the time interval(Dj0>15 days)
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [ $# -lt 8 -o $# -gt 10 ]; then
  echo " Error input parameters "
  echo "  initial_ensemble.sh start_year end_year start_month end_month Oy Om Od <Sdir>  "
  echo " .or. ~ start_year end_year start_month end_month Oy Om Od <Sdir> <Odir> "
  echo " .or. ~ start_year end_year start_month end_month Oy Om Od <Sdir> <Sdir_ice> <Odir> "
  exit 0
elif [ $# -eq 8 ]; then
  Srdir=$8  
  Ordir= "./" 
  Srdir_c=${Srdir}/cice
  Ordir_c=${Ordir}/cice
elif [ $# -eq 9 ]; then
  Srdir=$8  
  Ordir=$9 
  Srdir_c=${Srdir}/cice
  Ordir_c=${Ordir}/cice
elif [ $# -eq 10 ]; then
  Srdir=$8  
  Ordir=$10
  Srdir_c=$9
  Ordir_c=${Ordir}/cice
fi

echo ${Ordir}

Y1=$1
Y2=$2
mm1=$3
mm2=$4
Oyy=$5
Omm=$6
Odd=$7

echo "$Y1 $Y2 $mm1 $mm2 $Oyy $Omm $Odd"
echo "${Srdir_c}" 
echo "${Ordir_c}"
echo ""

# Step one: ensure the directoreis are existing
if [ ! -r ${Srdir} ]; then
   echo "missing the restart files under ${Srdir}"
   exit 0
fi
if [ ! -r ${Srdir_c} ]; then
   echo "missing the iced files under ${Srdir_c}"
   exit 0
fi
[ ! -r ${Ordir} ] && mkdir ${Ordir}
[ ! -r ${Ordir_c} ] && mkdir ${Ordir_c}

# Step two: checking the reasonability for the input parameters
if [ $Y1 -gt $Y2 -o $mm1 -lt 1 -o $mm2 -lt 1 -o $Omm -lt 1 -o $Omm -gt 12 -o $Odd -lt 1 -o $Odd -gt 31 ]; then
   echo "Wrong in the input 7 parameters!"
   exit 0
fi
if [ $mm2 -lt $mm1 ]; then
   echo "Wrong in the input parameters: $mm1 $mm2"
   exit 0
else 
   if [ $mm2 -lt 12 ]; then
     ((mm2 = mm2 + 1))
   fi
fi

kk=0     # before the label of the first member
#kk=22      # before the label of the first member
(( Ojdy = 1 + `datetojul ${Oyy} ${Omm} ${Odd} ${Oyy} 1 1 | tail -4c` ))

Dj0=15
for yy in `seq ${Y1} ${Y2}`; do
  Jdy1=$(datetojul ${yy} ${mm1} 1 ${yy} 1 1 | tail -4c) 
  Jdy2=$(datetojul ${yy} ${mm2} 1 ${yy} 1 1 | tail -4c) 
  echo ${Jdy1} ${Jdy2}
  if [ ${Jdy1} -lt 1 -a ${Jdy2} -lt 1 ]; then
    echo "Wrong in the input mm1(mm2) to fix the month window!"
    exit 0
  fi
  J1=0
  for jj in `seq ${Jdy1} ${Jdy2}`; do
     Sdate=$(jultodate ${jj} ${yy} 1 1)
     (( Sjdy = 1 + `datetojul ${Sdate:0:4} ${Sdate:4:2} ${Sdate:6:2} ${Sdate:0:4} 1 1 | tail -4c` ))
     Sfice="iced.${Sdate:0:4}-`echo 0${Sdate:4:2}|tail -3c`-`echo 0${Sdate:6:2}|tail -3c`-00000.nc"
     Sfhycom="restart.${Sdate:0:4}_`echo 00${Sjdy}|tail -4c`_00_0000"
     if [ -s ${Srdir}/${Sfhycom}.a -a -s ${Srdir}/${Sfhycom}.b -a ${Srdir_c}/${Sfice} ]; then
        (( Dj = jj - J1 ))
        if [ ${Dj} -gt ${Dj0} ]; then
           restart_hycomcice.sh ${Sdate:0:4} ${Sdate:4:2} ${Sdate:6:2} ${Oyy} ${Omm} ${Odd} ${Srdir} ${Ordir} 
           Fhycom="${Ordir}/restart.${Oyy}_`echo 00${Ojdy}|tail -4c`_00_0000"
           Fcice="${Ordir_c}/iced.${Oyy}-`echo 0${Omm}|tail -3c`-`echo 0${Odd}|tail -3c`-00000"
  
           #echo "${Fhycom}.a -a -s ${Fhycom}.b -a ${Fcice}.nc "
           if [ -s ${Fhycom}.a -a -s ${Fhycom}.b -a ${Fcice}.nc ]; then
              let kk=kk+1
              Msur="mem`echo 00${kk}|tail -4c`"
              mv ${Fhycom}.a ${Fhycom}_${Msur}.a
              mv ${Fhycom}.b ${Fhycom}_${Msur}.b
              mv ${Fcice}.nc ${Fcice}_${Msur}.nc
              (( J1 = jj ))
           fi
        fi
     fi
  done
done


