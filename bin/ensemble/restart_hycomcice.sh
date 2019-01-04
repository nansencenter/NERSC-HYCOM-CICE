#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# used to prepare the restart files of HYCOM_CICE according to the defined date
# input required:  previous date and source diretory (like $Sdir and $Sdir/cice if not defined)
# outputs:         target date which determine the filename and output directory (similar as above)
#
#  Note: the time steps in HYCOM and in CICE model should be double checked
#     T_hycom=600  and T_cice=900 (s) are defoult values in the script
# 
#  and the refered dates are 1901/1/1 and 1958/1/1
#    H_yy=1901  H_mm=1 H_dd=1
#    I_yy=1958  I_mm=1 I_dd=1
#
#  Source directory: initial restart
#  Sdir='/cluster/work/users/xiejp/work_2018/TP2_EXP/020_hycom';
#  Sdir_ice:  ${Sdir}/cice
#
#  Ouput directory: target restart
#  Odir='/cluster/work/users/xiejp/TOPAZ/TP2a0.10/expt_02.0/data';
#  Odir_ice:  ${Sdir}/cice
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Example:
#
#  restart_hycomcice.sh Py Pm Pd Oy Om Od 
#  .or.  restart_hycomcice.sh Py Pm Pd Oy Om Od <Sdir> <Odir>
#  .or.  restart_hycomcice.sh Py Pm Pd Oy Om Od <Sdir> <Sdir_ice> <Odir>
#  
# 2Jan 2019 by Jiping Xie  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

if [ $# -lt 6 -o $# -gt 9 ]; then
   echo " Error input parameters "
   echo "  restart_hycomcice.sh Py Pm Pd Oy Om Od "
   echo " .or.  restart_hycomcice.sh Py Pm Pd Oy Om Od <Sdir> <Odir>  "
   echo " .or.  restart_hycomcice.sh Py Pm Pd Oy Om Od <Sdir> <Sdir_ice> <Odir>  "
   exit 0
elif [ $# -eq 6 ]; then
   Sdir='./'
   Sdir_ice='./cice'
   Odir='./'
   Odir_ice='./cice'
elif [ $# -eq 7 ]; then
   Sdir=$7
   Odir='./'
   Sdir_ice=${Sdir}'/cice'
   Odir_ice=${Odir}'/cice'
elif [ $# -eq 8 ]; then
   Sdir=$7
   Odir=$8
   Sdir_ice=${Sdir}'/cice'
   Odir_ice=${Odir}'/cice'
elif [ $# -eq 9 ]; then
   Sdir=$7
   Sdir_ice=$8
   Odir=$9
   Odir_ice=${Odir}'/cice'
fi

# input parameters
Pyy=$1
Pmm=$2
Pdd=$3

Oyy=$4
Omm=$5
Odd=$6

# parameters predifined
T_hycom=600         #  hycom setting
H_yy=1901 
H_mm=1 
H_dd=1

T_cice=900          # cice setting
I_yy=1958  
I_mm=1 
I_dd=1

# first check for the input parameters (model day)
(( Pjdy = 1 + `datetojul ${Pyy} ${Pmm} ${Pdd} ${Pyy} 1 1 | tail -4c` ))
if [ ${Pjdy} -gt 366 -o ${Pjdy} -lt 0 ];then
   echo "Wrong input for initial date: " ${Pyy}-${Pmm}-${Pdd}
   exit 0
fi
(( Ojdy = 1 + `datetojul ${Oyy} ${Omm} ${Odd} ${Oyy} 1 1 | tail -4c` ))
if [ ${Ojdy} -gt 366 -o ${Ojdy} -lt 0 ];then
   echo "Wrong input for output date: " ${Oyy}-${Omm}-${Odd}
   exit 0
fi

#Sdir='/cluster/work/users/xiejp/work_2018/TP2_EXP/020_hycom';
#Sdir_ice='/cluster/work/users/xiejp/work_2018/TP2_EXP/020_cice';
# detecting the source directory
if [ ! -r ${Sdi_ice} ]; then
  Sdir_ice=${Sdir}/cice  
  if [ ! -r ${Sdir_ice} ]; then
     echo "Cannot found the source diretory: " ${Sdir}
     exit 0
  fi
fi

# check the availability about the restart file to be linked further
Pfice="iced.${Pyy}-`echo 0${Pmm}|tail -3c`-`echo 0${Pdd}|tail -3c`-00000.nc"
Pfhycom="restart.${Pyy}_`echo 00${Pjdy}|tail -4c`_00_0000"
if [ ! -s ${Sdir}/${Pfhycom}'.a' -o ! -s ${Sdir}/${Pfhycom}'.b' ]; then
   echo "missing the hycom restart file:  ${Sdir}/${Pfhycom}.[ab]"
   exit 0
fi
if [ ! -s ${Sdir_ice}/${Pfice} ]; then
   echo "missing the cice restart file:  ${Sdir_ice}/${Pfice}"
   exit 0
fi

Ofice="iced.${Oyy}-`echo 0${Omm}|tail -3c`-`echo 0${Odd}|tail -3c`-00000.nc"
Ofhycom="restart.${Oyy}_`echo 00${Ojdy}|tail -4c`_00_0000"
echo 'Input: ' ${Pfice} ' , ' ${Pfhycom}
echo 'Output: ' ${Ofice} ' and ' ${Ofhycom}

# Step one: copy the files
#Odir='/cluster/work/users/xiejp/TOPAZ/TP2a0.10/expt_02.0/data';
if [ ! -r ${Odir_ice} ]; then
   Odir_ice=${Odir}/cice;
   echo 'setup the out diretory ' ${Odir}
   mkdir ${Odir}
   mkdir ${Odir_ice}
fi
cp ${Sdir}/${Pfhycom}.a ${Odir}/${Ofhycom}.a
cp ${Sdir}/${Pfhycom}.b ${Odir}/${Ofhycom}.b
cp ${Sdir_ice}/${Pfice} ${Odir_ice}/${Ofice}

# Step two:  modify the information of nstep and dtime in .b file
echo ""
echo "modify the nstep and dtime in .b file ... "
(( Hstep = 86400 / ${T_hycom} ))
(( var_time = 1 + `datetojul ${Oyy} ${Omm} ${Odd} ${H_yy} ${H_mm} ${H_dd}` ))
(( var_nstep = ${Hstep} * ${var_time} ))
echo "HYCOM steps in one day: ${Hstep} so there are :  ${var_time} ${var_nstep} "
cd ${Odir}
Nl=`wc -l ${Ofhycom}.b | awk '{print $1}'`
sed -n "1,1p" ${Ofhycom}.b > ${Ofhycom}.tmp    # first line
# second line
Hline=$(sed -n '2,2p' ${Ofhycom}.b)
Bspa='   '
Fspa='     '
echo "${Hline%% = *} = ${Fspa}${var_nstep}${Bspa}${var_time}.0000000000${Hline#*.0000000000}" >> ${Ofhycom}.tmp
sed -n "3,${Nl}p" ${Ofhycom}.b >> ${Ofhycom}.tmp
[ -s ${Ofhycom}.tmp ] && mv ${Ofhycom}.tmp ${Ofhycom}.b

# Step three:  modify the information of nstep and dtime in cice-nc file
echo ""
echo "modify the nstep and dtime in cice-nc file ... "

(( ice_st = 86400 / ${T_cice} ))
(( ice_dy = `datetojul ${Oyy} ${Omm} ${Odd} ${I_yy} ${I_mm} ${I_dd}` ))
(( istep = ${ice_st} * ${ice_dy} ))
(( itime = ${istep} * ${T_cice} ))

echo "CIEE steps in one day: ${ice_st} so there are ${istep} ${itime}"
module load NCO/4.6.6-intel-2017a
cd ${Odir_ice}
ncatted -h -a istep1,global,o,i,${istep} ${Ofice} out1.nc
ncatted -O -h -a time,global,o,i,${itime} out1.nc
[ -s out1.nc ] && mv out1.nc ${Ofice} 

echo ""
echo "check the existence under $(pwd), and then try to run a simple case"
echo "normal endding"


