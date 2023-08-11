#!/bin/csh -x
#
set echo
#
# --- Form a HYCOM restart file from a HYCOM archive file.
#
set in_year = $1
set in_day = $2
set in_day = `printf '%03d' $in_day`

set out_year = $1
set out_day = $2
set out_day = `printf '%03d' $out_day`

set dir = `pwd`
set parentdir = `dirname $dir`
set foldername = `basename $parentdir`

setenv RID = $foldername #TP5a0.06
source $pwd/EXPT.src

setenv R  /cluster/work/users/achoth//${RID}
setenv D /cluster/work/users/achoth/${RID}/expt_${X}/data
setenv O /cluster/work/users/achoth/${RID}/expt_${X}/data/NewRestart


#
#mkdir -p ${O}
cd       ${D}
#
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s ${R}/topo/depth_${RID}_{T}.a regional.depth.a
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s ${R}/topo/depth_${RID}_{T}.b regional.depth.b
endif
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s ${R}/topo/regional.grid.a .
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s ${R}/topo/regional.grid.b .
endif
#
# ---  input archive file
# ---  input restart file: (tempelate) can be any restart file from any date
# --- output restart file
#
###AO ####/bin/rm -f  ${O}/restart.2022_270_00_0000.a ${O}/restart.2022_270_00_0000.b 
#aprun -n 1 ~/HYCOM-tools/archive/src/archv2restart <<E-o-D
#~/HYCOM-tools/archive/src/archv2restart <<E-o-D
/cluster/home/achoth/AchHycom/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL/archive/src/archv2restart <<E-o-D
${D}/archv.${in_year}_${in_day}_00.a
${D}/restart.${in_year}_${in_day}_00_0000.a
${O}/restart.${out_year}_${out_day}_00_0000.a
 081     'iexpt ' = experiment number x10 (000=from archive file)
  3     'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
800     'idm   ' = longitudinal array size
760     'jdm   ' = latitudinal  array size
  0     'kapref' = thermobaric ref. state (-1=input,0=none,1,2,3=constant)
 50     'kdm   ' = number of layers
 25.0   'thbase' = reference density (sigma units)
300.0   'baclin' = baroclinic time step (seconds), int. divisor of 86400
  0     'rmontg' = pbavg correction from relax.montg file  (0=F,1=T)
E-o-D
