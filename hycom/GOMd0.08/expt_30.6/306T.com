#!/bin/csh -x
#PBS -N XXX
#PBS -j oe
#PBS -o XXX.log
#PBS -W umask=027
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l place=scatter:excl
#PBS -l walltime=12:00:00
#PBS -l application=home-grown
#PBS -A NRLSS03755018
#PBS -q standard
#
set echo
set time = 1
set timestamp
C
C --- Create model interpolated surtmp for HYCOM.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284T.com > 284t098a.com
C
C --- Preamble, script keys on O/S name.
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
C   assumes /usr/5bin is before /bin and /usr/bin in PATH.
    breaksw
case 'Linux':
    which poe
    if (! $status) then
      setenv OS IDP
    endif
    breaksw
case 'OSF1':
    breaksw
case 'IRIX64':
    breaksw
case 'AIX':
    breaksw
case 'unicosmk':
    setenv ACCT `groups | awk '{print $1}'`
    breaksw
case 'unicos':
    setenv ACCT `newacct -l | awk '{print $4}'`
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
switch ($OS)
case 'SunOS':
case 'OSF1':
case 'AIX':
case 'IDP':
case 'unicos':
case 'unicosmk':
    if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget cp
      setenv pput cp
    endif
    breaksw
case 'IRIX64':
    setenv pget cp
    setenv pput cp
    breaksw
case 'Linux':
    setenv pget msfget
    setenv pput msfget
    breaksw
default:
    setenv pget cp
    setenv pput cp
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- X is data-set executable abbreviation, e.g. 1125_ec
C --- N is data-set name, e.g. ecmwf-reanal_ds111.6
C --- W is permanent native surtmp directory
C
setenv E 306
setenv P hycom/GOMd0.08/expt_30.6/data
setenv D ~/$P
setenv X kp_nc
setenv N ncep_cfsr/netcdf
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
      setenv W     /net/hermes/scrb/metzger/temp_ieee/$N
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    endif
    breaksw
case 'Linux':
case 'XT3':
    if (-e /work) then
#                  NRLSSC
      setenv S /work/${user}/$P
      setenv W /erdc2/metzger/force/$N
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv W ~/temp_ieee/$N
    endif
    breaksw
case 'OSF1':
#                 ERDC MSRC
    mkdir        /work/${user}
    chmod a+rx   /work/${user}
    setenv S     /work/${user}/$P
    setenv W     /u/home/metzger/temp_ieee/$N
    breaksw
case 'IRIX64':
    mkdir        /workspace/${user}
    chmod a+rx   /workspace/${user}
    setenv S     /workspace/${user}/$P
    setenv W1    /msas031/metzger/temp_ieee/$N1
    setenv W2    /msas031/metzger/temp_ieee/$N2
    breaksw
case 'AIX':
case 'IDP':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv W     /u/home/metzger/temp_ieee/$N
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/force/$N
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv W     /archive/navy/metzger/temp_ieee/$N
    endif
    breaksw
case 'unicos':
case 'unicosmk':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/GLBa0.72
    chgrp $ACCT  /tmp/${user}/GLBa0.72
    setenv S     /tmp/${user}/$P
    setenv W     /u/home/metzger/temp_ieee/$N
    breaksw
endsw
C
mkdir -p $S/ssta
cd       $S/ssta
C
C --- For whole year runs.
C ---   ymx number of years per model run.
C ---   Y01 initial model year of this run.
C ---   YXX is the last model year of this run, and the first of the next run.
C ---   A and B are identical, typically blank.
C --- For part year runs.
C ---   A is this part of the year, B is next part of the year.
C ---   Y01 initial model year of this run.
C ---   YXX year at end of this part year run.
C ---   ymx is 1.
C --- Note that these variables and the .awk generating script must
C ---  be consistant to get a good run script.
C
C --- For surtmp, only Y01 and A are used.
C
C --- One year spin-up run.
C
@ ymx =  1
C
setenv A "a"
setenv B "b"
setenv Y01 "096"
C
switch ("${B}")
case "${A}":
    setenv YXX `echo $Y01 $ymx | awk '{printf("%03d", $1+$2)}'`
    breaksw
case "a":
    setenv YXX `echo $Y01 | awk '{printf("%03d", $1+1)}'`
    breaksw
default:
    setenv YXX $Y01
endsw
C
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
C
C --- time limits.
C
if (-e ${D}/../${E}y${Y01}${A}.limits) then
  setenv TS `sed -e "s/-/ /g" ${D}/../${E}y${Y01}${A}.limits | awk '{print $1}'`
  setenv TM `cat              ${D}/../${E}y${Y01}${A}.limits | awk '{print $2}'`
else
# use "LIMITI" when starting a run after day zero.
# use "LIMITS9" (say) for a 9-day run.
  setenv TS `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $1}'`
  setenv TM `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $2}'`
endif
C
echo "TS =" $TS "TM =" $TM
C
setenv YC1 `echo $Y01 | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF071 cfsr-sea_${YC1}_01hr_surtmp.nc
setenv YCX `echo $YXX | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF072 cfsr-sea_${YCX}_01hr_surtmp.nc
C
C --- input files from file server.
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.b) then
  ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
if (-z regional.grid.a) then
  ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
C     
touch  $CDF071
if (-z $CDF071) then
  ${pget} ${W}/$CDF071 . &
endif
C
if ($CDF071 != $CDF072) then
  touch  $CDF072
  if (-z $CDF072) then
    ${pget} ${W}/$CDF072 . &
  endif
endif
C
C --- executable
C
/bin/cp ~wallcraf/hycom/ALL/force/src/${X} . &
wait
chmod ug+rx ${X}
ls -laFq
C
C --- NAMELIST input.
C
touch   fort.05in
/bin/rm fort.05in
cat <<E-o-D  > fort.05in
 &AFTITL
  CTITLE = '1234567890123456789012345678901234567890',
  CTITLE = 'CFSR-sea, 1hrly, deg C',
  CNAME  = 'surtmp',
 &END
 &AFTIME
  FSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
  PARMIN = -999.0,  !disable parmin
  PARMAX =  999.0,  !disable parmax
  PAROFF = -273.16, !K to degC
  TMPICE =   -1.79, !sea ice marker
 &END
 &AFFLAG
  IFFILE =   5,  !3:monthly; 5:actual day;
  INTERP =   0,  !0:piecewise-linear; 1:cubic-spline;
  INTMSK =   0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 &END
E-o-D
switch ($OS)
case 'unicos':
case 'unicos-lanl':
case 'unicosmk':
case 'unicos-t3d':
case 'sn6705':
case 'AIX':
case 'IDP':
C
C --- Fortran 90 NAMELIST delimiter
C
  /bin/rm -f fort.05i
  sed -e 's/&END/\//' -e 's/&end/\//' -e '/^#/d' < fort.05in > fort.05i
  breaksw
default:
C
C --- Fortran 77 NAMELIST delimiter
C
  /bin/rm -f fort.05i
  cp fort.05in fort.05i
  breaksw
endsw
C
C --- run the surtmp interpolation.
C
touch fort.10 fort.10a
/bin/rm -f fort.10 fort.10a
C
setenv FOR010A fort.10a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
case 'IDP':
    /bin/rm -f core
    touch core
    ./${X} < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./${X} < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicosmk':
    /bin/rm -f core
    touch core
    assign -V
    ./${X} < fort.05i
    if (! -z core)  debugview ${X} core
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./${X} < fort.05i
    if (! -z core)  debug -s ${X} core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output.
C
/bin/mv fort.10  surtmp_${Y01}${A}.b
/bin/mv fort.10a surtmp_${Y01}${A}.a
C
C --- CICE
C
if (-e ../cice) then
  cd ../cice
  setenv IDM 500
  setenv JDM 382
# incoming Qlw from netQlw and Tsur, just calculate surtmp4 here
  /bin/rm surtmp_${Y01}${A}.a surtmp2_${Y01}${A}.a surtmp4_${Y01}${A}.a
  hycom_expr ../ssta/surtmp_${Y01}${A}.a ONE ${IDM} ${JDM} 1.0 273.15 surtmp_${Y01}${A}.a >! surtmp_${Y01}${A}.b
  /bin/cp     surtmp_${Y01}${A}.a  surtmp_${Y01}${A}.A
  hycom_expr  surtmp_${Y01}${A}.a  surtmp_${Y01}${A}.A  ${IDM} ${JDM} 0.0 0.0 surtmp2_${Y01}${A}.a >! surtmp2_${Y01}${A}.b
  /bin/cp    surtmp2_${Y01}${A}.a surtmp2_${Y01}${A}.A
  hycom_expr surtmp2_${Y01}${A}.a surtmp2_${Y01}${A}.A ${IDM} ${JDM} 0.0 0.0 surtmp4_${Y01}${A}.a >! surtmp4_${Y01}${A}.b
endif
C
C  --- END OF JUST IN TIME SURTMP GENERATION SCRIPT.
C
