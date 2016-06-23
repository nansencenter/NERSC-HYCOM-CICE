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
date +"START  %c"
C
C --- Create model interpolated lwdflx for HYCOM.
C --- From NWP NRL .D files.
C
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284L.com > 284l098a.com
C
C --- Preamble, script keys on O/S name.
C
setenv SRC	~wallcraf/hycom/ALL/bin
setenv APRUN	""
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
case 'Linux':
    which yod
    if (! $status) then
      setenv OS XT3
    endif
    which aprun
    if (! $status) then
C --- XT4, XT5 or XC30
      setenv OS XT4
      setenv SRC	~wallcraf/hycom/ALLcnl/bin
      setenv APRUN	"aprun -n 1"
      sleep 300
    endif
    which dplace
    if (! $status) then
      setenv OS ICE
      source /usr/share/modules/init/csh
    endif
    which poe
    if (! $status) then
      setenv OS IDP
    endif
    breaksw
endsw
#
switch ($OS)
case 'SunOS':
C   assumes /usr/5bin is before /bin and /usr/bin in PATH.
    breaksw
case 'Linux':
case 'ICE':
case 'IDP':
    breaksw
case 'XT4':
    breaksw
case 'OSF1':
    breaksw
case 'IRIX64':
    breaksw
case 'AIX':
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
case 'Linux':
case 'ICE':
case 'IDP':
case 'XT4':
case 'OSF1':
case 'AIX':
case 'unicos':
    if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget cp
      setenv pput cp
    endif
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
C --- W is permanent native lwdflx directory
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
      setenv W     /net/hermes/scrb/metzger/force/$N
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/force/$N
    endif
    breaksw
case 'Linux':
    if (-e /external/fast) then
#                  NRLSSC
      setenv S /external/fast/${user}/$P
      setenv W ~/force/$N
    else if (-e /work) then
#                  ERDC
      setenv S /work/${user}/$P
      setenv W ~/force/$N
      mkdir -p         $S
      lfs setstripe -d $S
      lfs setstripe    $S -s 1048576 -i -1 -c 8
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/force/$N
    else
#                  Single Disk
      setenv S ~/$P/SCRATCH
      setenv W ~/force/$N
    endif
    breaksw
case 'ICE':
case 'IDP':
case 'XT3':
case 'XT4':
    if      (-e /work) then
#                  ERDC MSRC
      mkdir        /work/${user}
      chmod a+rx   /work/${user}
      setenv S     /work/${user}/$P
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/home/metzger/force/$N
    else if (-e /p/work1) then
#                  NAVO MSRC
      mkdir        /p/work1/${user}
      chmod a+rx   /p/work1/${user}
      setenv S     /p/work1/${user}/$P
      setenv W     /u/home/metzger/force/$N
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    endif
    mkdir -p         $S
    lfs setstripe -d $S
    lfs setstripe    $S -s 1048576 -i -1 -c 8
    breaksw
case 'OSF1':
#                 ERDC MSRC
    mkdir        /work/${user}
    chmod a+rx   /work/${user}
    setenv S     /work/${user}/$P
    setenv W     /u/home/metzger/force/$N
    breaksw
case 'IRIX64':
    mkdir        /workspace/${user}
    chmod a+rx   /workspace/${user}
    setenv S     /workspace/${user}/$P
    setenv W     /msas031/metzger/force/$N
    breaksw
case 'AIX':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv W     /u/home/metzger/force/$N
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
      setenv W     /archive/navy/metzger/force/$N
    endif
    breaksw
case 'unicos':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/GOMd0.08
    chgrp $ACCT  /tmp/${user}/GOMd0.08
    setenv S     /tmp/${user}/$P
    setenv W     /u/home/metzger/force/$N
    breaksw
endsw
C
mkdir -p $S/lrad
cd       $S/lrad
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
C --- For lwdflx, only Y01 and A are used.
C
C --- One year spin-up run.
C
@ ymx =  1
C
setenv A "j"
setenv B "k"
setenv Y01 "099"
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
touch *.nc
setenv YC1 `echo $Y01 | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF071 ./cfsr-sea_${YC1}_01hr_dlwsfc.nc
setenv YCX `echo $YXX | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF072 ./cfsr-sea_${YCX}_01hr_dlwsfc.nc
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
/bin/cp ${SRC}/../force/src/${X} . &
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
  CTITLE = 'CFSR, 1-hrly, ocn-only, W/m^2',
  CNAME  = 'dlwflx',
 &END
 &AFTIME
  FSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
  PARMIN = -9999.0,  !disable parmin
  PARMAX =  9999.0,  !disable parmax
  PAROFF =     0.0,  !no offset
 &END
 &AFFLAG
  IFFILE =   5,  !3:monthly; 5:actual day;
  INTERP =   0,  !0:bilinear; 1:cubic spline; 2:piecewise bessel; 3:piecewise bi-cubic;
  INTMSK =   0,  !0:no mask; 1:land/sea=0/1; 2:land/sea=1/0;
 &END
E-o-D
switch ($OS)
case 'unicos':
case 'IDP':
case 'AIX':
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
C --- run the lwdflx interpolation.
C
date +"lwdflx %c"
touch fort.10 fort.10a
/bin/rm -f fort.10 fort.10a
C
setenv FOR010A fort.10a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'ICE':
case 'IDP':
case 'OSF1':
case 'AIX':
    /bin/rm -f core
    touch core
    ${APRUN} ./${X} < fort.05i
    breaksw
case 'XT4':
    /bin/rm -f core
    touch core
    ${APRUN} ./${X} < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ${APRUN} ./${X} < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ${APRUN} ./${X} < fort.05i
    if (! -z core)  debug -s ${X} core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output.
C
/bin/mv fort.10  lwdflx_${Y01}${A}.b
/bin/mv fort.10a lwdflx_${Y01}${A}.a
C
if (-e ./SAVE) then
  ln lwdflx_${Y01}${A}.a ./SAVE/lwdflx_${Y01}${A}.a
  ln lwdflx_${Y01}${A}.b ./SAVE/lwdflx_${Y01}${A}.b
endif
C
C  --- END OF JUST IN TIME LDWFLX GENERATION SCRIPT.
C
date +"END    %c"
