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
C --- Create model interpolated winds for HYCOM. Use bi-cubic interpolation.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284W.com > 284w098a.com
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
case 'Linux':
    setenv pget msfget
    setenv pput msfget
    breaksw
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- N is data-set name, e.g. ec10m-reanal
C --- W is permanent native pcip directory
C
setenv E 306
setenv P hycom/GOMd0.08/expt_30.6/data
setenv D ~/$P
setenv N ncep_cfsr/netcdf
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
      setenv W     /net/hermes/scrb/metzger/wind_ieee/$N
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      setenv W     /u/b/metzger/wind_ieee/$N
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
      setenv W ~/wind_ieee/$N
    endif
    breaksw
case 'OSF1':
    mkdir        ~/scratch
    chmod a+rx   ~/scratch
    setenv S     ~/scratch/$P
    setenv W     /u/b/metzger/wind_ieee/$N
    breaksw
case 'IRIX64':
    mkdir        /workspace/${user}
    chmod a+rx   /workspace/${user}
    setenv S     /workspace/${user}/$P
    setenv W1    /msas031/metzger/wind_ieee/$N1
    setenv W2    /msas031/metzger/wind_ieee/$N2
    breaksw
case 'AIX':
case 'IDP':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv W     /u/b/metzger/wind_ieee/$N
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
      setenv W     /archive/navy/metzger/wind_ieee/$N
    endif
    breaksw
case 'unicos':
case 'unicosmk':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/ATLd0.08
    chgrp $ACCT  /tmp/${user}/ATLd0.08
    setenv S     /tmp/${user}/$P
    setenv W     /u/b/metzger/wind_ieee/$N
    breaksw
endsw
C
mkdir -p $S/wind
cd       $S/wind
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
C --- For winds, only Y01 and A are used.
C
C --- One year spin-up run.
C
@ ymx =  1
C
setenv A "a"
setenv B "b"
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
setenv YC1 `echo $Y01 | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF071 cfsr-sec2_${YC1}_01hr_strblk.nc
setenv YCX `echo $YXX | awk '{printf("%04d\n",$1+1900)}'`
setenv CDF072 cfsr-sec2_${YCX}_01hr_strblk.nc
C
C --- input files from file server.
C
touch      fort.44 fort.44a fort.45 fort.45a
if (-z fort.44) then
  ${pget} $D/../../force/offset/tauewd_zero.b fort.44  &
endif
if (-z fort.44a) then
  ${pget} $D/../../force/offset/tauewd_zero.a fort.44a &
endif
if (-z fort.45) then
  ${pget} $D/../../force/offset/taunwd_zero.b fort.45  &
endif
if (-z fort.45a) then
  ${pget} $D/../../force/offset/taunwd_zero.a fort.45a &
endif
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
/bin/cp ~wallcraf/hycom/ALL/force/src/wi_nc . &
wait
chmod ug+rx wi_nc
ls -laFq
C
C --- NAMELIST input.
C
touch   fort.05in
/bin/rm fort.05in
cat <<E-o-D  > fort.05in
 &WWTITL
  CTITLE = '123456789012345678901234567890123456789012345678901234567890',
  CTITLE = 'CFSR-sec2, 1hrly, MKS',
 /
 &WWTIME
  SPDMIN =   0.0,  !minimum allowed wind speed
  WSCALE =   1.0,  !scale factor to mks
  WSTART = ${TS},
  TSTART = ${TS},
  TMAX   = ${TM},
 /
 &WWFLAG
  IGRID  = 2,  !0=p; 1=u&v; 2=p
  ISPEED = 0,  !0:none; 1:const; 2:kara; 3:coare
  INTERP = 3,  !0:bilinear; 1:cubic spline; 2:piecewise bessel; 3:piecewise bi-cubic
  INTMSK = 0,  !0:no mask; 1:land/sea=0/1; 2:l/s=1/0;
  IFILL  = 3,  !0,1:tx&ty; 2,3:magnitude; 1,3:smooth; (intmsk>0 only)
  IOFILE = 0,  !0:single offset; 1:multiple offsets; 2:multi-off no .b check
  IWFILE = 4,  !1:ann/mon; 2:multi-file; 4:actual wind day
 /
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
C --- run the wind interpolation.
C
touch fort.10 fort.10a
/bin/rm -f fort.1[012] fort.1[012]a
C
setenv FOR010A fort.10a
setenv FOR011A fort.11a
setenv FOR012A fort.12a
C
setenv FOR044A fort.44a
setenv FOR045A fort.45a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
case 'IDP':
    /bin/rm -f core
    touch core
    ./wi_nc < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./wi_nc < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicosmk':
    /bin/rm -f core
    touch core
    assign -V
    ./wi_nc < fort.05i
    if (! -z core)  debugview wi_nc core
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./wi_nc < fort.05i
    if (! -z core)  debug -s wi_nc core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output.
C
/bin/mv fort.10  tauewd_${Y01}${A}.b
/bin/mv fort.10a tauewd_${Y01}${A}.a
/bin/mv fort.11  taunwd_${Y01}${A}.b
/bin/mv fort.11a taunwd_${Y01}${A}.a
if (-e fort.12) then
  /bin/mv fort.12  wndspd_${Y01}${A}.b
  /bin/mv fort.12a wndspd_${Y01}${A}.a
endif
C
C  --- END OF JUST IN TIME WIND GENERATION SCRIPT.
C
