#!/bin/csh -x
#PBS -N 373o
#PBS -j oe
#PBS -o 373o.log
#PBS -W umask=027 
#PBS -l application=hycom
#PBS -l mppwidth=8
#PBS -l mppnppn=8
#PBS -l walltime=0:30:00
#PBS -A NRLSS03755C3J
#PBS -q challenge
#
set echo
set time = 1
set timestamp
date
C
C --- Create 3hrly interannual forcing from monthly climo forcing.
C ---
C --- prebuild this script similar to a model script
# --- awk -f 284.awk y01=098 ab=a 284O.com > 284o098a.com
C
C --- Preamble, script keys on O/S name.
C
setenv OS `uname`
switch ($OS)
case 'SunOS':
C   assumes /usr/5bin is before /bin and /usr/bin in PATH.
    breaksw
case 'Linux':
    which aprun
    if (! $status) then
C --- XT4 or XT5
      setenv OS XT4
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
case 'Linux':
case 'XT4':
case 'OSF1':
case 'AIX':
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
endsw
C
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C --- N is data-set name, e.g. ec10m-reanal
C
setenv E 306
setenv P hycom/GOMd0.08/expt_30.6/data
setenv D ~/$P
setenv N RS_SST-mon
C
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
    else
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    endif
    breaksw
case 'Linux':
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
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    endif
    mkdir -p         $S
    lfs setstripe -d $S
    lfs setstripe    $S 1048576 -1 8
    breaksw
case 'OSF1':
    mkdir        ~/scratch
    chmod a+rx   ~/scratch
    setenv S     ~/scratch/$P
    breaksw
case 'IRIX64':
    mkdir        /workspace/${user}
    chmod a+rx   /workspace/${user}
    setenv S     /workspace/${user}/$P
    breaksw
case 'AIX':
    if (-e /gpfs/work) then
#                  ERDC MSRC
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
    endif
    breaksw
case 'unicos':
case 'unicosmk':
    mkdir        /tmp/${user}
    chmod a+rx   /tmp/${user}
    mkdir        /tmp/${user}/ATLd0.08
    chgrp $ACCT  /tmp/${user}/ATLd0.08
    setenv S     /tmp/${user}/$P
    breaksw
endsw
C
mkdir -p $S/ssto
cd       $S/ssto
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
setenv Y01 "011"
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
C --- use "LIMITI" when starting a run after day zero.
C
setenv TS `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $1}'`
setenv TM `echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} | awk '{print $2}'`
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
foreach f ( seatmp )
  touch ${f}.a
  if (-z ${f}.a) then
    ${pget} ${D}/../../force/${N}/${f}.a ${f}.a &
  endif
  touch ${f}.b
  if (-z ${f}.b) then
    ${pget} ${D}/../../force/${N}/${f}.b ${f}.b &
  endif
end
C
C --- executable
C
if (-e    ~wallcraf/hycom/ALL/force/src/time_interp) then
  /bin/cp ~wallcraf/hycom/ALL/force/src/time_interp . &
else
  /bin/cp ~wallcraf/hycom/ALLcnl/force/src/time_interp . &
endif
wait
chmod ug+rx time_interp
ls -laFq
C
C --- NAMELIST input.
C
touch   fort.05i
/bin/rm fort.05i
cat <<E-o-D  > fort.05i
 &AFTIME
  FINC   = 0.041666666667,  !1hrly
  FSTART = ${TS},  !TIME OF HEAT FLUX START              (DAYS)
  WSTART = ${TS},  !TIME OF WIND START                   (DAYS)
  TSTART = ${TS},  !TIME OF START OF CURRENT INTEGRATION (DAYS)
  TMAX   = ${TM},  !TIME OF END   OF CURRENT INTEGRATION (DAYS)
 /
 &AFFLAG
  INTERP =   3,  !cubic
 /
E-o-D
C
C --- run the wind extraction
C
touch      fort.10 fort.10a
/bin/rm -f fort.10 fort.10a
C
foreach f ( seatmp )
  setenv FOR010A fort.10a
C
C --- Input.
C
  setenv FOR020  ${f}.b
  setenv FOR020A ${f}.a
C
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
case 'AIX':
    /bin/rm -f core
    touch core
    ./time_interp < fort.05i
    breaksw
case 'XT4':
    /bin/rm -f core
    touch core
    aprun -n 1 ./time_interp < fort.05i
    breaksw
case 'IRIX64':
    /bin/rm -f core
    touch core
    assign -V
    ./time_interp < fort.05i
    assign -V
    assign -R
    breaksw
case 'unicosmk':
    /bin/rm -f core
    touch core
    assign -V
    ./time_interp < fort.05i
    if (! -z core)  debugview time_interp core
    assign -V
    assign -R
    breaksw
case 'unicos':
    /bin/rm -f core
    touch core
    assign -V
    ./time_interp < fort.05i
    if (! -z core)  debug -s time_interp core
    assign -V
    assign -R
    breaksw
endsw
C
C --- Output.
C
  /bin/mv fort.10  ${f}_${Y01}${A}.b
  /bin/mv fort.10a ${f}_${Y01}${A}.a
C
  if (-e ./SAVE) then
    ln ${f}_${Y01}${A}.a ./SAVE/${f}_${Y01}${A}.a
    ln ${f}_${Y01}${A}.b ./SAVE/${f}_${Y01}${A}.b
  endif
end
C
C  --- END OF JUST IN TIME SEATMP EXTRACTION SCRIPT.
C
date
#
#/usr/lpp/LoadL/full/bin/llq -w $LOADL_STEP_ID
