#! /bin/csh
#
set echo
set time = 1
set timestamp
C
C --- R is region name.
C --- T is topography number.
C --- K is number of layers.
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C
setenv R GOMd0.08
setenv T 02
setenv K 20
setenv E 306
setenv P hycom/${R}/expt_30.6/data
setenv D ~/$P
C
setenv OS `uname`
switch ($OS)
case 'Linux':
    which yod
    if (! $status) then
      setenv OS XT3
    endif
    which aprun
    if (! $status) then
C --- XT4 or XT5
      setenv OS XT4
    endif
    which dplace
    if (! $status) then
      setenv OS ICE
      source /usr/share/modules/init/csh
      module switch mpi/sgi_mpi-2.04 mpi/sgi_mpi-1.26
      module list
    endif
    which poe
    if (! $status) then
      setenv OS IDP
      module swap mpi mpi/intel/impi/4.1.0
      module load mkl
      module list
    endif
    breaksw
endsw
switch ($OS)
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    endif
    breaksw
case 'Linux':
    if (-e /export/a/$user) then
#              NRLSSC
      setenv S /export/a/${user}/$P
    else
#              Single Disk
      setenv S ~/$P/SCRATCH
    endif
    breaksw
case 'OSF1':
    if      (-e /work) then
#                  ERDC MSRC
      mkdir        /work/${user}
      chmod a+rx   /work/${user}
      setenv S     /work/${user}/$P
    else if (-e /workspace) then
#                  ASC MSRC
      mkdir        /workspace/${user}
      chmod a+rx   /workspace/${user}
      setenv S     /workspace/${user}/$P
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    endif
    breaksw
case 'IRIX64':
    if      (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    breaksw
case 'IDP':
case 'ICE':
case 'AIX':
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
#                  ARL MSRC
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
    endif
    breaksw
case 'unicos':
case 'unicosmk':
    if      (-e /work) then
#                  ERDC MSRC
      mkdir        /work/${user}
      chmod a+rx   /work/${user}
      mkdir        /work/${user}/$R
      chgrp $ACCT  /work/${user}/$R
      setenv S     /work/${user}/$P
    else
      mkdir        /tmp/${user}
      chmod a+rx   /tmp/${user}
      mkdir        /tmp/${user}/$R
      chgrp $ACCT  /tmp/${user}/$R
      setenv S     /tmp/${user}/$P
    endif
    breaksw
endsw
C
mkdir -p $S
cd       $S
mkdir -p SAVE
C
C --- For whole year runs.
C ---   Y01 initial model year of this run.
C ---   YXX is the last model year of this run, and the first of the next run.
C ---   A and B are identical, typically blank.
C --- For part year runs.
C ---   A is this part of the year, B is next part of the year.
C ---   Y01 is the start model year of this run.
C ---   YXX is the end   model year of this run, usually Y01.
C --- For a few hour/day run
C ---   A   is the start day and hour, of form "dDDDhHH".
C ---   B   is the end   day and hour, of form "dXXXhYY".
C ---   Y01 is the start model year of this run.
C ---   YXX is the end   model year of this run, usually Y01.
C --- Note that these variables are set by the .awk generating script.
C
setenv A "i"
setenv B "j"
setenv B2 "k"
setenv Y01 "099"
setenv YXX "099"
C
echo "Y01 =" $Y01 "YXX = " $YXX  "A =" ${A} "B =" ${B}
C
C --- local input files.
C
if (-e ${D}/../${E}y${Y01}${A}.limits) then
  /bin/cp ${D}/../${E}y${Y01}${A}.limits limits
else
#  use "LIMITI"  when starting a run after day zero.
#  use "LIMITS9" (say) for a 9-day run
  echo "LIMITS" | awk -f ${D}/../${E}.awk y01=${Y01} ab=${A} >! limits
endif
cat limits
C
if (-e ${D}/../ports.input_${Y01}${A}) then
  /bin/cp ${D}/../ports.input_${Y01}${A} ports.input
else
  /bin/cp ${D}/../ports.input ports.input
endif
C
if (-e ${D}/../tracer.input_${Y01}${A}) then
  /bin/cp ${D}/../tracer.input_${Y01}${A} tracer.input
else
  /bin/cp ${D}/../tracer.input tracer.input
endif
C
if (-e ${D}/../blkdat.input_${Y01}${A}) then
  /bin/cp ${D}/../blkdat.input_${Y01}${A} blkdat.input
else
  /bin/cp ${D}/../blkdat.input blkdat.input
endif
C
C --- check that iexpt from blkdat.input agrees with E from this script.
C
setenv EB `grep "'iexpt ' =" blk* | awk '{printf("%03d", $1)}'`
if ($EB != $E) then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
endif
#C
#C --- turn on detailed debugging.
#C
#touch PIPE_DEBUG
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
switch ($OS)
case 'SunOS':
case 'IDP':
case 'ICE':
#case 'Linux':
case 'OSF1':
case 'IRIX64':
case 'AIX':
case 'unicos':
case 'unicosmk':
    if      (-e ~wallcraf/bin/pget_navo) then
      setenv pget ~wallcraf/bin/pget_navo
      setenv pput ~wallcraf/bin/pput_navo
    else if (-e ~wallcraf/bin/pget) then
      setenv pget ~wallcraf/bin/pget
      setenv pput ~wallcraf/bin/pput
    else
      setenv pget /bin/cp
      setenv pput /bin/cp
    endif
    breaksw
default:
    setenv pget /bin/cp
    setenv pput /bin/cp
endsw
C
C --- input files from file server.
C
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
   ${pget} ${D}/../../topo/depth_${R}_${T}.a regional.depth.a &
endif
if (-z regional.depth.b) then
   ${pget} ${D}/../../topo/depth_${R}_${T}.b regional.depth.b &
endif
C
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
   ${pget} ${D}/../../topo/regional.grid.a regional.grid.a &
endif
if (-z regional.grid.b) then
   ${pget} ${D}/../../topo/regional.grid.b regional.grid.b &
endif
C
touch  forcing_cfsr_${Y01}.tar
if (-z forcing_cfsr_${Y01}.tar) then
   ${pget} ${D}/forcing_cfsr_${Y01}.tar . &
endif
C
#setenv SW ""
#setenv SW "kpar"
setenv SW "chl"
if ($SW != "") then
  touch  forcing.${SW}.a
  touch  forcing.${SW}.b
  if (-z forcing.${SW}.a) then
     ${pget} ${D}/../../force/seawifs/${SW}.a forcing.${SW}.a &
  endif
  if (-z forcing.${SW}.b) then
     ${pget} ${D}/../../force/seawifs/${SW}.b forcing.${SW}.b &
  endif
endif
C
touch relax.rmu.a relax.saln.a relax.temp.a relax.intf.a
touch relax.rmu.b relax.saln.b relax.temp.b relax.intf.b
if (-z relax.rmu.a) then
   ${pget} ${D}/../../relax/${E}/relax_rmu.a relax.rmu.a  &
endif
if (-z relax.rmu.b) then
   ${pget} ${D}/../../relax/${E}/relax_rmu.b relax.rmu.b  &
endif
if (-z relax.saln.a) then
   ${pget} ${D}/../../relax/${E}/relax_sal.a relax.saln.a &
endif
if (-z relax.saln.b) then
   ${pget} ${D}/../../relax/${E}/relax_sal.b relax.saln.b &
endif
if (-z relax.temp.a) then
   ${pget} ${D}/../../relax/${E}/relax_tem.a relax.temp.a &
endif
if (-z relax.temp.b) then
   ${pget} ${D}/../../relax/${E}/relax_tem.b relax.temp.b &
endif
if (-z relax.intf.a) then
   ${pget} ${D}/../../relax/${E}/relax_int.a relax.intf.a &
endif
if (-z relax.intf.b) then
   ${pget} ${D}/../../relax/${E}/relax_int.b relax.intf.b &
endif
C
touch relax.rmutr.a relax.trcr.a
touch relax.rmutr.b relax.trcr.b
if (-z relax.rmutr.a) then
   ${pget} ${D}/../../relax/${E}/tracer_rmu.a relax.rmutr.a  &
endif
if (-z relax.rmutr.b) then
   ${pget} ${D}/../../relax/${E}/tracer_rmu.b relax.rmutr.b  &
endif
if (-z relax.trcr.a) then
   ${pget} ${D}/../../relax/${E}/relax_trc.a relax.trcr.a &
endif
if (-z relax.trcr.b) then
   ${pget} ${D}/../../relax/${E}/relax_trc.b relax.trcr.b &
endif
C
touch iso.top.a
touch iso.top.b
if (-z iso.top.a) then
   ${pget} ${D}/../../relax/${E}/iso_top.a iso.top.a  &
endif
if (-z iso.top.b) then
   ${pget} ${D}/../../relax/${E}/iso_top.b iso.top.b  &
endif
C
touch thkdf4.a
touch thkdf4.b
if (-z thkdf4.a) then
   ${pget} ${D}/../../relax/${E}/thkdf4.a thkdf4.a  &
endif
if (-z thkdf4.b) then
   ${pget} ${D}/../../relax/${E}/thkdf4.b thkdf4.b  &
endif
C
C --- restart input
C
touch   restart_in.a restart_in.b restart_out.a restart_out.b restart_out1.a restart_out1.b
if (-z restart_in.b) then
  setenv RI "       0.00"
else
  setenv RI `head -2 restart_in.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
endif
if (-z restart_out.b) then
  setenv RO "       0.00"
else
  setenv RO `head -2 restart_out.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
endif
if (-z restart_out1.b) then
  setenv R1 "       0.00"
else
  setenv R1 `head -2 restart_out1.b | tail -1 | awk  '{printf("%11.2f\n", $5)}'`
endif
setenv LI `awk  '{printf("%11.2f\n", $1)}' limits`
C
if (`echo $LI | awk '{if ($1 <= 0.0) print 1; else print 0}'`) then
C --- no restart needed
  /bin/rm restart_in.a   restart_in.b
  /bin/rm restart_out.a  restart_out.b
  /bin/rm restart_out1.a restart_out1.b
else if (`echo $LI $RI | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C --- restart is already in restart_in
  /bin/rm restart_out.a  restart_out.b
  /bin/rm restart_out1.a restart_out1.b
else if (`echo $LI $RO | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C --- restart is in restart_out
  /bin/mv restart_out.a  restart_in.a
  /bin/mv restart_out.b  restart_in.b
  /bin/rm restart_out1.a restart_out1.b
else if (`echo $LI $R1 | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C ---   restart is in restart_out1
  /bin/mv restart_out1.a restart_in.a
  /bin/mv restart_out1.b restart_in.b
  /bin/rm restart_out.a  restart_out.b
else
C ---   get restart from permenant storage
  /bin/rm restart_in.a   restart_in.b
  /bin/rm restart_out.a  restart_out.b
  /bin/rm restart_out1.a restart_out1.b
  ${pget} ${D}/restart_${Y01}${A}.a restart_in.a &
  ${pget} ${D}/restart_${Y01}${A}.b restart_in.b &
endif
C
C --- Nesting input archive files for 3 months.
C
mkdir -p nest
C
if (-e ./nest) then
  setenv NA ${Y01}${A}
  setenv NB ${YXX}${B}
  setenv NB2 ${YXX}${B2}
  cd ./nest
  touch rmu.a rmu.b
  if (-z rmu.a) then
     ${pget} ${D}/../../relax/${E}/nest_rmu.a rmu.a &
  endif
  if (-z rmu.b) then
     ${pget} ${D}/../../relax/${E}/nest_rmu.b rmu.b &
  endif
  touch   arch.dummy.b
  /bin/rm arch*.[ab]
  touch archv_${NA}.tar
  if (-z archv_${NA}.tar) then
    ${pget} ${D}/nest/archv_${NA}.tar archv_${NA}.tar
  endif
  touch archv_${NB}.tar
  if (-z archv_${NB}.tar) then
    ${pget} ${D}/nest/archv_${NB}.tar archv_${NB}.tar &
  endif
  touch archv_${NB2}.tar
  if (-z archv_${NB2}.tar) then
    ${pget} ${D}/nest/archv_${NB2}.tar archv_${NB2}.tar &
  endif
  cd ..
endif
C
C --- let all file copies complete.
C
wait
C
C --- untar the forcing files
C
tar xvf forcing_cfsr_${Y01}.tar
C
C --- zero file length means no rivers.
C
if (-z forcing.rivers.a) then
   /bin/rm forcing.rivers.[ab]
endif
/bin/ls -laFq
C
if (-e ./nest) then
  ls -laFq nest
endif
