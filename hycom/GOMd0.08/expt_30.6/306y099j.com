#! /bin/csh
#
# --- check that the C comment command is available.
#
C >& /dev/null
if (! $status) then
  if (-e ${home}/hycom/ALL/bin/C) then
    set path = ( ${path} ${home}/hycom/ALL/bin )
  else
    echo "Please put the command hycom/ALL/bin/C in your path"
  endif
endif
#
set echo
set time = 1
set timestamp
C
C --- Experiment GOMd0.08 - 30.X
C --- 20 layer HYCOM on 0.08 degree Gulf of Mexico region
C
C --- 30.0 - CFSR forcing, 10m winds and pressure.
C --- 30.1 - CFSR forcing, 10m winds, no pressure input.
C --- 30.2 - twin of 30.0 with 2.2.97
C --- 30.3 - twin of 30.2 with 2.2.97 and flxflg=6.
C --- 30.4 - twin of 30.0 with 2.2.97W
C --- 30.5 - twin of 30.1 with 2.2.97W
C --- 30.6 - twin of 30.3 with 2.2.98 and downward lw and sw
C
C --- Preamble, script keys on O/S name.
C
C --- Set parallel configuration, see ../README/README.expt_parallel.
C --- NOMP = number of OpenMP threads, 0 for no OpenMP, 1 for inactive OpenMP
C --- NMPI = number of MPI    tasks,   0 for no MPI
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
      module swap mpi mpi/intel/impi/4.1.3
      module load mkl
      module list
    endif
    breaksw
endsw
switch ($OS)
case 'SunOS':
case 'Linux':
case 'OSF1':
    setenv NOMP 0
    setenv NMPI 16
    breaksw
case 'IDP':
case 'ICE':
case 'XT3':
case 'XT4':
    setenv NOMP 0
    setenv NMPI 32
    breaksw
case 'IRIX64':
case 'AIX':
    setenv NOMP 0
    setenv NMPI 32
    breaksw
case 'unicos':
    setenv ACCT `newacct -l | awk '{print $4}'`
    setenv NOMP 0
    setenv NMPI 0
    breaksw
default:
    echo 'Unknown Operating System: ' $OS
    exit (1)
endsw
C
C --- modify NOMP and NMPI based on batch limits
C
if ( $?LSB_MCPU_HOSTS ) then
# LSF batch system
# if ( $?LSB_INITIAL_NUM_PROCESSORS) then
#   setenv NCPU $LSB_INITIAL_NUM_PROCESSORS
# else
#   setenv NCPU `echo $LSB_MCPU_HOSTS | awk '{print $2+$4+$6+$8+$10+$12}'`
# endif
# if      ($NMPI == 0) then
#   setenv NOMP $NCPU
# else if ($NOMP == 0) then
#   setenv NMPI $NCPU
# else
#   setenv NMPI `echo $NCPU $NOMP | awk '{print int($1/$2)}'`
# endif
else if ( $?GRD_TOTAL_MPI_TASKS ) then
# GRD batch system
  if      ($NMPI == 0) then
    echo "error - NMPI=0, but running in a MPI batch queue"
    exit
  else
    setenv NMPI $GRD_TOTAL_MPI_TASKS
  endif
else if ( $?NSLOTS ) then
# codine or GRD batch system
  if      ($NMPI == 0) then
    setenv NOMP $NSLOTS
  else if ($NOMP == 0) then
    setenv NMPI $NSLOTS
  else
    setenv NMPI `echo $NSLOTS $NOMP | awk '{print int($1/$2)}'`
  endif
endif
echo "NOMP is " $NOMP " and NMPI is " $NMPI
C
C --- R is region name.
C --- V is source code version number.
C --- T is topography number.
C --- K is number of layers.
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C
setenv R GOMd0.08
setenv V 2.2.98
setenv T 02
setenv K relo
setenv E 306
setenv P hycom/${R}/expt_30.6/data
setenv D ~/$P
C
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
    else if (-e /scr) then
#                  NAVO MSRC
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
    else
#              Single Disk
      setenv S ~/$P/SCRATCH
    endif
    breaksw
case 'IDP':
case 'ICE':
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
    else
#                  Single Disk
      setenv S     ~/$P/SCRATCH
    endif
    mkdir -p         $S
    lfs setstripe -d $S
    lfs setstripe    $S 1048576 -1 8
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
case 'AIX':
    if      (-e /gpfs/work) then
#                  ERDC MSRC, under PBS
      mkdir        /gpfs/work/${user}
      chmod a+rx   /gpfs/work/${user}
      setenv S     /gpfs/work/${user}/$P
      setenv POE  pbspoe
    else if (-e /scr) then
#                  NAVO MSRC, under LoadLeveler or LSF
      mkdir        /scr/${user}
      chmod a+rx   /scr/${user}
      setenv S     /scr/${user}/$P
      if ($?LSB_JOBINDEX) then
        setenv POE mpirun.lsf
      else
        setenv POE poe
      endif
    else
#                  ARL MSRC, under GRD
      mkdir        /usr/var/tmp/${user}
      chmod a+rx   /usr/var/tmp/${user}
      setenv S     /usr/var/tmp/${user}/$P
      setenv POE  grd_poe
    endif
    breaksw
case 'unicos':
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
setenv A "j"
setenv B "k"
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
if (-e    ${D}/../ports_z.input_${Y01}${A}) then
  /bin/cp ${D}/../ports_z.input_${Y01}${A} ports_z.input
  /bin/cp ${D}/../ports_u.input_${Y01}${A} ports_u.input
  /bin/cp ${D}/../ports_v.input_${Y01}${A} ports_v.input
  /bin/cp ${D}/../ports_a.input_${Y01}${A} ports_a.input
else
  /bin/cp ${D}/../ports_z.input ports_z.input
  /bin/cp ${D}/../ports_u.input ports_u.input
  /bin/cp ${D}/../ports_v.input ports_v.input
  /bin/cp ${D}/../ports_a.input ports_a.input
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
if (-e      ${D}/../archs.input_${Y01}${A}) then
  /bin/cp   ${D}/../archs.input_${Y01}${A} archs.input
else if (-e ${D}/../archs.input) then
  /bin/cp   ${D}/../archs.input archs.input
endif
C
if (-e      ${D}/../profile.input_${Y01}${A}) then
  /bin/cp   ${D}/../profile.input_${Y01}${A} profile.input
else if (-e ${D}/../profile.input) then
  /bin/cp   ${D}/../profile.input profile.input
else
  touch profile.input
endif
if (! -z profile.input) then
  if (-e       ./ARCHP) then
    /bin/mv -f ./ARCHP ./ARCHP_$$
  endif
  mkdir ./ARCHP
endif
C
if (-e ./cice) then
  if (-e ${D}/../ice_in_${Y01}${A}) then
    /bin/cp ${D}/../ice_in_${Y01}${A} ice_in
  else
    /bin/cp ${D}/../ice_in ice_in
  endif
endif
if ($NMPI != 0) then
# setenv NPATCH `echo $NMPI | awk '{printf("%04d", $1)}'`
  setenv NPATCH `echo $NMPI | awk '{printf("%03d", $1)}'`
  /bin/rm -f patch.input
  /bin/cp ${D}/../../topo/partit/depth_${R}_${T}.${NPATCH}  patch.input
# /bin/cp ${D}/../../topo/partit/depth_${R}_${T}.${NPATCH}u patch.input
C
  /bin/rm -f archt.input
  if (-e ${D}/../archt.input_${Y01}${A}) then
    /bin/cp ${D}/../archt.input_${Y01}${A} archt.input
  else if (-e ${D}/../archt.input) then
    /bin/cp ${D}/../archt.input archt.input
  else
    touch archt.input
  endif
  if (! -z archt.input) then
    if (-e       ./ARCHT) then
      /bin/mv -f ./ARCHT ./ARCHT_$$
    endif
    mkdir ./ARCHT
    switch ($OS)
    case 'ICE':
    case 'XT3':
    case 'XT4':
      lfs setstripe ./ARCHT 1048576 -1 8
      breaksw
    endsw
    cd    ./ARCHT
    /bin/cp ../archt.input  .
    /bin/cp ../patch.input  .
    /bin/ln ../regional.*.? .
    ${home}/hycom/ALL/topo/src/topo_subset < archt.input
    cat patch.subreg | xargs mkdir
    /bin/rm -f regional.*.?
    cd ..
  endif
endif
C
C --- check that iexpt from blkdat.input agrees with E from this script.
C
setenv EB `grep "'iexpt ' =" blk* | awk '{printf("%03d", $1)}'`
if ($EB != $E) then
  cd $D/..
  /bin/mv -f LIST LIST_BADRUN
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
#case 'Linux':
case 'OSF1':
case 'AIX':
case 'unicos':
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
case 'IRIX64':
    setenv pget /bin/cp
    setenv pput /bin/cp
    breaksw
case 'IDP':
case 'ICE':
case 'XT3':
case 'XT4':
#   setenv pget /bin/cp
    setenv pget ~wallcraf/bin/pget
    setenv pput ~wallcraf/bin/pput
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
if (-e ./cice) then
C
C --- CICE non-synoptic files.
C
  touch regional.cice.r cice.prec_lanl_12.r cice.rhoa_ncar85-88_12.r
  if (-z regional.cice.r) then
     ${pget} ${D}/../../topo/regional.cice.r regional.cice.r &
  endif
  if (-z cice.prec_lanl_12.r) then
     ${pget} ${D}/../../force/cice/cice.prec_lanl_12.r cice.prec_lanl_12.r &
  endif
  if (-z cice.rhoa_ncar85-88_12.r) then
     ${pget} ${D}/../../force/cice/cice.rhoa_ncar85-88_12.r cice.rhoa_ncar85-88_12.r &
  endif
endif
C
if (! -e ./wind) then
C
C --- Climatological atmospheric forcing.
C
  setenv FN era40-sec_1978-2002_mn6hra
  touch forcing.tauewd.a forcing.taunwd.a forcing.wndspd.a
  touch forcing.radflx.a forcing.shwflx.a forcing.vapmix.a forcing.precip.a
  touch forcing.airtmp.a forcing.seatmp.a forcing.surtmp.a forcing.mslprs.a
  touch forcing.tauewd.b forcing.taunwd.b forcing.wndspd.b
  touch forcing.radflx.b forcing.shwflx.b forcing.vapmix.b forcing.precip.b
  touch forcing.airtmp.b forcing.seatmp.b forcing.surtmp.b forcing.mslprs.b
  if (-z forcing.tauewd.a) then
     ${pget} ${D}/../../force/${FN}/tauewd.a      forcing.tauewd.a &
  endif
  if (-z forcing.tauewd.b) then
     ${pget} ${D}/../../force/${FN}/tauewd.b      forcing.tauewd.b &
  endif
  if (-z forcing.taunwd.a) then
     ${pget} ${D}/../../force/${FN}/taunwd.a      forcing.taunwd.a &
  endif
  if (-z forcing.taunwd.b) then
     ${pget} ${D}/../../force/${FN}/taunwd.b      forcing.taunwd.b &
  endif
  if (-z forcing.wndspd.a) then
     ${pget} ${D}/../../force/${FN}/wndspd.a      forcing.wndspd.a &
  endif
  if (-z forcing.wndspd.b) then
     ${pget} ${D}/../../force/${FN}/wndspd.b      forcing.wndspd.b &
  endif
  if (-z forcing.vapmix.a) then
     ${pget} ${D}/../../force/${FN}/vapmix.a      forcing.vapmix.a &
  endif
  if (-z forcing.vapmix.b) then
     ${pget} ${D}/../../force/${FN}/vapmix.b      forcing.vapmix.b &
  endif
  setenv AO ""
# setenv AO "_037c"
  if (-z forcing.airtmp.a) then
     ${pget} ${D}/../../force/${FN}/airtmp${AO}.a forcing.airtmp.a &
  endif
  if (-z forcing.airtmp.b) then
     ${pget} ${D}/../../force/${FN}/airtmp${AO}.b forcing.airtmp.b &
  endif
# setenv PO "+2mmweek"
  setenv PO ""
# setenv PO "_zero"
  if (-z forcing.precip.a) then
     ${pget} ${D}/../../force/${FN}/precip${PO}.a forcing.precip.a &
  endif
  if (-z forcing.precip.b) then
     ${pget} ${D}/../../force/${FN}/precip${PO}.b forcing.precip.b &
  endif
  setenv FO ""
# setenv FO "-s14w"
  if (-z forcing.radflx.a) then
     ${pget} ${D}/../../force/${FN}/radflx${FO}.a forcing.radflx.a &
  endif
  if (-z forcing.radflx.b) then
     ${pget} ${D}/../../force/${FN}/radflx${FO}.b forcing.radflx.b &
  endif
  if (-z forcing.shwflx.a) then
     ${pget} ${D}/../../force/${FN}/shwflx${FO}.a forcing.shwflx.a &
  endif
  if (-z forcing.shwflx.b) then
     ${pget} ${D}/../../force/${FN}/shwflx${FO}.b forcing.shwflx.b &
  endif
  if (-z forcing.surtmp.a) then
     ${pget} ${D}/../../force/${FN}/surtmp.a      forcing.surtmp.a &
  endif
  if (-z forcing.surtmp.b) then
     ${pget} ${D}/../../force/${FN}/surtmp.b      forcing.surtmp.b &
  endif
  setenv FS $FN
# setenv FS PF_SST-mn6hr
  if (-z forcing.seatmp.a) then
     ${pget} ${D}/../../force/${FS}/seatmp.a      forcing.seatmp.a &
  endif
  if (-z forcing.seatmp.b) then
     ${pget} ${D}/../../force/${FS}/seatmp.b      forcing.seatmp.b &
  endif
endif
C
C --- time-invarent wind stress offset
C
setenv OFS ""
#setenv OFS "_era40-nogaps"
if ($OFS != "") then
  touch  forcing.ofstrs.a
  touch  forcing.ofstrs.b
  if (-z forcing.ofstrs.a) then
     ${pget} ${D}/../../force/offset/ofstrs${OFS}.a forcing.ofstrs.a &
  endif
  if (-z forcing.ofstrs.b) then
     ${pget} ${D}/../../force/offset/ofstrs${OFS}.b forcing.ofstrs.b &
  endif
endif
C
C --- time-invarent heat flux offset
C
setenv OFF ""
#setenv OFF "_${E}"
if ($OFF != "") then
  touch  forcing.offlux.a
  touch  forcing.offlux.b
  if (-z forcing.offlux.a) then
     ${pget} ${D}/../../force/offset/offlux${OFF}.a forcing.offlux.a &
  endif
  if (-z forcing.offlux.b) then
     ${pget} ${D}/../../force/offset/offlux${OFF}.b forcing.offlux.b &
  endif 
endif 
C
touch  forcing.rivers.a
touch  forcing.rivers.b
if (-z forcing.rivers.a) then
   ${pget} ${D}/../../force/rivers/rivers_${T}.a forcing.rivers.a &
endif
if (-z forcing.rivers.b) then
   ${pget} ${D}/../../force/rivers/rivers_${T}.b forcing.rivers.b &
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
#C
#touch tbaric.a
#touch tbaric.b
#if (-z tbaric.a) then
#   ${pget} ${D}/../../relax/${E}/tbaric.a tbaric.a  &
#endif
#if (-z tbaric.b) then
#   ${pget} ${D}/../../relax/${E}/tbaric.b tbaric.b  &
#endif
C
setenv XS ""
#setenv XS "010"
if ($XS != "") then
  touch  relax.ssh.a
  if (-z relax.ssh.a) then
     ${pget} ${D}/../../relax/SSH/relax_ssh_${XS}.a relax.ssh.a &
  endif
  touch  relax.ssh.b
  if (-z relax.ssh.b) then
     ${pget} ${D}/../../relax/SSH/relax_ssh_${XS}.b relax.ssh.b &
  endif
endif
C
setenv XR ""
#setenv XR "glb"
if ($XR != "") then
  touch  relax.sssrmx.a
  if (-z relax.sssrmx.a) then
     ${pget} ${D}/../../relax/SSSRMX/sssrmx_${XR}.a relax.sssrmx.a &
  endif
  touch  relax.sssrmx.b
  if (-z relax.sssrmx.b) then
     ${pget} ${D}/../../relax/SSSRMX/sssrmx_${XR}.b relax.sssrmx.b &
  endif
endif
#C
#touch diwlat.a
#touch diwlat.b
#if (-z diwlat.a) then
#   ${pget} ${D}/../../relax/${E}/diwlat.a diwlat.a  &
#endif
#if (-z diwlat.b) then
#   ${pget} ${D}/../../relax/${E}/diwlat.b diwlat.b  &
#endif
#C
#touch iso.sigma.a
#touch iso.sigma.b
#if (-z iso.sigma.a) then
#   ${pget} ${D}/../../relax/${E}/iso_sigma.a iso.sigma.a  &
#endif
#if (-z iso.sigma.b) then
#   ${pget} ${D}/../../relax/${E}/iso_sigma.b iso.sigma.b  &
#endif
#C
#touch iso.top.a
#touch iso.top.b
#if (-z iso.top.a) then
#   ${pget} ${D}/../../relax/${E}/iso_top.a iso.top.a  &
#endif
#if (-z iso.top.b) then
#   ${pget} ${D}/../../relax/${E}/iso_top.b iso.top.b  &
#endif
#C
#touch thkdf4.a
#touch thkdf4.b
#if (-z thkdf4.a) then
#   ${pget} ${D}/../../relax/${E}/thkdf4.a thkdf4.a  &
#endif
#if (-z thkdf4.b) then
#   ${pget} ${D}/../../relax/${E}/thkdf4.b thkdf4.b  &
#endif
#C
#touch veldf2.a
#touch veldf2.b
#if (-z veldf2.a) then
#   ${pget} ${D}/../../relax/${E}/veldf2.a veldf2.a  &
#endif
#if (-z veldf2.b) then
#   ${pget} ${D}/../../relax/${E}/veldf2.b veldf2.b  &
#endif
#C
#touch veldf4.a
#touch veldf4.b
#if (-z veldf4.a) then
#   ${pget} ${D}/../../relax/${E}/veldf4.a veldf4.a  &
#endif
#if (-z veldf4.b) then
#   ${pget} ${D}/../../relax/${E}/veldf4.b veldf4.b  &
#endif
C
setenv TZ ""
#setenv TZ "_10mm"
C
if ($TZ != "") then
  touch  cb.a
  if (-z cb.a) then
     ${pget} ${D}/../../relax/DRAG/cb_${T}${TZ}.a cb.a  &
  endif
  touch  cb.b
  if (-z cb.b) then
     ${pget} ${D}/../../relax/DRAG/cb_${T}${TZ}.b cb.b  &
  endif
endif
#C
#setenv TT ""
##setenv TT ".lim5"
#C
#touch tidal.tensor.a
#touch tidal.tensor.b
#if (-z tidal.tensor.a) then
#   ${pget} ${D}/../../relax/DRAG/tidal.Ntensor.${T}${TT}.a tidal.tensor.a  &
#endif
#if (-z tidal.tensor.b) then
#   ${pget} ${D}/../../relax/DRAG/tidal.Ntensor.${T}${TT}.b tidal.tensor.b  &
#endif
#C
#touch tidal.rh.a
#touch tidal.rh.b
#if (-z tidal.rh.a) then
#   ${pget} ${D}/../../relax/${E}/tidal.rh.a tidal.rh.a  &
#endif
#if (-z tidal.rh.b) then
#   ${pget} ${D}/../../relax/${E}/tidal.rh.b tidal.rh.b  &
#endif
C
setenv TS ""
#setenv TS $T
if ($TS != "") then
  touch  tidal.sal.a
  if (-z tidal.sal.a) then
     ${pget} ${D}/../../relax/SAL/tpxo8aM2_salQtide_${TS}.a tidal.sal.a &
  endif
  touch  tidal.sal.b
  if (-z tidal.sal.b) then
     ${pget} ${D}/../../relax/SAL/tpxo8aM2_salQtide_${TS}.b tidal.sal.b &
  endif
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
  /bin/rm -f restart_in.a   restart_in.b
  /bin/rm -f restart_out.a  restart_out.b
  /bin/rm -f restart_out1.a restart_out1.b
else if (`echo $LI $RI | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C --- restart is already in restart_in
  /bin/rm -f restart_out.a  restart_out.b
  /bin/rm -f restart_out1.a restart_out1.b
else if (`echo $LI $RO | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C --- restart is in restart_out
  /bin/mv -f restart_out.a  restart_in.a
  /bin/mv -f restart_out.b  restart_in.b
  /bin/rm -f restart_out1.a restart_out1.b
else if (`echo $LI $R1 | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C ---   restart is in restart_out1
  /bin/mv -f restart_out1.a restart_in.a
  /bin/mv -f restart_out1.b restart_in.b
  /bin/rm -f restart_out.a  restart_out.b
else
C ---   get restart from permenant storage
  /bin/rm -f restart_in.a   restart_in.b
  /bin/rm -f restart_out.a  restart_out.b
  /bin/rm -f restart_out1.a restart_out1.b
  ${pget} ${D}/restart_${Y01}${A}.a restart_in.a &
  ${pget} ${D}/restart_${Y01}${A}.b restart_in.b &
endif
if (-e ./cice) then
C
C --- CICE restart input
C
touch   cice.restart_in cice.restart_out
if (-z cice.restart_in) then
  setenv RI "       0.00"
else
  setenv RI `cice_stat cice.restart_in  | awk  '{printf("%11.2f\n", $4)}'`
endif
if (-z cice.restart_out) then
  setenv RO "       0.00"
else
  setenv RO `cice_stat cice.restart_out | awk  '{printf("%11.2f\n", $4)}'`
endif
setenv LI `awk  '{printf("%11.2f\n", $1)}' limits`
C
if (`echo $LI $RI | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C --- restart is already in cice.restart_in
  /bin/rm -f cice.restart_out
else if (`echo $LI $RO | awk '{if ($1-0.1 < $2 && $1+0.1 > $2) print 1; else print 0}'`) then
C --- restart is in cice.restart_out
  /bin/mv cice.restart_out  cice.restart_in
else
C ---   get restart from permenant storage
  /bin/rm -f cice.restart_in
  /bin/rm -f cice.restart_out
  ${pget} ${D}/cice.restart_${Y01}${A} cice.restart_in &
endif
echo "cice.restart_in" >! cice.restart_file
endif
C
C --- model executable
C
if      ($NMPI == 0 && $NOMP == 0) then
  setenv TYPE one
else if ($NMPI == 0) then
  setenv TYPE omp
else if ($NOMP == 0) then
  if ( ! $?TYPE ) then
    setenv TYPE mpi
  endif
else
  setenv TYPE ompi
endif
if (-e ./cice) then
  setenv TYPE cice
  setenv HEXE hycom_cice
else
  setenv HEXE hycom
endif
/bin/cp ${D}/../../src_${V}_${K}_${TYPE}/${HEXE} . &
C
C --- summary printout
C
touch   summary_out
/bin/mv summary_out summary_old
C
C --- heat transport output
C
touch   flxdp_out.a flxdp_out.b
/bin/mv flxdp_out.a flxdp_old.a
/bin/mv flxdp_out.b flxdp_old.b
C
touch   ovrtn_out
/bin/mv ovrtn_out ovrtn_old
C
C --- clean up old archive files, typically from batch system rerun.
C
mkdir KEEP
touch archv.dummy.b
foreach f (arch*.{a,b,txt})
  /bin/mv $f KEEP/$f
end
C
C --- Nesting input archive files.
C
if (-e ./nest) then
  cd ./nest
  touch rmu.a rmu.b
  if (-z rmu.a) then
     ${pget} ${D}/../../relax/${E}/nest_rmu.a rmu.a &
  endif
  if (-z rmu.b) then
     ${pget} ${D}/../../relax/${E}/nest_rmu.b rmu.b &
  endif
  touch   arch.dummy.b
  /bin/rm -f arch*.[ab]
  touch  archv_${Y01}${A}.tar
  if (-z archv_${Y01}${A}.tar) then
    ${pget} ${D}/nest/archv_${Y01}${A}.tar archv_${Y01}${A}.tar
  endif
  tar xvf archv_${Y01}${A}.tar
  cd ..
endif
C
C --- let all file copies complete.
C
wait
C
C --- zero file length means no rivers.
C
if (-z forcing.rivers.a) then
   /bin/rm -f forcing.rivers.[ab]
endif
C
C --- Just in time atmospheric forcing.
C
if (-e ./mslp && ! -e ./wind) then
C
C --- mslprs forcing without wind forcing.
C --- Check to see if mslprs files exist, if not make them and wait.
C
  /bin/rm -f forcing.mslprs.a
  /bin/rm -f forcing.mslprs.b
  /bin/rm -f ./mslp/${E}m${Y01}${A}
  if (-e     ./mslp/mslprs_${Y01}${A}.a && \
      -e     ./mslp/mslprs_${Y01}${A}.b    ) then
    /bin/ln  ./mslp/mslprs_${Y01}${A}.a forcing.mslprs.a
    /bin/ln  ./mslp/mslprs_${Y01}${A}.b forcing.mslprs.b
  else
    cd ./mslp
    touch ${E}m${Y01}${A}
    /bin/rm -f ${E}m${Y01}${A}.com ${E}m${Y01}${A}.log
    awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}M.com > ${E}m${Y01}${A}.com
    csh ${E}m${Y01}${A}.com >& ${E}m${Y01}${A}.log
    cd ..
    /bin/ln  ./mslp/mslprs_${Y01}${A}.a forcing.mslprs.a
    /bin/ln  ./mslp/mslprs_${Y01}${A}.b forcing.mslprs.b
  endif
C
C --- If the mslprs for the next segment does not exist,
C --- interpolate them to model grid while current segment is running.
C
  if (-e ./mslp/mslprs_${YXX}${B}.a && \
      -e ./mslp/mslprs_${YXX}${B}.b    ) then
C
C --- next segments mslprs already exists.
C
  else
    cd ./mslp
    /bin/rm -f ${E}m${YXX}${B}.com ${E}m${YXX}${B}.log
    awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}M.com > ${E}m${YXX}${B}.com
    csh ${E}m${YXX}${B}.com >& ${E}m${YXX}${B}.log &
    cd ..
  endif
endif
C
if (-e ./wind) then
  if (! -e ./flux) then
    echo './flux must exist if ./wind does'
    exit
  endif
  if (-e ./wind/ZERO) then
C
C --- zero-length tau-wd and wndspd, use wnd-wd instead
C
    touch ./wind/tauewd_${Y01}${A}.a
    touch ./wind/tauewd_${Y01}${A}.b
    touch ./wind/taunwd_${Y01}${A}.a
    touch ./wind/taunwd_${Y01}${A}.b
    touch ./wind/wndspd_${Y01}${A}.a
    touch ./wind/wndspd_${Y01}${A}.b
    touch ./wind/tauewd_${YXX}${B}.a
    touch ./wind/tauewd_${YXX}${B}.b
    touch ./wind/taunwd_${YXX}${B}.a
    touch ./wind/taunwd_${YXX}${B}.b
    touch ./wind/wndspd_${YXX}${B}.a
    touch ./wind/wndspd_${YXX}${B}.b
  endif
C
C --- Check to see if wind and flux files exist, if not make them and wait.
C
  /bin/rm -f forcing.tauewd.a forcing.taunwd.a forcing.wndspd.a
  /bin/rm -f forcing.tauewd.b forcing.taunwd.b forcing.wndspd.b
  /bin/rm -f ./wind/${E}w${Y01}${A}
  if (-e     ./wind/tauewd_${Y01}${A}.a && \
      -e     ./wind/tauewd_${Y01}${A}.b && \
      -e     ./wind/taunwd_${Y01}${A}.a && \
      -e     ./wind/taunwd_${Y01}${A}.b    ) then
    /bin/ln  ./wind/tauewd_${Y01}${A}.a forcing.tauewd.a
    /bin/ln  ./wind/taunwd_${Y01}${A}.a forcing.taunwd.a
    /bin/ln  ./wind/tauewd_${Y01}${A}.b forcing.tauewd.b
    /bin/ln  ./wind/taunwd_${Y01}${A}.b forcing.taunwd.b
    if (-e     ./wind/wndspd_${Y01}${A}.a && \
        -e     ./wind/wndspd_${Y01}${A}.b    ) then
      /bin/ln  ./wind/wndspd_${Y01}${A}.a forcing.wndspd.a
      /bin/ln  ./wind/wndspd_${Y01}${A}.b forcing.wndspd.b
    endif
  else
    cd ./wind
    touch ${E}w${Y01}${A}
    /bin/rm -f ${E}w${Y01}${A}.com ${E}w${Y01}${A}.log
    awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}W.com > ${E}w${Y01}${A}.com
    csh ${E}w${Y01}${A}.com >& ${E}w${Y01}${A}.log &
    cd ..
  endif
  if (-e ./wspd) then
    /bin/rm -f forcing.wndspd.a
    /bin/rm -f forcing.wndspd.b
    /bin/rm -f ./wspd/${E}s${Y01}${A}
    if (-e     ./wspd/wndspd_${Y01}${A}.a && \
        -e     ./wspd/wndspd_${Y01}${A}.b    ) then
      /bin/ln  ./wspd/wndspd_${Y01}${A}.a forcing.wndspd.a
      /bin/ln  ./wspd/wndspd_${Y01}${A}.b forcing.wndspd.b
    else
      cd ./wspd
      touch ${E}s${Y01}${A}
      /bin/rm -f ${E}s${Y01}${A}.com ${E}s${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}S.com > ${E}s${Y01}${A}.com
      csh ${E}s${Y01}${A}.com >& ${E}s${Y01}${A}.log &
      cd ..
    endif
  endif
  if (-e ./wvel) then
    if (-e ./wind/ZERO) then
      /bin/rm -f forcing.wndewd.a forcing.wndnwd.a
      /bin/rm -f forcing.wndewd.b forcing.wndnwd.b
      /bin/rm -f ./wvel/${E}v${Y01}${A}
      if (-e     ./wvel/wndewd_${Y01}${A}.a && \
          -e     ./wvel/wndewd_${Y01}${A}.b && \
          -e     ./wvel/wndnwd_${Y01}${A}.a && \
          -e     ./wvel/wndnwd_${Y01}${A}.b    ) then
        /bin/ln  ./wvel/wndewd_${Y01}${A}.a forcing.wndewd.a
        /bin/ln  ./wvel/wndnwd_${Y01}${A}.a forcing.wndnwd.a
        /bin/ln  ./wvel/wndewd_${Y01}${A}.b forcing.wndewd.b
        /bin/ln  ./wvel/wndnwd_${Y01}${A}.b forcing.wndnwd.b
      else
        cd ./wvel
        touch ${E}v${Y01}${A}
        /bin/rm -f ${E}v${Y01}${A}.com ${E}v${Y01}${A}.log
        awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}V.com > ${E}v${Y01}${A}.com
        csh ${E}v${Y01}${A}.com >& ${E}v${Y01}${A}.log &
        cd ..
      endif
    else
      /bin/rm -f forcing.wndspd.a
      /bin/rm -f forcing.wndspd.b
      /bin/rm -f ./wvel/${E}v${Y01}${A}
      if (-e     ./wvel/wndspd_${Y01}${A}.a && \
          -e     ./wvel/wndspd_${Y01}${A}.b    ) then
        /bin/ln  ./wvel/wndspd_${Y01}${A}.a forcing.wndspd.a
        /bin/ln  ./wvel/wndspd_${Y01}${A}.b forcing.wndspd.b
      else
        cd ./wvel
        touch ${E}v${Y01}${A}
        /bin/rm -f ${E}v${Y01}${A}.com ${E}v${Y01}${A}.log
        awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}V.com > ${E}v${Y01}${A}.com
        csh ${E}v${Y01}${A}.com >& ${E}v${Y01}${A}.log &
        cd ..
      endif
    endif
  endif
  if (-e ./grad) then
    /bin/rm -f ./grad/${E}g${Y01}${A}
    if (-e     ./grad/glbrad_${Y01}${A}.a && \
        -e     ./grad/glbrad_${Y01}${A}.b    ) then
C
C     this segments glbrad already exists
C
      touch ./grad/${E}g${Y01}${A}
    else
      cd ./grad
      touch ${E}g${Y01}${A}
      /bin/rm -f ${E}g${Y01}${A}.com ${E}g${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}G.com > ${E}g${Y01}${A}.com
      csh ${E}g${Y01}${A}.com >& ${E}g${Y01}${A}.log &
      cd ..
    endif
  endif
  if (-e ./lrad) then
    /bin/rm -f forcing.radflx.a
    /bin/rm -f forcing.radflx.b
    /bin/rm -f forcing.shwflx.a
    /bin/rm -f forcing.shwflx.b
    /bin/rm -f ./lrad/${E}l${Y01}${A}
    if (-e     ./lrad/lwdflx_${Y01}${A}.a && \
        -e     ./lrad/lwdflx_${Y01}${A}.b    ) then
C
C     this segments lwdflx already exists
C
      touch ./lrad/${E}l${Y01}${A}
    else
      cd ./lrad
      touch ${E}l${Y01}${A}
      /bin/rm -f ${E}l${Y01}${A}.com ${E}l${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}L.com > ${E}l${Y01}${A}.com
      csh ${E}l${Y01}${A}.com >& ${E}l${Y01}${A}.log &
      cd ..
    endif
  endif
  if (-e ./mslp) then
    /bin/rm -f forcing.mslprs.a
    /bin/rm -f forcing.mslprs.b
    /bin/rm -f ./mslp/${E}m${Y01}${A}
    if (-e     ./mslp/mslprs_${Y01}${A}.a && \
        -e     ./mslp/mslprs_${Y01}${A}.b    ) then
      /bin/ln  ./mslp/mslprs_${Y01}${A}.a forcing.mslprs.a
      /bin/ln  ./mslp/mslprs_${Y01}${A}.b forcing.mslprs.b
    else
      cd ./mslp
      touch ${E}m${Y01}${A}
      /bin/rm -f ${E}m${Y01}${A}.com ${E}m${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}M.com > ${E}m${Y01}${A}.com &
      csh ${E}m${Y01}${A}.com >& ${E}m${Y01}${A}.log
      cd ..
    endif
  endif
  if (-e ./ssta) then
    /bin/rm -f forcing.surtmp.a
    /bin/rm -f forcing.surtmp.b
    /bin/rm -f ./ssta/${E}p${Y01}${A}
    if (-e     ./ssta/surtmp_${Y01}${A}.a && \
        -e     ./ssta/surtmp_${Y01}${A}.b    ) then
      /bin/ln  ./ssta/surtmp_${Y01}${A}.a forcing.surtmp.a
      /bin/ln  ./ssta/surtmp_${Y01}${A}.b forcing.surtmp.b
    else
      cd ./ssta
      touch ${E}t${Y01}${A}
      /bin/rm -f ${E}t${Y01}${A}.com ${E}t${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}T.com > ${E}t${Y01}${A}.com
      csh ${E}t${Y01}${A}.com >& ${E}t${Y01}${A}.log &
      cd ..
    endif
  endif
  if (-e ./ssto) then
    /bin/rm -f forcing.seatmp.a
    /bin/rm -f forcing.seatmp.b
    /bin/rm -f ./ssto/${E}p${Y01}${A}
    if (-e     ./ssto/seatmp_${Y01}${A}.a && \
        -e     ./ssto/seatmp_${Y01}${A}.b    ) then
      /bin/ln  ./ssto/seatmp_${Y01}${A}.a forcing.seatmp.a
      /bin/ln  ./ssto/seatmp_${Y01}${A}.b forcing.seatmp.b
    else
      cd ./ssto
      touch ${E}o${Y01}${A}
      /bin/rm -f ${E}o${Y01}${A}.com ${E}o${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}O.com > ${E}o${Y01}${A}.com
      csh ${E}o${Y01}${A}.com >& ${E}o${Y01}${A}.log &
      cd ..
    endif
  endif
  if (-e ./pcip) then
    /bin/rm -f forcing.precip.a
    /bin/rm -f forcing.precip.b
    /bin/rm -f ./pcip/${E}p${Y01}${A}
    if (-e     ./pcip/precip_${Y01}${A}.a && \
        -e     ./pcip/precip_${Y01}${A}.b    ) then
      /bin/ln  ./pcip/precip_${Y01}${A}.a forcing.precip.a
      /bin/ln  ./pcip/precip_${Y01}${A}.b forcing.precip.b
    else
      cd ./pcip
      touch ${E}p${Y01}${A}
      /bin/rm -f ${E}p${Y01}${A}.com ${E}p${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}P.com > ${E}p${Y01}${A}.com
      csh ${E}p${Y01}${A}.com >& ${E}p${Y01}${A}.log &
      cd ..
    endif
  endif
  if (-e ./flxt) then
    /bin/rm -f forcing.radflx.a
    /bin/rm -f forcing.radflx.b
    /bin/rm -f ./flxt/${E}p${Y01}${A}
    if (-e     ./flxt/totflx_${Y01}${A}.a && \
        -e     ./flxt/totflx_${Y01}${A}.b    ) then
      /bin/ln  ./flxt/totflx_${Y01}${A}.a forcing.radflx.a
      /bin/ln  ./flxt/totflx_${Y01}${A}.b forcing.radflx.b
    else
      cd ./flxt
      touch ${E}q${Y01}${A}
      /bin/rm -f ${E}p${Y01}${A}.com ${E}p${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}Q.com > ${E}q${Y01}${A}.com
      csh ${E}q${Y01}${A}.com >& ${E}q${Y01}${A}.log &
      cd ..
    endif
  endif
  if (-e ./flux) then
    /bin/rm -f forcing.airtmp.a forcing.shwflx.a forcing.vapmix.a
    /bin/rm -f forcing.airtmp.b forcing.shwflx.b forcing.vapmix.b
    if (! -e ./pcip) then
      /bin/rm -f forcing.precip.a
      /bin/rm -f forcing.precip.b
      touch ./flux/precip_${Y01}${A}.a
      touch ./flux/precip_${Y01}${A}.b
    endif
    if (! -e ./flxt && ! -e ./lrad) then
      /bin/rm -f forcing.radflx.a
      /bin/rm -f forcing.radflx.b
      touch ./flux/radflx_${Y01}${A}.a
      touch ./flux/radflx_${Y01}${A}.b
    endif
    if (! -e ./lrad) then
      touch ./flux/shwflx_${Y01}${A}.a
      touch ./flux/shwflx_${Y01}${A}.b
    endif
    /bin/rm -f ./flux/${E}f${Y01}${A}
    if (-e     ./flux/airtmp_${Y01}${A}.a && \
        -e     ./flux/airtmp_${Y01}${A}.b && \
        -e     ./flux/vapmix_${Y01}${A}.a && \
        -e     ./flux/vapmix_${Y01}${A}.b    ) then
      /bin/ln  ./flux/airtmp_${Y01}${A}.a forcing.airtmp.a
      /bin/ln  ./flux/airtmp_${Y01}${A}.b forcing.airtmp.b
      /bin/ln  ./flux/vapmix_${Y01}${A}.a forcing.vapmix.a
      /bin/ln  ./flux/vapmix_${Y01}${A}.b forcing.vapmix.b
      if (! -e ./pcip) then
        /bin/ln ./flux/precip_${Y01}${A}.b forcing.precip.b
        /bin/ln ./flux/precip_${Y01}${A}.a forcing.precip.a
      endif
      if (! -e ./flxt && ! -e ./lrad) then
        /bin/ln ./flux/radflx_${Y01}${A}.a forcing.radflx.a
        /bin/ln ./flux/radflx_${Y01}${A}.b forcing.radflx.b
      endif
      if (! -e ./lrad) then
        /bin/ln ./flux/shwflx_${Y01}${A}.a forcing.shwflx.a
        /bin/ln ./flux/shwflx_${Y01}${A}.b forcing.shwflx.b
      endif
    else
      cd ./flux
      touch ${E}f${Y01}${A}
      /bin/rm -f ${E}f${Y01}${A}.com ${E}f${Y01}${A}.log
      awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}F.com > ${E}f${Y01}${A}.com
      csh ${E}f${Y01}${A}.com >& ${E}f${Y01}${A}.log &
      cd ..
    endif
  endif
  wait
  if (-e    ./wind/${E}w${Y01}${A}) then
    /bin/ln ./wind/tauewd_${Y01}${A}.a forcing.tauewd.a
    /bin/ln ./wind/taunwd_${Y01}${A}.a forcing.taunwd.a
    /bin/ln ./wind/tauewd_${Y01}${A}.b forcing.tauewd.b
    /bin/ln ./wind/taunwd_${Y01}${A}.b forcing.taunwd.b
    if (-e    ./wind/wndspd_${Y01}${A}.a && \
        -e    ./wind/wndspd_${Y01}${A}.b    ) then
      /bin/ln ./wind/wndspd_${Y01}${A}.a forcing.wndspd.a
      /bin/ln ./wind/wndspd_${Y01}${A}.b forcing.wndspd.b
    endif
  endif
  if (-e ./wvel) then
    if (-e ./wvel/${E}v${Y01}${A}) then
      if (-e ./wind/ZERO) then
        /bin/ln  ./wvel/wndewd_${Y01}${A}.a forcing.wndewd.a
        /bin/ln  ./wvel/wndnwd_${Y01}${A}.a forcing.wndnwd.a
        /bin/ln  ./wvel/wndewd_${Y01}${A}.b forcing.wndewd.b
        /bin/ln  ./wvel/wndnwd_${Y01}${A}.b forcing.wndnwd.b
      else
        /bin/ln ./wvel/wndspd_${Y01}${A}.a forcing.wndspd.a
        /bin/ln ./wvel/wndspd_${Y01}${A}.b forcing.wndspd.b
      endif
    endif
  endif
  if (-e ./wspd) then
    if (-e ./wspd/${E}s${Y01}${A}) then
      /bin/ln ./wspd/wndspd_${Y01}${A}.a forcing.wndspd.a
      /bin/ln ./wspd/wndspd_${Y01}${A}.b forcing.wndspd.b
    endif
  endif
  if (-e ./mslp) then
    if (-e ./mslp/${E}m${Y01}${A}) then
      /bin/ln ./mslp/mslprs_${Y01}${A}.a forcing.mslprs.a
      /bin/ln ./mslp/mslprs_${Y01}${A}.b forcing.mslprs.b
    endif
  endif
  if (-e ./ssta) then
    if (-e ./ssta/${E}t${Y01}${A}) then
      /bin/ln ./ssta/surtmp_${Y01}${A}.a forcing.surtmp.a
      /bin/ln ./ssta/surtmp_${Y01}${A}.b forcing.surtmp.b
    endif
  endif
  if (-e ./ssto) then
    if (-e ./ssto/${E}t${Y01}${A}) then
      /bin/ln ./ssto/seatmp_${Y01}${A}.a forcing.seatmp.a
      /bin/ln ./ssto/seatmp_${Y01}${A}.b forcing.seatmp.b
    endif
  endif
  if (-e ./pcip) then
    if (-e ./pcip/${E}p${Y01}${A}) then
      /bin/ln ./pcip/precip_${Y01}${A}.a forcing.precip.a
      /bin/ln ./pcip/precip_${Y01}${A}.b forcing.precip.b
    endif
  endif
  if (-e ./lrad) then
    if (-e ./grad/${E}g${Y01}${A}) then
      /bin/ln ./grad/glbrad_${Y01}${A}.a forcing.shwflx.a
      /bin/ln ./grad/glbrad_${Y01}${A}.b forcing.shwflx.b
    endif
    if (-e ./lrad/${E}l${Y01}${A}) then
      /bin/ln ./lrad/lwdflx_${Y01}${A}.a forcing.radflx.a
      /bin/ln ./lrad/lwdflx_${Y01}${A}.b forcing.radflx.b
    endif
  endif
  if (-e ./flxt) then
    if (-e ./flxt/${E}q${Y01}${A}) then
      /bin/ln ./flxt/totflx_${Y01}${A}.a forcing.radflx.a
      /bin/ln ./flxt/totflx_${Y01}${A}.b forcing.radflx.b
    endif
  endif
  if (-e ./flux) then
    if (-e    ./flux/${E}f${Y01}${A}) then
      /bin/ln ./flux/airtmp_${Y01}${A}.a forcing.airtmp.a
      /bin/ln ./flux/airtmp_${Y01}${A}.b forcing.airtmp.b
      /bin/ln ./flux/vapmix_${Y01}${A}.a forcing.vapmix.a
      /bin/ln ./flux/vapmix_${Y01}${A}.b forcing.vapmix.b
      if (! -e ./pcip) then
        /bin/ln ./flux/precip_${Y01}${A}.b forcing.precip.b
        /bin/ln ./flux/precip_${Y01}${A}.a forcing.precip.a
      endif
      if (! -e ./flxt && ! -e ./lrad) then
        /bin/ln ./flux/radflx_${Y01}${A}.a forcing.radflx.a
        /bin/ln ./flux/radflx_${Y01}${A}.b forcing.radflx.b
      endif
      if (! -e ./lrad) then
        /bin/ln ./flux/shwflx_${Y01}${A}.a forcing.shwflx.a
        /bin/ln ./flux/shwflx_${Y01}${A}.b forcing.shwflx.b
      endif
    endif
  endif
C
C --- CICE.
C
  if (-e ./cice) then
    if (-e ./cice/lwdflx_${Y01}${A}.r) then
C
C --- this segments lwdflx already exists.
C
    else
      cd ./cice
     switch ($OS)
      case 'XT4':
      case 'XC30':
        setenv SRC      ~wallcraf/hycom/ALLcnl/bin
        setenv APRUN    "aprun -n 1"
        breaksw
      default:
        setenv SRC      ~wallcraf/hycom/ALL/bin
        setenv APRUN    ""
      endsw
      setenv IDM 258
      setenv JDM 175
      setenv JDA 175
      if (-e ../lrad) then
        /bin/rm -f lwdflx_${Y01}${A}.?
        /bin/ln ../lrad/lwdflx_${Y01}${A}.a lwdflx_${Y01}${A}.A
        /bin/cp ../lrad/lwdflx_${Y01}${A}.b lwdflx_${Y01}${A}.B
      else
        /bin/rm -f lwdflx_${Y01}${A}.?
# ---   CICE: emissivity=0.95, StefanBoltzman=567.e-10, 0.95*567.e-10=538.65e-10
        ${APRUN} ${SRC}/hycom_expr netQlw_${Y01}${A}.a surtmp4_${Y01}${A}.a ${IDM} ${JDM} 1.0 538.65e-10 lwdflx_${Y01}${A}.A >! lwdflx_${Y01}${A}.B
# ---   NWP: emissivity=1.00, StefanBoltzman=567.e-10
########${APRUN} ${SRC}/hycom_expr netQlw_${Y01}${A}.a surtmp4_${Y01}${A}.a ${IDM} ${JDM} 1.0 567.00e-10 lwdflx_${Y01}${A}.A >! lwdflx_${Y01}${A}.B
      endif
      touch  ../forcing.offlux.a
      if (-z ../forcing.offlux.a) then
        /bin/mv -f lwdflx_${Y01}${A}.A lwdflx_${Y01}${A}.a
        /bin/mv -f lwdflx_${Y01}${A}.B lwdflx_${Y01}${A}.b
      else
# ---   add offlux to lwdflx.
        ${APRUN} ${SRC}/hycom_expr lwdflx_${Y01}${A}.A ../forcing.offlux.a ${IDM} ${JDM} 1.0 1.0 lwdflx_${Y01}${A}.a repeat >! lwdflx_${Y01}${A}.b
      endif
      ${APRUN} ${SRC}/hycom2raw8 lwdflx_${Y01}${A}.a ${IDM} ${JDM} 1 1 ${IDM} ${JDA} lwdflx_${Y01}${A}.r >! lwdflx_${Y01}${A}.B
      if (-e ./SAVE) then
        foreach f ( lwdflx_${Y01}${A}.[rB] )
          ln ${f} ./SAVE/${f}
        end
      endif
      cd ..
    endif
    /bin/rm cice/*${Y01}${A}.[Aab]
    if (-e ./grad) then
      /bin/rm -f                         cice.glbrad.r
      /bin/ln ./cice/glbrad_${Y01}${A}.r cice.glbrad.r
    else
      /bin/rm -f                         cice.netrad.r
      /bin/ln ./cice/netrad_${Y01}${A}.r cice.netrad.r
    endif
    foreach t ( airtmp lwdflx vapmix wndewd wndnwd )
      /bin/rm -f                       cice.${t}.r
      /bin/ln ./cice/${t}_${Y01}${A}.r cice.${t}.r
    end
  endif
C
C --- If the winds or fluxes for the next segment dont exist, 
C --- interpolate them to model grid while current segment is running.
C
  if (-e ./wind/tauewd_${YXX}${B}.a && \
      -e ./wind/tauewd_${YXX}${B}.b && \
      -e ./wind/taunwd_${YXX}${B}.a && \
      -e ./wind/taunwd_${YXX}${B}.b    ) then
C
C --- next segments winds already exist.
C
  else
    cd ./wind
    /bin/rm -f ${E}w${YXX}${B}.com ${E}w${YXX}${B}.log
    awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}W.com > ${E}w${YXX}${B}.com
    csh ${E}w${YXX}${B}.com >& ${E}w${YXX}${B}.log &
    cd ..
  endif
  if (-e ./wspd) then
    if (-e ./wspd/wndspd_${YXX}${B}.a && \
        -e ./wspd/wndspd_${YXX}${B}.b    ) then
C
C ---   next segments wndspd already exists.
C
    else
      cd ./wspd
      /bin/rm -f ${E}s${YXX}${B}.com ${E}s${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}S.com > ${E}s${YXX}${B}.com
      csh ${E}s${YXX}${B}.com >& ${E}s${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./wvel) then
    if (-e ./wind/ZERO) then
      if (-e     ./wvel/wndewd_${YXX}${B}.a && \
          -e     ./wvel/wndewd_${YXX}${B}.b && \
          -e     ./wvel/wndnwd_${YXX}${B}.a && \
          -e     ./wvel/wndnwd_${YXX}${B}.b    ) then
C
C ---     next segments wnd-wd already exists.
C
      else
        cd ./wvel
        /bin/rm -f ${E}v${YXX}${B}.com ${E}v${YXX}${B}.log
        awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}V.com > ${E}v${YXX}${B}.com
        csh ${E}v${YXX}${B}.com >& ${E}v${YXX}${B}.log &
        cd ..
      endif
    else
      if (-e ./wvel/wndspd_${YXX}${B}.a && \
          -e ./wvel/wndspd_${YXX}${B}.b    ) then
C
C ---     next segments wndspd already exists.
C
      else
        cd ./wvel
        /bin/rm -f ${E}v${YXX}${B}.com ${E}v${YXX}${B}.log
        awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}V.com > ${E}v${YXX}${B}.com
        csh ${E}v${YXX}${B}.com >& ${E}v${YXX}${B}.log &
        cd ..
      endif
    endif
  endif
  if (-e ./grad) then
    if (-e ./grad/glbrad_${YXX}${B}.a && \
        -e ./grad/glbrad_${YXX}${B}.b    ) then
C
C ---   next segments glbrad already exists.
C
    else
      cd ./grad
      /bin/rm -f ${E}g${YXX}${B}.com ${E}g${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}G.com > ${E}g${YXX}${B}.com
      csh ${E}g${YXX}${B}.com >& ${E}g${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./lrad) then
    if (-e ./lrad/lwdflx_${YXX}${B}.a && \
        -e ./lrad/lwdflx_${YXX}${B}.b    ) then
C
C ---   next segments lwdflx already exists.
C
    else
      cd ./lrad
      /bin/rm -f ${E}l${YXX}${B}.com ${E}l${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}L.com > ${E}l${YXX}${B}.com
      csh ${E}l${YXX}${B}.com >& ${E}l${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./mslp) then
    if (-e ./mslp/mslprs_${YXX}${B}.a && \
        -e ./mslp/mslprs_${YXX}${B}.b    ) then
C
C ---   next segments mslprs already exists.
C
    else
      cd ./mslp
      /bin/rm -f ${E}m${YXX}${B}.com ${E}m${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}M.com > ${E}m${YXX}${B}.com
      csh ${E}m${YXX}${B}.com >& ${E}m${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./ssta) then
    if (-e ./ssta/surtmp_${YXX}${B}.a && \
        -e ./ssta/surtmp_${YXX}${B}.b    ) then
C
C ---   next segments surtmp already exists.
C
    else
      cd ./ssta
      /bin/rm -f ${E}t${YXX}${B}.com ${E}t${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}T.com > ${E}t${YXX}${B}.com
      csh ${E}t${YXX}${B}.com >& ${E}t${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./ssto) then
    if (-e ./ssto/seatmp_${YXX}${B}.a && \
        -e ./ssto/seatmp_${YXX}${B}.b    ) then
C
C ---   next segments seatmp already exists.
C
    else
      cd ./ssto
      /bin/rm -f ${E}o${YXX}${B}.com ${E}o${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}O.com > ${E}o${YXX}${B}.com
      csh ${E}o${YXX}${B}.com >& ${E}o${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./pcip) then
    if (-e ./pcip/precip_${YXX}${B}.a && \
        -e ./pcip/precip_${YXX}${B}.b    ) then
C
C ---   next segments pcip already exists.
C
    else
      cd ./pcip
      /bin/rm -f ${E}p${YXX}${B}.com ${E}p${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}P.com > ${E}p${YXX}${B}.com
      csh ${E}p${YXX}${B}.com >& ${E}p${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./flxt) then
    if (-e ./flxt/totflx_${YXX}${B}.a && \
        -e ./flxt/totflx_${YXX}${B}.b    ) then
C
C ---   next segments flxt already exist.
C
    else
      cd ./flxt
      /bin/rm -f ${E}q${YXX}${B}.com ${E}q${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}Q.com > ${E}q${YXX}${B}.com
      csh ${E}q${YXX}${B}.com >& ${E}q${YXX}${B}.log &
      cd ..
    endif
  endif
  if (-e ./flux) then
    if (! -e ./pcip) then
      touch ./flux/precip_${YXX}${B}.a
      touch ./flux/precip_${YXX}${B}.b
    endif
    if (! -e ./flxt) then
      touch ./flux/radflx_${YXX}${B}.a
      touch ./flux/radflx_${YXX}${B}.b
    endif
    if (-e ./lrad) then
      touch ./flux/radflx_${YXX}${B}.a
      touch ./flux/radflx_${YXX}${B}.b
      touch ./flux/shwflx_${YXX}${B}.a
      touch ./flux/shwflx_${YXX}${B}.b
    endif
    if (-e ./flux/airtmp_${YXX}${B}.a && \
        -e ./flux/airtmp_${YXX}${B}.b && \
        -e ./flux/precip_${YXX}${B}.a && \
        -e ./flux/precip_${YXX}${B}.b && \
        -e ./flux/radflx_${YXX}${B}.a && \
        -e ./flux/radflx_${YXX}${B}.b && \
        -e ./flux/shwflx_${YXX}${B}.a && \
        -e ./flux/shwflx_${YXX}${B}.b && \
        -e ./flux/vapmix_${YXX}${B}.a && \
        -e ./flux/vapmix_${YXX}${B}.b    ) then
C
C ---   next segments fluxes already exist.
C
    else
      cd ./flux
      /bin/rm -f ${E}f${YXX}${B}.com ${E}f${YXX}${B}.log
      awk -f $D/../${E}.awk y01=${YXX} ab=${B} $D/../${E}F.com > ${E}f${YXX}${B}.com
      csh ${E}f${YXX}${B}.com >& ${E}f${YXX}${B}.log &
      cd ..
    endif
  endif
endif
C
C --- Nesting input archive files for next segment.
C
if (-e ./nest) then
  cd ./nest
  touch  archv_${YXX}${B}.tar
  if (-z archv_${YXX}${B}.tar) then
    ${pget} ${D}/nest/archv_${YXX}${B}.tar archv_${YXX}${B}.tar &
  endif
  cd ..
endif
C
chmod ug+x ${HEXE}
/bin/ls -laFq
C
if (-e ./nest) then
  ls -laFq nest
endif
C
C ---  Check to make sure restart file is there
C
if (`echo $LI | awk '{print ($1 > 0.0)}'` && -z restart_in.a) then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
endif
C
if ($NMPI == 0) then
C
C --- run the model, without MPI or SHMEM
C
if ($NOMP == 0) then
  setenv NOMP 1
endif
C
switch ($OS)
case 'SunOS':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    env OMP_NUM_THREADS=$NOMP ./${HEXE}
    breaksw
case 'IDP':
case 'ICE':
case 'Linux':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    env OMP_NUM_THREADS=$NOMP MPSTKZ=8M ./${HEXE}
    breaksw
case 'OSF1':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    env OMP_NUM_THREADS=$NOMP ./${HEXE}
    breaksw
case 'IRIX64':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv FILENV .assign
    assign -R
    assign -s sbin u:18
    assign -V
    env OMP_NUM_THREADS=$NOMP ./${HEXE}
    assign -V
    assign -R
    breaksw
case 'AIX':
C
C   --- $NOMP CPUs/THREADs, if compiled for IBM OpenMP.
C
    /bin/rm -f core
    touch core
    setenv SPINLOOPTIME     500
    setenv YIELDLOOPTIME    500
    setenv XLSMPOPTS       "parthds=${NOMP}:spins=0:yields=0"
    ./${HEXE}
    breaksw
#case 'AIX':
#C
#C   --- $NOMP CPUs/THREADs, if compiled for KAI OpenMP.
#C
#    /bin/rm -f core
#    touch core
#    env OMP_NUM_THREADS=$NOMP ./${HEXE}
#    breaksw
case 'unicos':
C
C   --- $NOMP CPUs/THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    assign -V
    env OMP_NUM_THREADS=$NOMP ./${HEXE}
    if (! -z core)  debug -s ${HEXE} core
    assign -V
    assign -R
    breaksw
endsw
else
C
C --- run the model, with MPI or SHMEM and perhaps also with OpenMP.
C
touch patch.input
if (-z patch.input) then
C
C --- patch.input is always required for MPI or SHMEM.
C
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
  exit
endif
C
switch ($OS)
case 'SunOS':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv OMP_NUM_THREADS $NOMP
#   mpirun -np $NMPI ./${HEXE}
    pam ./${HEXE}
    breaksw
case 'Linux':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv OMP_NUM_THREADS $NOMP
    mpirun -np $NMPI ./${HEXE}
    breaksw
case 'ICE':
    if ($NOMP == 0) then
        limit stacksize unlimited
        setenv MPI_DSM_DISTRIBUTE   yes
        setenv MPI_BUFS_PER_HOST    768
        setenv MPI_BUFS_PER_PROC    128
        setenv MPI_GROUP_MAX        128
        time mpiexec_mpt -np $NMPI ./${HEXE}
    else
        limit stacksize unlimited
        setenv OMP_NUM_THREADS      $NOMP
        setenv OMP_STACKSIZE        127M
        setenv MPI_DSM_DISTRIBUTE   yes
        setenv MPI_BUFS_PER_HOST    768
        setenv MPI_BUFS_PER_PROC    128
        setenv MPI_GROUP_MAX        128
        time mpiexec_mpt -np $NMPI omplace -nt $NOMP ./${HEXE}
    endif
    breaksw
case 'IDP':
# ---   debugging
#   setenv I_MPI_DEBUG                      5
# ---   From "Using IntelMPI on Discover"
# ---   https://modelingguru.nasa.gov/docs/DOC-1670
    setenv I_MPI_DAPL_SCALABLE_PROGRESS     1
    setenv I_MPI_DAPL_RNDV_WRITE            1
    setenv I_MPI_JOB_STARTUP_TIMEOUT        10000
    setenv I_MPI_HYDRA_BRANCH_COUNT         512
    setenv DAPL_ACK_RETRY            7
    setenv DAPL_ACK_TIMER            23
    setenv DAPL_RNR_RETRY            7
    setenv DAPL_RNR_TIMER            28
# ---   intel scaling suggestions
    setenv DAPL_CM_ARP_TIMEOUT_MS    8000
    setenv DAPL_CM_ARP_RETRY_COUNT   25
    setenv DAPL_CM_ROUTE_TIMEOUT_MS  20000
    setenv DAPL_CM_ROUTE_RETRY_COUNT 15
    setenv DAPL_MAX_CM_RESPONSE_TIME 20
    setenv DAPL_MAX_CM_RETRIES       15
    if ($NOMP == 0) then
        setenv OMP_NUM_THREADS      1
        mpirun ./${HEXE}
    else
        setenv OMP_NUM_THREADS      $NOMP
        mpirun ./${HEXE}
    endif
    breaksw
case 'XT3':
C
C   --- $NMPI MPI tasks.
C
    /bin/rm -f core
    touch core
    setenv NO_STOP_MESSAGE
#   work-around for a lustre bug: pause for 20 seconds
### sleep 20
    setenv MPICH_RANK_REORDER_METHOD	1
    setenv MPI_COLL_OPT_ON		1
#   setenv IOBUF_PARAMS '%stdout,%stderr,*'
    sleep 120
    time yod -small_pages -np $NMPI ./${HEXE}
    breaksw
case 'XT4':
C
C   --- $NMPI MPI tasks.
C
    /bin/rm -f core
    touch core
    setenv NO_STOP_MESSAGE
#   work-around for a lustre bug: pause for 20 seconds
#   sleep 20
    setenv MPICH_RANK_REORDER_METHOD	1
    setenv MPI_COLL_OPT_ON		1
    time aprun -n $NMPI ./${HEXE}
    breaksw
case 'OSF1':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv OMP_NUM_THREADS $NOMP
#   mpirun -np $NMPI ./${HEXE}
    time prun -n $NMPI ./${HEXE}
    breaksw
case 'IRIX64':
if ($TYPE == "shmem") then
C
C   --- $NMPI SHMEM tasks
C
    /bin/rm -f core
    touch core
    setenv FILENV .assign
    assign -R
    assign -s sbin u:18
    assign -V
    setenv OMP_NUM_THREADS	1
    setenv SMA_DSM_TOPOLOGY	free
    setenv SMA_DSM_VERBOSE	1
    setenv SMA_VERSION		1
    env NPES=$NMPI ./${HEXE}
    assign -V
    assign -R
    breaksw
else
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for OpenMP.
C
    /bin/rm -f core
    touch core
    setenv FILENV .assign
    assign -R
    assign -s sbin u:18
    assign -V
    setenv OMP_NUM_THREADS	$NOMP
    setenv MPI_DSM_VERBOSE	1
    setenv MPI_REQUEST_MAX	8192
    mpirun -np $NMPI ./${HEXE}
#   mpirun -np $NMPI ./${HEXE} < /dev/null
    assign -V
    assign -R
    breaksw
endif
case 'AIX':
C
C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for IBM OpenMP.
C
    /bin/rm -f core
    touch core
    setenv SPINLOOPTIME		500
    setenv YIELDLOOPTIME	500
    setenv XLSMPOPTS		"parthds=${NOMP}:spins=0:yields=0"
    setenv MP_SHARED_MEMORY	yes
    setenv MP_SINGLE_THREAD	yes
#   setenv MP_SINGLE_THREAD	no
    setenv MP_EAGER_LIMIT	65536
#   setenv MP_EUILIB		us
#   list where the MPI job will run
#   env MP_LABELIO=YES $POE hostname
    if      (-e /site/bin/launch) then
      setenv MEMORY_AFFINITY    MCM
      setenv UL_MODE            PURE_MPI
      setenv UL_TARGET_CPU_LIST AUTO_SELECT
      time $POE /site/bin/launch ./${HEXE}
    else
      time $POE ./${HEXE}
    endif
    breaksw
#case 'AIX':
#C
#C   --- $NMPI MPI tasks and $NOMP THREADs, if compiled for KAI OpenMP.
#C
#    /bin/rm -f core
#    touch core
#    setenv OMP_NUM_THREADS	$NOMP
#    setenv MP_SHARED_MEMORY	yes
#    setenv MP_SINGLE_THREAD	yes
#    setenv MP_EAGER_LIMIT	65536
#    setenv MP_EUILIB		us
#    setenv MP_EUIDEVICE		css0
##   list where the MPI job will run
#    env MP_LABELIO=YES $POE hostname
#    time $POE ./${HEXE}
#    breaksw
default:
    echo "This O/S," $OS ", is not configured for MPI/SHMEM"
    exit (1)
endsw
endif
C
touch   PIPE_DEBUG
/bin/rm PIPE_DEBUG
C
C --- archive output in a separate tar directory
C
touch archv.dummy.a archv.dummy.b archv.dummy.txt
touch archs.dummy.a archs.dummy.b archs.dummy.txt
touch archm.dummy.a archm.dummy.b archm.dummy.txt
touch arche.dummy.a arche.dummy.b arche.dummy.txt
touch cice.dummy.nc
C
if (-e ./SAVE) then
  foreach t ( v s m e )
    foreach f (arch${t}.*.a)
      /bin/ln ${f} SAVE/${f}
    end
    foreach f (arch${t}.*.b)
      /bin/ln ${f} SAVE/${f}
    end
    foreach f (arch${t}.*.txt)
      /bin/ln ${f} SAVE/${f}
    end
  end
  foreach f (cice.*.nc)
    /bin/ln -f ${f} SAVE/${f}
  end
endif
C
foreach t ( v s m e )
  mkdir ./tar${t}_${Y01}${A}
switch ($OS)
case 'XT3':
case 'XT4':
  lfs setstripe ./tar${t}_${Y01}${A} 1048576 -1 8
  breaksw
endsw
  foreach f (arch${t}.*.a)
    /bin/mv ${f} ./tar${t}_${Y01}${A}/${E}_${f}
  end
  foreach f (arch${t}.*.b)
    /bin/mv ${f} ./tar${t}_${Y01}${A}/${E}_${f}
  end
  foreach f (arch${t}.*.txt)
    /bin/mv ${f} ./tar${t}_${Y01}${A}/${E}_${f}
  end
  date
end
foreach f (cice.*.nc)
  /bin/mv ${f} ./tarc_${Y01}${A}/${E}_${f}
end
C 
if (! -z archt.input) then
  if (-e ./tart_${Y01}${A}) then
    /bin/mv ./tart_${Y01}${A} ./tart_${Y01}${A}_$$
  endif
  /bin/mv ./ARCHT ./tart_${Y01}${A}
endif
C
C --- add last day to next months tar directory, for actual day months only
C
setenv DL `awk  '{printf("%15.2f\n", $2)}' limits`
setenv DA `echo 3 1.0 1.0 $DL $DL | ~/hycom/ALL/bin/hycom_nest_dates | head -1`
foreach t ( v s c )
  mkdir ./tar${t}_${YXX}${B}
  ln -f ./tar${t}_${Y01}${A}/${E}_arch?.${DA}.* ./tar${t}_${YXX}${B}
end
C
C --- build and run or submit the tar script
C
awk -f $D/../${E}.awk y01=${Y01} ab=${A} $D/../${E}A.com >! \
       tar_${Y01}${A}.com
csh tar_${Y01}${A}.com >&! tar_${Y01}${A}.log &
#llsubmit ./tar_${Y01}${A}.com
#~wallcraf/bin/q_navo tar_${Y01}${A}.com
C
C --- heat transport statistics output
C
if (-e flxdp_out.a) then
  ${pput} flxdp_out.a ${S}/flxdp_${Y01}${A}.a
endif
if (-e flxdp_out.b) then
  ${pput} flxdp_out.b ${S}/flxdp_${Y01}${A}.b
endif
if (-e ovrtn_out) then
  ${pput} ovrtn_out ${S}/ovrtn_${Y01}${A}
endif
C
C --- restart output
C
if (-e restart_out.a) then
  ${pput} restart_out.a ${S}/restart_${YXX}${B}.a
endif
if (-e restart_out.b) then
  ${pput} restart_out.b ${S}/restart_${YXX}${B}.b
endif
endif
if (-e ./cice) then
C
C --- CICE restart output, assumes single-month runs
C
  /bin/mv cice.restart*01 cice.restart_out
  if (-e cice.restart_out) then
    ${pput} cice.restart_out ${S}/cice.restart_${YXX}${B}
  endif
endif
C
if (-e ./wind) then
#C
#C --- Delete just in time wind and flux files.
#C
#  touch summary_out
#  tail -1 summary_out
#  tail -1 summary_out | grep -c "^normal stop"
#  if ( `tail -1 summary_out | grep -c "^normal stop"` == 1 ) then
#    /bin/rm -f ./wind/*_${Y01}${A}.[ab]
#    /bin/rm -f ./wspd/*_${Y01}${A}.[ab]
#    /bin/rm -f ./wvel/*_${Y01}${A}.[ab]
#    /bin/rm -f ./grad/*_${Y01}${A}.[ab]
#    /bin/rm -f ./flux/*_${Y01}${A}.[ab]
#    /bin/rm -f ./flxt/*_${Y01}${A}.[ab]
#    /bin/rm -f ./lrad/*_${Y01}${A}.[ab]
#    /bin/rm -f ./mslp/*_${Y01}${A}.[ab]
#    /bin/rm -f ./pcip/*_${Y01}${A}.[ab]
#    /bin/rm -f ./ssta/*_${Y01}${A}.[ab]
#    /bin/rm -f ./ssto/*_${Y01}${A}.[ab]
#    /bin/rm -f ./cice/*_${Y01}${A}.[rB]
#  endif
C
  if (-e ./wind/${E}w${Y01}${A}.com) then
    /bin/mv ./wind/${E}w${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./wspd/${E}s${Y01}${A}.com) then
    /bin/mv ./wspd/${E}s${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./wvel/${E}v${Y01}${A}.com) then
    /bin/mv ./wvel/${E}v${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./grad/${E}g${Y01}${A}.com) then
    /bin/mv ./grad/${E}g${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./flux/${E}f${Y01}${A}.com) then
    /bin/mv ./flux/${E}f${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./flxt/${E}q${Y01}${A}.com) then
    /bin/mv ./flxt/${E}q${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./lrad/${E}l${Y01}${A}.com) then
    /bin/mv ./lrad/${E}l${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./mslp/${E}m${Y01}${A}.com) then
    /bin/mv ./mslp/${E}m${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./pcip/${E}p${Y01}${A}.com) then
    /bin/mv ./pcip/${E}p${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./ssta/${E}t${Y01}${A}.com) then
    /bin/mv ./ssta/${E}t${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./ssto/${E}o${Y01}${A}.com) then
    /bin/mv ./ssto/${E}o${Y01}${A}.{com,log} $D/..
  endif
C
C --- Wait for wind and flux interpolation of next segment.
C
  wait
C
  if (-e ./wind/${E}w${YXX}${B}.com) then
    /bin/mv ./wind/${E}w${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./wspd/${E}s${Y01}${A}.com) then
    /bin/mv ./wspd/${E}s${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./wvel/${E}v${Y01}${A}.com) then
    /bin/mv ./wvel/${E}v${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./grad/${E}g${Y01}${A}.com) then
    /bin/mv ./grad/${E}g${Y01}${A}.{com,log} $D/..
  endif
  if (-e ./flux/${E}f${YXX}${B}.com) then
    /bin/mv ./flux/${E}f${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./flxt/${E}q${YXX}${B}.com) then
    /bin/mv ./flxt/${E}q${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./lrad/${E}l${YXX}${B}.com) then
    /bin/mv ./lrad/${E}l${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./mslp/${E}m${YXX}${B}.com) then
    /bin/mv ./mslp/${E}m${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./pcip/${E}p${YXX}${B}.com) then
    /bin/mv ./pcip/${E}p${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./ssta/${E}t${YXX}${B}.com) then
    /bin/mv ./ssta/${E}t${YXX}${B}.{com,log} $D/..
  endif
  if (-e ./ssto/${E}o${YXX}${B}.com) then
    /bin/mv ./ssto/${E}o${YXX}${B}.{com,log} $D/..
  endif
endif
C
C --- wait for nesting .tar file.
C
if (-e ./nest) then
  wait
endif
C
C --- HYCOM error stop is implied by the absence of a normal stop.
C
touch summary_out
tail -1 summary_out
tail -1 summary_out | grep -c "^normal stop"
if ( `tail -1 summary_out | grep -c "^normal stop"` == 0 ) then
  cd $D/..
  /bin/mv LIST LIST_BADRUN
  echo "BADRUN" > LIST
endif
C
C --- wait for tar bundles to complete
C
wait
C
C  --- END OF MODEL RUN SCRIPT
C
