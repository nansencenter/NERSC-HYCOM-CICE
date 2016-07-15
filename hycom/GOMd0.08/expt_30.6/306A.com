#!/bin/csh
#$ -S /bin/csh
#$ -P ONRDC10855C79
#$ -cwd
#$ -M metzger@nrlssc.navy.mil
#$ -m abe
#$ -j y
#$ -l ibmp3
#$ -l 4hr
#
set echo
set time = 1
set timestamp
C
C --- tar archive files in a batch job
C
setenv OS `uname`
C
C --- pget, pput "copy" files between scratch and permanent storage.
C --- Can both be cp if the permanent filesystem is mounted locally.
C
switch ($OS)
case 'SunOS':
#case 'Linux':
#case 'IRIX64':
case 'OSF1':
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
C --- E is expt number.
C --- P is primary path.
C --- D is permanent directory.
C --- S is scratch   directory, must not be the permanent directory.
C
setenv E 306
setenv P hycom/GOMd0.08/expt_30.6/data
setenv D ~/$P
C
switch ($OS)
case 'HP-UX':
    setenv S     /scratch3/${user}/$P
    breaksw
case 'OSF1':
#                ASC MSRC
    setenv S     /workspace/${user}/$P
    breaksw
case 'SunOS':
    if (-e /net/hermes/scrb) then
#                  NRLSSC
      setenv S     /net/hermes/scrb/${user}/$P
    else
#                  NAVO MSRC
      setenv S     /scr/${user}/$P
    endif
    breaksw
case 'IRIX64':
    setenv S     /scr/${user}/$P
    breaksw
case 'AIX':
    if      (-e /gpfs/work) then
#                  ERDC MSRC, under PBS
      setenv S     /gpfs/work/${user}/$P
      setenv POE  pbspoe
    else if (-e /scr) then
#                  NAVO MSRC, under LoadLeveler
      setenv S     /scr/${user}/$P
      setenv POE  poe
    else if (-e /scratch/tempest) then
#                  MHPCC DC, under LoadLeveler
      setenv S     /scratch/tempest/users/${user}/$P
      setenv POE  poe
    else
#                  ARL MSRC, under GRD
      setenv S     /usr/var/tmp/${user}/$P
      setenv POE  grd_poe
    endif
    breaksw
case 'unicos':
case 'unicosmk':
case 'sn6705':
    setenv S     /work/${user}/$P
    breaksw
endsw
C
setenv A " "
setenv Y01 "001"
C
C --- run in the tar directory.
C
chmod 750 $S/tar_${Y01}${A}
cd        $S/tar_${Y01}${A}
C
/bin/rm  -f  ../${E}_archv_${Y01}${A}.tar ${E}_archv.dummy.*
/bin/tar cvf ../${E}_archv_${Y01}${A}.tar .
C
C --- clean up in the primary scratch directory.
C
cd ..
/bin/rm -r tar_${Y01}${A}
${pput} ${E}_archv_${Y01}${A}.tar $D/${E}_archv_${Y01}${A}.tar
