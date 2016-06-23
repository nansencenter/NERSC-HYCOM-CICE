#!/bin/csh
set echo
#
# --- build new expt files from old.
# --- some files will need additional manual editing.
#
# --- customized version for the GOMd0.08 2.? series of expts.
#
# DO = old experiment directory name
#  O = old experiment number
# DN = new experiment directory name
#  N = new experiment number
#  R = region name.
#
setenv DO expt_30.5
setenv  O 305
setenv DN expt_30.6
setenv  N 306
setenv  R GOMd0.08
#
setenv RO $R
setenv  D ../../${RO}/${DO}
#
foreach t ( .com A.com F.com G.com M.com O.com P.com S.com T.com V.com W.com y099l.limits pbs.com )
  sed -e "s/setenv E .*${O}.*/setenv E ${N}/g" -e "s/${DO}/${DN}/g"  -e "s/${RO}/${R}/g" ${D}/${O}${t} >! ${N}${t}
end
#
foreach f ( cp_all.com )
  sed -e "s/setenv E .*${O}.*/setenv E ${N}/g" -e "s/${DO}/${DN}/g"  -e "s/${RO}/${R}/g" ${D}/${f} >! ${f}
end
#
cp ${D}/${O}.awk ${N}.awk
#
cp ${D}/*.input* .
#
cp ${D}/dummy*.com .
#
# --- experiment run sequence:
#
mlist 099 099 1 i j k
cp LIST LIST++
#
# --- create data directory (local and on archive):
#
mkdir data
krsh newton.navo.hpc.mil mkdir -p /u/home/$user/hycom/${R}/${DN}/data
#
# --- nesting for 2.0, forcing from 30.0, everything else from expt 06.0:
#
krsh newton.navo.hpc.mil ln -s "../../expt_02.0/data/nest" /u/home/${user}/hycom/GOMd0.08/${DN}/data/nest
krsh newton.navo.hpc.mil ln "/u/home/${user}/hycom/GOMd0.08/expt_30.0/data/forcing_cfsr_099.tar" /u/home/${user}/hycom/GOMd0.08/${DN}/data
krsh newton.navo.hpc.mil ln -s "060" /u/home/${user}/hycom/GOMd0.08/relax/$N
krsh newton.navo.hpc.mil ln "/u/home/${user}/hycom/GOMd0.08/expt_06.0/data/restart_099i.?" /u/home/${user}/hycom/GOMd0.08/${DN}/data
