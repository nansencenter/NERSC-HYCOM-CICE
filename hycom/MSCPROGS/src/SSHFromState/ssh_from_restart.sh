#!/bin/bash
# Program

#prog=/home/fimm/nersc/knutali/Progs/TOOLS/SSHFromState/ssh_from_restart 
# This assumes that the work routine is in the same dir as the wrapper script
prog="${0%.sh}"

host=$(hostname | sed "s/\..*//")
xthost=$(xthostname)
if [ "$host" == "fimm" ] ; then
   export LD_LIBRARY_PATH=/local/Matlab-2007a/bin/glnxa64/:$LD_LIBRARY_PATH
#elif [ "$xthost" == "hexagon" ] ; then
#   #
#   #No matlab dir (yet)
#   echo "hexagon"
fi

$prog $@
