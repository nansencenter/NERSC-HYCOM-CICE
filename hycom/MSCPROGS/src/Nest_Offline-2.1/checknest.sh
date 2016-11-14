#!/bin/bash
# script calls checknest.sh with correctly set libraries


prog=$( echo $0 | sed "s/\.sh$//")
echo $prog
export LD_LIBRARY_PATH=/local/Matlab-2007a/bin/glnxa64/:$LD_LIBRARY_PATH

$prog $@

