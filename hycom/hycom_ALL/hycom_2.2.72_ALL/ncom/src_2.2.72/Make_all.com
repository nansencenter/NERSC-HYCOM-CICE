#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_all.com >& Make_all.log
#
# --- make all ncom executables
#
# --- set ARCH to the correct value for this machine.
#
source Make_all.src
#
printenv ARCH
#
if (! -e ../../config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
# --- ncom programs
#
#foreach m ( ncomc2archv ncom2archv grid2ncom )
foreach m ( ncomc2archv grid2ncom )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end
