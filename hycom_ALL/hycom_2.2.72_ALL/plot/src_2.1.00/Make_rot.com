#
#set echo
#
# --- Usage:  ./Make_rot.com >& Make_rot.log
#
# --- make rot plot executables
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
# --- standard plot programs
#
foreach m ( hycomproc_rot fieldproc_rot )
  make -f Makefile_rot ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end
