#
#set echo
#
# --- Usage:  ./Make_all.com >& Make_all.log
#
# --- make all topo setup executables
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
foreach m ( sub_topog sub_topog_panam subregion isub_topog isubregion isuba_arche isuba_field isuba_gmapi isuba_topog isubaregion half_topog sub_grid )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
  if (-e /usr/bin/ldedit) then
#   try to set medium pages on POWER5+ and POWER6
    /usr/bin/ldedit -bdatapsize=64K -bstackpsize=64K ${m}
  endif                                                  
end
