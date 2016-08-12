#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_all.com >& Make_all.log
#
# --- make all archive executables (except netCDF)
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
# --- archive modifying programs
#
foreach m ( hycomarchv micomarchv conv_archv hybgen_archv mrgl_archv ncoda_archv ncoda_archv_inc remap_archv remapi_archv trim_archv archt2archv archv2data2d archm2data2d archv2data2t archv2data3z archm2data3z archv2datasf archv2datasfl archv2datasfz archv2ncombc archv2restart archm2restart field2data restart2archv )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
  if (-e /usr/bin/ldedit) then
#   try to set medium pages on POWER5+
    /usr/bin/ldedit -bdatapsize=64K -bstackpsize=64K ${m}
  endif
end
