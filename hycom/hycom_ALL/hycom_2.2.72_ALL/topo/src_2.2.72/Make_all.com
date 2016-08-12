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
foreach m ( malt2h m2h p2h bathy_02min bathy_05min hudson landsea_02min landsea_05min 1d 2d batrop clip diff edit flat grid_360 hgram landfill landmask latitude lonlat lonlat_2d lpanam map mapsub mercator merge modify onesea onesea-b onesea_fill panam partit partit_noshrink partit_arctic partit_arctic_ns ppmX ports rotated rough shrink slope smallsea smooth smooth_skip subset tiles zcells zthin )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end
