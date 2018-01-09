#
#set echo
#
# --- Usage:  ./Make_all.com >& Make_all.log
#
# --- make all relax setup executables, except those needing netCDF.
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
#foreach m ( iso_density relax_flat relax_tracer rmu rmu2 tracer_const z_archive z_const z_levitus z_modify sst_modas sst_rs z_modas )
foreach m ( iso_density relax_flat relax_tracer rmu rmu2 tracer_const z_archive z_const z_levitus z_modify )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end
#
#foreach s ( 0 2 4 )
foreach s ( 0 2 )
  foreach m ( relax_sig${s} relaxv_sig${s} relax_zon${s} bottom_sig${s} mxlay_sig${s} )
    if (-e ${m}.F) then
      make ${m} ARCH=${ARCH} >&! Make_${m}
      if ($status) then
        echo "Make failed:" ${m}
      else
        echo "Make worked:" ${m}
      endif
    endif
  end
end
