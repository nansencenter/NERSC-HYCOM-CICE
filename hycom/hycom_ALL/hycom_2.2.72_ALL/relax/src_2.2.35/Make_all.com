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
#foreach m ( bottom iso_density relax_flat_rivers relax_tracer relaxi relaxv rmu rmu2 tracer_const z_archive z_const z_levitus z_medatlas z_modify sst_modas sst_rs z_modas )
foreach m ( iso_density relax_flat_rivers relax_tracer relaxi relaxi_nemo  relaxv rmu rmu2 tracer_const z_archive z_const z_levitus z_medatlas z_modify sst_pf sst_pf_4km z_woa_dynamic relaxi_archvz nemo_archvz  nemo_archvz_modified )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end
