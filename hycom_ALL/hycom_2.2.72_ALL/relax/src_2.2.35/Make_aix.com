#
#set echo
#
# --- use 64K pages on POWER5 and POWER6
#
foreach f ( iso_density relax_flat_rivers relax_tracer relaxi relaxv rmu rmu2 tracer_const z_archive z_const z_levitus z_medatlas z_modify sst_pf sst_pf_4km ssh_modas sst_modas sst_gdem3 sst_gdem4 sst_woa z_gdem3 z_gdem4 z_woa_tracer )
  ldedit -bdatapsize=64K -bstackpsize=64K $f
end
