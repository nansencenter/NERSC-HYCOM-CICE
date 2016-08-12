#
#set echo
#
# --- use 64K pages on POWER5 and POWER6
#
foreach f ( ap ap_nc aphf_add aphf_climo aphf_diurnal aphf_extend aphf_flcorr aphf_meanfit aphf_monthly aphf_offset aphf_scale aphf_tacorr force2nc kp kp_const kp_nc kphf_const kphf_table nrl2nc off_diff off_zero pcip_riv_hf pcip_riv_mon pcip_zero time_interp w_const wi wi_curl wi_magstress wi_meanfit wi_nc )
  ldedit -bdatapsize=64K -bstackpsize=64K $f
end 
