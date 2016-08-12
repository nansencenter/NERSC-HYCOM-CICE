#
# --- delete old HYCOM related executables.
#
set echo
#
/bin/rm -f *.o
/bin/rm -f *.a
/bin/rm -f *.mod
#
foreach f ( echo2 endian )
  touch       ${f}
  /bin/rm -f  ${f}
end
foreach f ( cice_restart cice_stat hycom_palette lonlat_dist hycom_alat hycom_archm_dates hycom_archv_dates hycom_depth hycom_depth_40 hycom_nest_dates hycom_profile+sig hycom_profile+thstar hycom_profile2pcm hycom_profile2z hycom_profile2zi hycom_profile_argo hycom_profile_hybgen+ hycom_profile_hybgen hycom_profile_locsig hycom_profile_mld hycom_profile_remap hycom_sigma hycom_ts hycom_wind_date hycom_wind_ymdh hycom_ymdh_wind hycom_yoflat sigma0_to_sigma2 sigma2_to_sigma0 ts_to_sigma z2zi zi2z hycom_date_wind hycom_profile2plm hycom_profile_hybgen_34 hycom_profile_hybgen_35 hycom_subset_xy hycom_dp0k hycom_dp0k_cm hycom_dp0k_sigma )
  touch       ${f}
  /bin/rm -f  ${f}
end
foreach f ( hycom_crosscorr hycom_crosscorr_lag hycom_join unf42hycom unf82hycom hycom2raw hycom2raw8 hycom_1st_isopyc hycom_arctic hycom_arctic_ok hycom_bandmask hycom_binning hycom_binning_fld hycom_bouflx hycom_clip hycom_count hycom_eddy_center hycom_expr hycom_extract hycom_fill hycom_halfsm hycom_histogram hycom_ij2lonlat hycom_islands hycom_larger hycom_lonlat2ij hycom_lonlat2xy hycom_mask hycom_mass hycom_mean hycom_meanfit hycom_median hycom_meridional hycom_meridional_lon hycom_mixlay hycom_mixlay_old hycom_mxthrd hycom_NaN hycom_print hycom_range hycom_range_ij hycom_rivers hycom_rotate hycom_runmean hycom_sample hycom_sample_list hycom_sea_ok hycom_shift hycom_skill hycom_slopefit hycom_smooth hycom_stericssh hycom_subset hycom_superset hycom_thirdsm hycom_tidelat hycom_triple hycom_void hycom_xy2lonlat hycom_zonal hycom_zonal_lat ascii2hycom raw2hycom raw82hycom hycom_2d_ok hycom_autocorr hycom_autocorr_lag hycom_boxmean hycom_boxtime hycom_index_sort hycom_mask_ok hycom_mass_corr hycom_newzi hycom_quadlsq hycom_regression hycom_sstice hycom_profile_list hycom_botfric hycom_boxsmooth hycom_diflat hycom_merge hycom_sample_xy hycom_scatter hycom_tidebody hycom_vmean hycom_xward )
  touch       ${f}
  /bin/rm -f  ${f}
end
foreach f ( hycom_binning_nc hycom_force2nc hycom_profile2z_nc hycom_profile2s_nc )
  touch       ${f}
  /bin/rm -f  ${f}
end
#
foreach OS ( SunOS64 SunOS LinuxA OSF1 Linux Linux64 LinuxIFC LinuxICE LinuxGF XT5 IRIX64 AIX unicosmk )
  foreach f ( echo2 endian )
    if ( -e ${f}_${OS} ) then
      /bin/rm -f  ${f}_${OS}
    endif
  end
  foreach f ( clim_stat wind_stat wind_stat_check wind_stat_range wind_stat_range2 wind_stat_range5 wind_stat_raw hycom_sigma wind_stat_nc wind_stat_range_nc )
    if ( -e ${f}_${OS} ) then
      /bin/rm -f  ${f}_${OS}
    endif
  end
  foreach f ( cice_restart cice_stat hycom_palette lonlat_dist hycom_alat hycom_archm_dates hycom_archv_dates hycom_depth hycom_depth_40 hycom_nest_dates hycom_profile+sig hycom_profile+thstar hycom_profile2pcm hycom_profile2z hycom_profile2zi hycom_profile_argo hycom_profile_hybgen+ hycom_profile_hybgen hycom_profile_locsig hycom_profile_mld hycom_profile_remap hycom_sigma hycom_ts hycom_wind_date hycom_wind_ymdh hycom_ymdh_wind hycom_yoflat sigma0_to_sigma2 sigma2_to_sigma0 ts_to_sigma z2zi zi2z hycom_date_wind hycom_profile2plm hycom_profile_hybgen_34 hycom_profile_hybgen_35 hycom_subset_xy hycom_dp0k hycom_dp0k_cm hycom_dp0k_sigma )
    if ( -e ${f}_${OS} ) then
      /bin/rm -f  ${f}_${OS}
    endif
  end
  foreach f ( hycom_crosscorr hycom_crosscorr_lag hycom_join unf42hycom unf82hycom hycom2raw hycom2raw8 hycom_1st_isopyc hycom_arctic hycom_arctic_ok hycom_bandmask hycom_binning hycom_binning_fld hycom_bouflx hycom_clip hycom_count hycom_eddy_center hycom_expr hycom_extract hycom_fill hycom_halfsm hycom_histogram hycom_ij2lonlat hycom_islands hycom_larger hycom_lonlat2ij hycom_lonlat2xy hycom_mask hycom_mass hycom_mean hycom_meanfit hycom_median hycom_meridional hycom_meridional_lon hycom_mixlay hycom_mixlay_old hycom_mxthrd hycom_NaN hycom_print hycom_range hycom_range_ij hycom_rivers hycom_rotate hycom_runmean hycom_sample hycom_sample_list hycom_sea_ok hycom_shift hycom_skill hycom_slopefit hycom_smooth hycom_stericssh hycom_subset hycom_superset hycom_thirdsm hycom_tidelat hycom_triple hycom_void hycom_xy2lonlat hycom_zonal hycom_zonal_lat ascii2hycom raw2hycom raw82hycom hycom_2d_ok hycom_autocorr hycom_autocorr_lag hycom_boxmean hycom_boxtime hycom_index_sort hycom_mask_ok hycom_mass_corr hycom_newzi hycom_quadlsq hycom_regression hycom_sstice hycom_profile_list hycom_botfric hycom_boxsmooth hycom_diflat hycom_merge hycom_sample_xy hycom_scatter hycom_tidebody hycom_vmean hycom_xward )
    if ( -e ${f}_${OS} ) then
      /bin/rm -f  ${f}_${OS}
    endif
  end
  foreach f ( hycom_binning_nc hycom_scrip_nc hycom_force2nc hycom_profile2z_nc hycom_profile2s_nc hycom_seaice_nc )
    if ( -e ${f}_${OS} ) then
      /bin/rm -f  ${f}_${OS}
    endif
  end
end
