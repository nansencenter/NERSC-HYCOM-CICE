#
# --- create HYCOM related executables.
#
set echo
#
#setenv OS `/bin/uname`
setenv OS `uname`
if ($OS == "Linux") then
  if (`/bin/uname -m` == "x86_64") then
	setenv OS Linux64
  endif
# setenv OS LinuxIFC
# setenv OS LinuxICE
# setenv OS LinuxGF
# setenv OS XT5
endif
#if ($OS == "SunOS") then
#  setenv OS SunOS64
#endif
if ($OS == "UNICOS/mp") then
  setenv OS X1
endif
#
# --- the following are extracted from hycom/ALL/config/*_setup
#
switch ($OS)
case 'Linux64':
#       compile for 64-bit AMD64
	setenv FC	"pgf90"
	setenv FFLAGS	"-g -fast -byteswapio -tp k8-64 -mcmodel=medium -Mnolarge_arrays"
	setenv FLIBS	""
	setenv CC	"gcc"
	setenv CFLAGS	"-O -march=k8 -m64 -mcmodel=medium"
	breaksw
case 'LinuxICE':
#       compile for SGI Altix ICE, Intel compiler
	setenv FC	"ifort"
	setenv FFLAGS	"-g -O3 -fp-model source -convert big_endian"
	setenv FLIBS	"-shared-intel"
	setenv CC	"icc"
	setenv CFLAGS	"-O"
	breaksw
case 'LinuxGF':
#       compile for gfortran
	setenv FC	"gfortran"
	setenv FFLAGS	"-fPIC -fno-second-underscore -fconvert=big-endian -O"
	setenv FLIBS	""
	setenv CC	"gcc"
	setenv CFLAGS	"-fPIC -fno-second-underscore -O"
	breaksw
case 'XT5':
#       compile for XT5 via aprun
        setenv FC       "env XTPE_INFO_MESSAGE_OFF=1 ftn"
        setenv FFLAGS   "-g -fast -byteswapio -tp barcelona-64 -mcmodel=medium -Mnolarge_arrays"
        setenv FLIBS    ""
        setenv CC       "env XTPE_INFO_MESSAGE_OFF=1 cc"
        setenv CFLAGS   "-g -O -tp barcelona-64 -mcmodel=medium -Mnolarge_arrays"
        breaksw
case 'LinuxIFC':
#       compile for Pentium 4, Intel compiler
	setenv FC	"ifort"
	setenv FFLAGS	"-g -tpp7 -O3 -convert big_endian"
	setenv FLIBS	"-Vaxlib"
	setenv CC	"gcc"
	setenv CFLAGS	"-O -march=pentium4 -m32"
	breaksw
case 'Linux':
#       compile for Pentium 4 (also 32-bit AMD64)
	setenv FC	"pgf90"
	setenv FFLAGS	"-g -fast -byteswapio -tp p7"
	setenv FLIBS	"-Mlfs"
	setenv CC	"gcc"
	setenv CFLAGS	"-O -march=pentium4 -m32"
	breaksw
case 'AIX':
	setenv FC	"xlf95"
	setenv FFLAGS	"-qfixed -O3 -qstrict -qarch=pwr3 -qtune=pwr3 -qcache=auto -qspillsize=32000 -q64"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-O -q64 -DAIX"
	breaksw
case 'IRIX64':
	setenv FC	"f90"
	setenv FFLAGS	"-g3 -64 -O3 -macro_expand"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"g3 -64 -O3"
	breaksw
case 'SunOS64':
	setenv FC	"f95"
	setenv FFLAGS	"-g -fast -xarch=native64 -nodepend -xvector=no -O2 -xpp=cpp"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-g -fast -xarch=native64"
	breaksw
case 'SunOS':
	setenv FC	"f95"
	setenv FFLAGS	"-g -fast -nodepend -xvector=no -O2 -xpp=cpp"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-g -fast"
	breaksw
case 'OSF1':
	setenv FC	"f90"
	setenv FFLAGS	"-g3 -fpe1 -fast -O5 -convert big_endian -assume byterecl"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	"-g3 -fast"
	breaksw
case 'unicosmk':
	setenv FC	"f90"
	setenv FFLAGS	"-X 1 -V -f fixed -O scalar2,unroll2,pipeline1,vector3 -d p -M 801"
	setenv FLIBS	""
	setenv CC	"cc"
	setenv CFLAGS	""
	breaksw
case 'X1':
	setenv FC	"ftn"
	setenv FFLAGS	"-Ossp"
	setenv FLIBS	"x1_sys.o"
	setenv CC	"cc"
	setenv CFLAGS	"-hssp -UCRAY"
	$FC $FFLAGS -c x1_sys.f
	breaksw
default:
	echo 'Unknown Operating System: ' $OS
	exit (1)
endsw
#
# --- *.c programs
#
foreach f ( echo2 endian )
  if ( ! -e ${f}_${OS} ) then
    $CC $CFLAGS ${f}.c -o ${f}_${OS}
  else if ( -f `find ${f}.c -prune -newer ${f}_${OS}` ) then
    $CC $CFLAGS ${f}.c -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
# --- *.f programs
#
foreach f ( clim_stat wind_stat wind_stat_check wind_stat_range wind_stat_range2 wind_stat_range5 wind_stat_raw hycom_sigma )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
end
foreach f ( cice_restart cice_stat hycom_palette lonlat_dist hycom_alat hycom_archm_dates hycom_archv_dates hycom_depth hycom_depth_40 hycom_nest_dates hycom_profile+sig hycom_profile+thstar hycom_profile2pcm hycom_profile2z hycom_profile2zi hycom_profile_argo hycom_profile_hybgen+ hycom_profile_hybgen hycom_profile_locsig hycom_profile_mld hycom_profile_remap hycom_sigma hycom_ts hycom_wind_date hycom_wind_ymdh hycom_ymdh_wind hycom_yoflat sigma0_to_sigma2 sigma2_to_sigma0 ts_to_sigma z2zi zi2z hycom_date_wind hycom_profile2plm hycom_profile_hybgen_34 hycom_profile_hybgen_35 hycom_subset_xy hycom_dp0k hycom_dp0k_cm hycom_dp0k_sigma )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f $FLIBS -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
# --- *.F programs, may need hycom_endian_io.o and/or parse.o.
#
$FC $FFLAGS -c hycom_endian_io.F
$CC $CFLAGS -c parse.c
#
foreach f ( hycom_crosscorr hycom_crosscorr_lag hycom_join unf42hycom unf82hycom hycom2raw hycom2raw8 hycom_1st_isopyc hycom_arctic hycom_arctic_ok hycom_bandmask hycom_binning hycom_binning_fld hycom_bouflx hycom_clip hycom_count hycom_eddy_center hycom_expr hycom_extract hycom_fill hycom_halfsm hycom_histogram hycom_ij2lonlat hycom_islands hycom_larger hycom_lonlat2ij hycom_lonlat2xy hycom_mask hycom_mass hycom_mean hycom_meanfit hycom_median hycom_meridional hycom_meridional_lon hycom_mixlay hycom_mixlay_old hycom_mxthrd hycom_NaN hycom_print hycom_range hycom_range_ij hycom_rivers hycom_rotate hycom_runmean hycom_sample hycom_sample_list hycom_sea_ok hycom_shift hycom_skill hycom_slopefit hycom_smooth hycom_stericssh hycom_subset hycom_superset hycom_thirdsm hycom_tidelat hycom_triple hycom_void hycom_xy2lonlat hycom_zonal hycom_zonal_lat ascii2hycom raw2hycom raw82hycom hycom_2d_ok hycom_autocorr hycom_autocorr_lag hycom_boxmean hycom_boxtime hycom_index_sort hycom_mask_ok hycom_mass_corr hycom_newzi hycom_quadlsq hycom_regression hycom_sstice hycom_botfric hycom_boxsmooth hycom_diflat hycom_merge hycom_sample_xy hycom_scatter hycom_tidebody hycom_vmean hycom_xward )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_endian_io.o parse.o -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_endian_io.o parse.o -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
# --- archive reading programs, may additionally need hycom_profile_lib.o
#
$FC $FFLAGS -c hycom_profile_lib.F
#
foreach f ( hycom_profile_list )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F $FLIBS hycom_profile_lib.o hycom_endian_io.o parse.o -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
