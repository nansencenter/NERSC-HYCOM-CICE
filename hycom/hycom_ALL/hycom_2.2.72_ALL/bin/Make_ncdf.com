#
# --- create HYCOM related netCDF executables.
# --- run after Make_all.com
#
set echo
#
# --- set NCDF to the root directory for netCDF version 3.5 or later.
# --- available from: http://www.unidata.ucar.edu/packages/netcdf/
#
source ../Make_ncdf.src
#
# --- softlink to netCDF module and library (and typesizes.mod if needed)
#
/bin/rm -f netcdf.mod libnetcdf.a
/bin/rm -f typesizes.mod
#
#ln -s ${NCDF}/include/*.mod   .
#ln -s ${NCDF}/lib/libnetcdf.a .
#
setenv OS `/bin/uname`
if ($OS == "Linux") then
  if (`/bin/uname -m` == "alpha") then
	setenv OS LinuxA
  endif
  if (`/bin/uname -m` == "x86_64") then
	setenv OS Linux64
  endif
# setenv OS Linux_Fram
 setenv OS LinuxGF_NC
# setenv OS XT5
endif
#if ($OS == "SunOS") then
#  setenv OS SunOS64
#endif
#
# --- the following are extracted from hycom/ALL/config/*_setup
#
switch ($OS)
case 'Linux64':
#       compile for 64-bit AMD64, probably also ok for Intel EM64T (Nocona)
	setenv FC	"pgf90"
	setenv FFLAGS	"-g -fast -byteswapio -tp k8-64 -mcmodel=medium -Mnolarge_arrays"
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
	setenv CC	"gcc"
	setenv CFLAGS	"-O -march=pentium4 -m32"
	breaksw
case 'AIX':
	setenv FC	"xlf95"
	setenv FFLAGS	"-qfixed -O3 -qstrict -qarch=pwr3 -qtune=pwr3 -qcache=auto -qspillsize=32000 -q64"
	setenv CC	"cc"
	setenv CFLAGS	"-O -DAIX"
	breaksw
case 'IRIX64':
	setenv FC	"f90"
	setenv FFLAGS	"-g3 -64 -O3 -macro_expand"
	setenv CC	"cc"
	setenv CFLAGS	"g3 -64 -O3"
	breaksw
case 'SunOS64':
	setenv FC	"f95"
	setenv FFLAGS	"-g -fast -xarch=native64 -nodepend -xvector=no -O2 -xpp=cpp"
	setenv CC	"cc"
	setenv CFLAGS	"-g -fast -xarch=native64"
	breaksw
case 'SunOS':
	setenv FC	"f95"
	setenv FFLAGS	"-g -fast -nodepend -xvector=no -O2 -xpp=cpp"
	setenv CC	"cc"
	setenv CFLAGS	"-g -fast"
	breaksw
case 'OSF1':
	setenv FC	"f90"
	setenv FFLAGS	"-g3 -fpe1 -fast -O5 -convert big_endian -assume byterecl"
	setenv CC	"cc"
	setenv CFLAGS	"-g3 -fast"
	breaksw
case 'unicosmk':
	setenv FC	"f90"
	setenv FFLAGS	"-X 1 -V -f fixed -O scalar2,unroll2,pipeline1,vector3 -d p -M 801"
	setenv CC	"cc"
	setenv CFLAGS	""
	breaksw
case 'LinuxGF_NC':
#       compile for gfortran
	setenv FC	"gfortran"
	#setenv FFLAGS	"-I/cluster/software/easybuild/software/netCDF-Fortran/4.4.4-foss-2016b/include -fPIC -fno-second-underscore -fconvert=big-endian -O"
	setenv FFLAGS	"-I/cluster/software/GCCcore/6.3.0/include/c++/6.3.0 -fPIC -fno-second-underscore -fconvert=big-endian -O"
	setenv FLIBS	"-L/cluster/software/netCDF-Fortran/4.4.4-foss-2017a-HDF5-1.8.18/lib -lfftw3 -lnetcdff -lnetcdf "
	setenv CC	"gcc"
	setenv CFLAGS	"-fPIC -fno-second-underscore -O"
	breaksw
case 'Linux_Fram':
#       compile for gfortran 
	setenv FC	"gfortran"
	setenv FFLAGS	"-I/cluster/software/GCCcore/6.3.0/include/c++/6.3.0 -fPIC -fno-second-underscore -fconvert=big-endian -O"
	setenv FLIBS	"-lfftw3 -lnetcdff -lnetcdf "
	setenv CC	"icc"
	setenv CFLAGS	"-fPIC -fno-second-underscore -O"
	breaksw
default:
	echo 'Unknown Operating System: ' $OS
	exit (1)
endsw
#
foreach f ( hycom_binning_nc hycom_scrip_nc )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F hycom_endian_io.o parse.o ${FLIBS} ${EXTRANCDF} -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F hycom_endian_io.o parse.o ${FLIBS} ${EXTRANCDF} -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
foreach f ( hycom_profile2z_nc hycom_profile2s_nc hycom_seaice_nc )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.F hycom_profile_lib.o hycom_endian_io.o parse.o ${FLIBS} ${EXTRANCDF} -o ${f}_${OS}
  else if ( -f `find ${f}.F -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.F hycom_profile_lib.o hycom_endian_io.o parse.o ${FLIBS} ${EXTRANCDF} -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
  /bin/ln -s  ${f}_${OS} ${f}
end
#
foreach f ( wind_stat_nc wind_stat_range_nc )
  if ( ! -e ${f}_${OS} ) then
    $FC $FFLAGS ${f}.f ${FLIBS} ${EXTRANCDF} -o ${f}_${OS}
  else if ( -f `find ${f}.f -prune -newer ${f}_${OS}` ) then
    $FC $FFLAGS ${f}.f ${FLIBS} ${EXTRANCDF} -o ${f}_${OS}
  else
    echo "${f}_${OS} is already up to date"
  endif
  touch       ${f}
  /bin/rm -f  ${f}
  chmod a+rx  ${f}_${OS}
end
