#!/bin/csh
#
#set echo
#
# --- Usage:  ./Make_ncdf.com >& Make_ncdf.log
#
# --- make all netCDF force executables
#
# --- set NCDF to the root directory for netCDF version 3.5.
# --- available from: http://www.unidata.ucar.edu/packages/netcdf/
#
source ../../Make_ncdf.src
#
# --- set ARCH to the correct value for this machine.
#
source Make_all.src
#
echo "NCDF = " $NCDF
echo "ARCH = " $ARCH
#
if (! -e ../../config/${ARCH}_setup) then
  echo "ARCH = " $ARCH "  is not supported"
  exit 1
endif
#
# --- softlink to netCDF module and library (and typesizes.mod for OSF1 only)
#
/bin/rm -f netcdf.mod libnetcdf.a
/bin/rm -f netcdf.inc
/bin/rm -f typesizes.mod
#
ln -s ${NCDF}/include/*.mod   .
ln -s ${NCDF}/include/*.inc   .
ln -s ${NCDF}/lib/libnetcdf.a .
#
# --- netCDF force programs
#
foreach m ( force2nc nrl2nc ap_nc kp_nc tp_nc wi_nc )
  make ${m} ARCH=${ARCH} >&! Make_${m}
  if ($status) then
    echo "Make failed:" ${m}
  else
    echo "Make worked:" ${m}
  endif
end
