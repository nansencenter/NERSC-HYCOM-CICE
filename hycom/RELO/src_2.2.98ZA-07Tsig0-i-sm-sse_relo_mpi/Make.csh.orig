#!/bin/csh
#
set echo
cd $cwd
#
# --- Usage:  ./Make.com >& Make.log
#
# --- make hycom with TYPE from this directory's name (src_*_$TYPE).
# --- assumes dimensions.h is correct for $TYPE.
#
# --- set ARCH to the correct value for this machine.
# --- ARCH that start with A are for ARCTIC patch regions
#
#module swap compiler compiler/intel/12.1.3
#module swap mpi      mpi/intel/ibmpe
#module list
#setenv ARCH intelsse-pe-sm-relo
setenv ARCH xt4
setenv TYPE esmf

#
#module swap compiler compiler/intel/12.1.3
#module swap mpi      mpi/intel/impi/4.1.3
#module list
#setenv ARCH intelsse-impi-sm-relo
#setenv ARCH intelsse-impi-sm-SD-relo
#
#module switch PrgEnv-cray PrgEnv-intel
#module list
#setenv ARCH xc30-intel-relo
#
#setenv ARCH intel-pgi-SD-relo
#
setenv TYPE `echo $cwd | awk -F"_" '{print $NF}'`
echo "ARCH = " $ARCH "  TYPE = " $TYPE
#
if (! -e ../config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
endif
#
# --- some machines require gmake
#
#gmake ARCH=$ARCH TYPE=$TYPE hycom
make hycom ARCH=$ARCH TYPE=$TYPE
make esmf ARCH=$ARCH TYPE=$TYPE
#
if ( $ARCH == "Asp5" || $ARCH == "sp5") then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
endif
if ( $ARCH == "Asp6" || $ARCH == "sp6") then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
endif
if ( $ARCH == "Asp6-nofl" || $ARCH == "sp6-nofl") then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
endif
