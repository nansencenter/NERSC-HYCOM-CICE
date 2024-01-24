#!/bin/csh
#
# --- Usage:  ./Make_cice.com >& Make_cice.log
#
# --- make cice (ESMF HYCOM component) with TYPE=cice.
# --- this directory's name must be src_*_cice.
# --- assumes dimensions.h is correct for TYPE=cice (i.e. for mpi).
#
# --- set ARCH to the correct value for this machine.
# --- ARCH that start with A are for ARCTIC patch regions
#
#setenv ARCH alphaL
#setenv ARCH alpha
#setenv ARCH amd64
#setenv ARCH intel
#setenv ARCH o2k
#setenv ARCH sp3
#setenv ARCH sp4
#setenv ARCH sun64
#setenv ARCH sun
#setenv ARCH t3e
#setenv ARCH xt3
#
#setenv ARCH xt4
#
#KAL - ARCH must now be set before running Make_cice.csh script. It is of type
#KAL - xt4.hexagon ...
setenv TYPE cice
#setenv ARCH Linux.sisu.intel   
setenv CICE_FLAG $2
setenv ARCH $1
echo "(1) ARCH = $ARCH"
echo "if error: make sure to set the correct ARCH in Make_cice.csh"
echo "Make_cice.csh: Environment variable TYPE=$TYPE"
echo "Make_cice.csh: Environment variable ARCH=$ARCH"



if     ($TYPE != "cice") then
  echo "TYPE must be cice to invoke cice make target"
  exit 1
endif
#
if (! -e ../config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
endif
#
# --- cice needs additional environment variables.
#
#if ($TYPE == "cice") then
#  switch ($ARCH)
#  case 'sp5':
#    setenv BEI_HOME /site/BEI
#    setenv ESMF_DIR ${BEI_HOME}/esmf/4.0.0rp2
#    breaksw
#  case 'o2k':
#    setenv BEI_HOME /usr/local/usp/BEI
#    setenv ESMF_DIR ${BEI_HOME}/esmf/4.0.0rp2
#    breaksw
#  case 'xt3':
#    setenv BEI_HOME /usr/local/usp/BEI
#    setenv ESMF_DIR ${BEI_HOME}/esmf/4.0.0rp2
#    breaksw
#  default:
#    echo "TYPE = cice  needs BEI_HOME and ESMF_DIR"
#    exit (1)
#  endsw
#endif
#setenv ESMF_DIR /home/nersc/knutali/opt/esmf_5_2_0rp3-nonetcdf/

# --- KAL. Touch this file to make sure it exists. It may be empty, but the makefile will look for it
touch ./hycom_feature_flags
# setup cpp flags
setenv NERSC_FLAG "-DNERSC_HYCOM_CICE -DNERSC_USE_ESMF -DNERSC_ATM_CPL -DNERSC_saltflux -DNERSC_T2F"
#
# --- make HYCOM component, and update hycom_cice
#
# --- force a relink, because CICE is not in the dependencies
/bin/rm hycom_cice
if ($CICE_FLAG == 0) then
	echo "only hycom"
      	make ARCH=$ARCH TYPE=$TYPE CICE_FLAG=$CICE_FLAG hycom
      	echo "Replacing HYCOM with HYCOM_CICE ..."
      	mv hycom hycom_cice
else
	make ARCH=$ARCH TYPE=$TYPE CICE_DIR=../CICE/ hycom_cice_nersc
endif
# --- some machines require gmake
#gmake ARCH=$ARCH TYPE=$TYPE hycom_cice
exit $status
