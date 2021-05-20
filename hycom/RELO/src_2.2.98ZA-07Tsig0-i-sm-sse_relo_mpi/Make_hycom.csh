#!/bin/csh
#
# --- Usage:  ./Make_hycom.com >& Make_hycom.log
#
# --- set ARCH to the correct value for this machine.
#

echo $ICEFLG
setenv ARCH $1
setenv ICEFLG $2

if     ($ICEFLG != 0) then
    setenv TYPE cice
    echo "Make_cice.csh: Environment variable TYPE=$TYPE"
    echo "Make_cice.csh: Environment variable ARCH=$ARCH"
    if     ($TYPE != "cice") then
        echo "TYPE must be cice to invoke cice make target"
        exit 1
    endif
#
else
    setenv TYPE hycom
    echo "Make_hycom.csh: Environment variable TYPE=$TYPE"
    echo "Make_hycom.csh: Environment variable ARCH=$ARCH"
endif
if (! -e ../config/${ARCH}_${TYPE}) then
    echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
    exit 1
endif


touch ./hycom_feature_flags
#
# --- make HYCOM component, and update hycom_cice
#
# --- force a relink, because CICE is not in the dependencies
if ( $ICEFLG == 0 ) then
    /bin/rm hycom
    echo "only hycom"
    make ARCH=$ARCH TYPE=$TYPE hycom
    mv -f hycom hycom_oasis
else
    echo "HYCOM-CICE"
    /bin/rm hycom_cice
    make ARCH=$ARCH TYPE=$TYPE CICE_DIR=./CICE/ hycom_cice
endif
# --- some machines require gmake
exit $status

