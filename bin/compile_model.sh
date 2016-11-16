#!/bin/bash
# Script for quickly setting up model source code and compiling it. 


# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
export EDIR=$(pwd)
source ${BINDIR}/common_functions.sh || { echo "Could not source ${BINDIR}/common_functions.sh" ; exit 1 ; }
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ./EXPT.src || { echo "Could not source ./EXPT.src" ; exit 1 ; }
sourcedir=$BASEDIR/../hycom/RELO/src_${V}ZA-07Tsig0-i-sm-sse_relo_mpi/ # Yes, we always use this version
sourceconfdir=$BASEDIR/../hycom/RELO/config/



usage="
   This script must be run in an experiment directory as it needs to read EXPT.src

   On first invocation this script will fetch hycom source from $sourcedir,
   place it in the build directory, set up and compile the model.

   The script will also set up the equation of state to use in hycom, depending on 
   the setting of SIGVER in EXPT.src.

   NB: On subsequent calls, this script will compile the model, but not update with code
   from $sourcedir, unless you provide the -u flag

   Why choose this behaviour? 
   --------------------------
   This behaviour makes it possible to modify and test code in the build dir without
   affecting the code in other experiments, or the code in  $sourcedir .
   When you are satisfied with the code changes in build, the code can be brought into the main
   code in $sourcedir, then pushed to the code repository

   Example:
      $(basename $0) [ -u ]  [ -m mpi_library ]  compiler

   arguments   :
      compiler : compiler to use for. Currently supported: ifort gfortran pgi

   optional arguments :
      -u              : update code in build dir from NHCROOT/hycom/RELO/
      -m mpi_library  : on some machines you need to specify what mpi library to use


   Examples:
      $(basename $0) pgi
         will compile the model code inside the local build/ directory. Portland compilers are used.
         If the build dir does not exist,it will be created and code synced from $sourcedir. 

      $(basename $0) -u  pgi
         will compile the model code inside the local build/ directory. Portland compilers are used.
         Code will ALWAYS be synced from $sourcedir, so local changes will be overwritten.
"

options=$(getopt -o m:u -- "$@")
[ $? -eq 0 ] || {
    echo "$usage"
    echo "Incorrect options provided"
    exit 1
}
mpilib=""
eval set -- "$options"
while true; do
    case "$1" in
    -m)
       shift;
       mpilib=$1
        ;;
    -u)
       update="update"
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done


#Check arguments 
if [ $# -gt 0 ] ; then
   compiler=$1
else 
   echo "$usage"
   echo "Need to provide compiler to script, options: gfortran pgi ifort "
   exit 1
fi
echo "$(basename $0) : compiler=$compiler"
echo "$(basename $0) : mpilib=$mpilib"


# Check ARCH 
ARCH=$(uname -s)
ARCHLOWER=$(echo $ARCH | tr '[:upper:]' '[:lower:]' )
if [[ "$ARCH" == "Linux" ]] ;then
   true
else 
   echo "$usage"
   echo
   echo
   echo "Unsupported ARCH=$ARCH"
   exit 2
fi
echo "$(basename $0) : ARCH=$ARCH"

# SITE deduced from hostname. 
unames=$(uname -s)
# Hardcoded cases - hexagon
if [ "${unamen:0:7}" == "hexagon" ] ; then
   SITE="hexagon"
   MACROID=$ARCH.$SITE.$compiler

# Generic case. SITE is empty
elif [ "${unames:0:5}" == "Linux" ] ; then
   SITE=""
   if [ -z "${mpilib}" ] ; then
      echo "mpilib must be set on input running on generic linux machine (-m option)"
      exit 4
   fi
   MACROID=$ARCH.$compiler.$mpilib

else
   echo "Unknown SITE. uname -n gives $unamen"
   exit 3
fi
echo "$(basename $0) : SITE=$SITE"




# Deduce ESMF dir from SITE and possibly ARCH
if [[ -n "${ESMF_DIR}" ]] &&  [[ -n "${ESMF_MOD_DIR}" ]] && [[ -n "${ESMF_LIB_DIR}" ]] ; then
   echo "Using preset ESMF_DIR    =$ESMF_DIR"
   echo "Using preset ESMF_MOD_DIR=$ESMF_DIR"
   echo "Using preset ESMF_LIB_DIR=$ESMF_DIR"

# If site is given, use hardcoded settings for this machine
elif [ "$SITE" == "hexagon" ] ; then
   module load craype-interlagos
   if [[ -z "${ESMF_DIR}" ]] ; then
      # use as default
      export ESMF_DIR=/home/nersc/knutali/opt/esmf_5_2_0rp3-nonetcdf/
   fi
   export ESMF_MOD_DIR=${ESMF_DIR}/mod/modO/Unicos.$compiler.64.mpi.default/
   export ESMF_LIB_DIR=${ESMF_DIR}/lib/libO/Unicos.$compiler.64.mpi.default/


# If site is not given, try to use a generic setup. MAcro names composed of compiler name and mpi lib name (openmpi, mpich, lam, etc etc(
elif [ "${unames:0:5}" == "Linux" -a "$SITE" == "" ] ; then
   if [ -z "${ESMF_DIR}" ] ; then
      echo "ESMF_DIR must be set before calling script when running on generic linux machine"
      exit 4
   fi
   export ESMF_MOD_DIR=${ESMF_DIR}/mod/modO/Linux.$compiler.64.$mpilib.default/
   export ESMF_LIB_DIR=${ESMF_DIR}/lib/libO/Linux.$compiler.64.$mpilib.default/

else 
   echo "Dont know where ESMF_DIR is located on this machine"
   exit 4
fi
echo "$(basename $0) : ESMF_DIR=$ESMF_DIR"
echo "$(basename $0) : ESMF_MOD_DIR=$ESMF_MOD_DIR"
echo "$(basename $0) : ESMF_LIB_DIR=$ESMF_LIB_DIR"

# Get some useful info from blkdat.input
THFLAG=$(blkdat_get blkdat.input thflag)
IDM=$(blkdat_get blkdat.input idm)
JDM=$(blkdat_get blkdat.input jdm)

#Get statement function include file from SIGVER (in EXPT.src)
if [ $SIGVER -eq 1 ] ; then
   TERMS=7
   MYTHFLAG=0
elif [ $SIGVER -eq 2 ] ; then
   TERMS=7
   MYTHFLAG=2
elif [ $SIGVER -eq 3 ] ; then
   TERMS=9
   MYTHFLAG=0
elif [ $SIGVER -eq 4 ] ; then
   TERMS=9
   MYTHFLAG=2
elif [ $SIGVER -eq 5 ] ; then
   TERMS=17
   MYTHFLAG=0
elif [ $SIGVER -eq 6 ] ; then
   TERMS=17
   MYTHFLAG=2
elif [ $SIGVER -eq 7 ] ; then
   TERMS=12
   MYTHFLAG=0
elif [ $SIGVER -eq 8 ] ; then
   TERMS=12
   MYTHFLAG=2
else
   echo "SIGVER = $SIGVER"
   echo "only 7,9,12 and 17 term eq of state is supported (SIGVER=1 through 8) is supported"
   exit 1
fi
TERMS2=$(echo 0$TERMS | tail -c3)
echo "SIGVER      = $SIGVER .There are $TERMS terms in equation of state"

# Sanity check
if [ $THFLAG -ne $MYTHFLAG ] ; then
   echo "thflag in blkdat.input ($THFLAG) does not match thflag from eq of state ($MYTHFLAG)"
   exit 1
fi

# Copy code to expt dir
targetdir=$(source_dir $V $TERMS $THFLAG)
targetdir=$EDIR/build/$targetdir
targetconfdir=$EDIR/build/config/
if [ ! -d $EDIR/build/ ] ; then 
   mkdir $EDIR/build
   echo "build dir $EDIR/build not found. Setting it up with repo code from $sourcedir"
   rsync -avhL $sourcedir/ $targetdir/
   rsync -avhL $sourceconfdir/ $targetconfdir/
else 
   if [ "$update" == "update" ] ; then
      echo "build dir $EDIR/build found. Updating code in that subdirectory"
      rsync -avhL $sourcedir/ $targetdir/
      rsync -avhL $sourceconfdir/ $targetconfdir/
   else 
      echo "build dir $EDIR/build found. Using code in that subdirectory [not updated with code in $sourcedir]"
   fi
fi

# Copy hycom feature flag in expt dir if present
[ -f $EDIR/hycom_feature_flags  ] && cp $EDIR/hycom_feature_flags $targetdir

# Set up correct eq of state for hycom
stmt=stmt_fns_SIGMA${MYTHFLAG}_${TERMS}term.h
cd $targetdir
echo "Now setting up stmt_fns.h in $targetdir"
rm stmt_fns.h
ln -s ALT_CODE/$stmt stmt_fns.h

# 1) Compile CICE
cd $targetdir/CICE/
env RES=gx3 GRID=${IDM}x${JDM} SITE=$SITE MACROID=$MACROID ./comp_ice.esmf
res=$?
if [ $res -ne 0 ] ; then 
   echo
   echo "Error when compiling CICE, see above "
   exit $res
fi

# Create hycom objects and final hycom_cice executable. 
cd $targetdir
echo "Now compiling hycom_cice in $targetdir." 
env ARCH=$MACROID csh Make_cice.csh 
res=$?
if [ $res -ne 0 ] ; then 
   echo
   echo "Error when compiling HYCOM, see above "
   exit $res
fi
echo "Success"
