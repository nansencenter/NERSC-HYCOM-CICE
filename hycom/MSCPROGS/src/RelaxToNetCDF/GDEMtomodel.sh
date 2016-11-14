# Use 
PROGDIR="$(dirname $0)/"
PROGNAME="$(basename $0)"
MONTH=$1

if [ $# -ne 1 ] ; then
   echo "Script to convert from hycom GDEM climatology to netcdf "
   echo "file at a given projection, supplied by proj.in"
   echo
   echo "Usage :  $PROGNAME monthnumber"
   exit 1
fi



##################################################################################
# First -- Convert GDEM (in hycom format) to z levels on horizontal model grid
##################################################################################
KSIGMA=$(cat blkdat.input | egrep "'thflag'" | cut -f 1 | tr -d " ")
INTERP=1 # Interpolation method: 1=piecewice linear , 2=Cubic spline
ITEST=0
JTEST=0
ARCH=intel
ZINT="/home/nersc/knutali/Progs/HYCOM_UTILITY_ROUTINES/HYCOM_2.1.03/ALL/relax/src/z_gdem.$ARCH"

if [ ! -f $ZINT ] ; then
   echo "Cant find interpolation routine $ZINT"
   exit
fi


if [ "${HYCOMCLIM_PATH}" == "" ] ; then
   echo "Specify path to hycom climatology in environment variable HYCOMCLIM_PATH"
   exit
fi

if [ ! -d "${HYCOMCLIM_PATH}" ] ; then
   echo "${HYCOMCLIM_PATH} is not a directory"
   exit
fi

# This is done to avoid any overwriting in the link below
touch *.d ; rm *.d

# Link climatology files in here
ln -s ${HYCOMCLIM_PATH}/* .


ICTYPE=1 # Input file type -- 1= Rho0&T, 2=Rho0&S
while [ $ICTYPE -le 2 ] ; do

   month2=$(echo "0$MONTH" | tail -c3)
   [ -f "r_m${month2}.d" ] && ln -s r_m${month2}.d fort.71
   [ -f "t_m${month2}.d" ] && ln -s t_m${month2}.d fort.72
   [ -f "s_m${month2}.d" ] && ln -s s_m${month2}.d fort.73

   $ZINT << EOF
      &AFTITL CTITLE='GDEM' /
      &AFFLAG ICTYPE=$ICTYPE , KSIGMA=$KSIGMA , INTERP=$INTERP , MONTH=$MONTH , ITEST=$ITEST , JTEST=$JTEST /
EOF

   # Move output files to new names
   #echo "|$KSIGMA|"
   [ -s  fort.71 ] && rm fort.71
   [ -s  fort.72 ] && rm fort.72
   [ -s  fort.73 ] && rm fort.73
   [ ! -f fort.010a ] && echo "$0: Error Line No $LINENO" && exit
   [ ! -f fort.10   ] && echo "$0: Error Line No $LINENO" && exit
   [ ! -f fort.011a ] && echo "$0: Error Line No $LINENO" && exit
   [ ! -f fort.11   ] && echo "$0: Error Line No $LINENO" && exit
   [ ! -f fort.012a ] && echo "$0: Error Line No $LINENO" && exit
   [ ! -f fort.12   ] && echo "$0: Error Line No $LINENO" && exit

   # for $ICTYPE 1 , use temp fields
   if [ $ICTYPE -eq 1 ] ; then
      mv fort.010a temp_sig${KSIGMA}_m${month2}.a 
      mv fort.10   temp_sig${KSIGMA}_m${month2}.b
   elif [ $ICTYPE -eq 2 ] ; then
      mv fort.011a saln_sig${KSIGMA}_m${month2}.a
      mv fort.11   saln_sig${KSIGMA}_m${month2}.b
   fi

   [ -f fort.010a ] && rm fort.010a
   [ -f fort.10   ] && rm fort.10
   [ -f fort.011a ] && rm fort.011a
   [ -f fort.11   ] && rm fort.11
   [ -f fort.012a ] && rm fort.012a
   [ -f fort.12   ] && rm fort.12

   let ICTYPE=$ICTYPE+1
done




##################################################################################
# The temp_sig${KSIGMA}_m${month2}.a  etc files now contain values at fixed z-levels.
# Find what indices to read, depending on input
##################################################################################

# This little stub gets depth levels in the hycom GDEM files (the .d ones),
# and returns them in hclim_depthlevels
[ -f hclim_depthlevels ] && rm hclim_depthlevels
$PROGDIR/hclimlevels


# This program requires two inputs, month, and depth level. The depth level
# is matched to correct record in data file, using hclim_depthlevels.
#
# The data is read, and dumped in a cf-1.0 compliant netcdf file - one per month
# TODO:  use depthlevels.in for this ?
depthlevels="30.0 50.0  100.0 200.0 400.0 700.0 1000.0 1500.0 2000.0 2600.0 3000.0"

$PROGDIR/hclimtonc $MONTH $depthlevels

