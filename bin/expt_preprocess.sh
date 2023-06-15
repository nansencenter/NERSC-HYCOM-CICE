#!/bin/bash
set -x

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
export BINDIR=$(cd $(dirname $0) && pwd)/
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src "; exit 1; }
source $BINDIR/common_functions.sh  || { echo "Could not source common_functions.sh "; exit 1; }

# KAL - new input start and end date times in ISO format
if [ $# -lt 2 ] ; then
   tellerror " Need start and stop times as input"
   exit 1
else 
   starttime=$1
   endtime=$2
   initstr=""
   if [ $# -eq 3 ] ; then
      initstr="$3"
   fi
fi
echo "Start time is $starttime"
echo "Stop  time is $endtime"

# Parse starttime and endtime
if [[ $starttime =~  ([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2}) ]] ; then
   start_year=${BASH_REMATCH[1]}
   echo ${start_year}
   start_month=${BASH_REMATCH[2]}
   start_day=${BASH_REMATCH[3]}
   start_hour=${BASH_REMATCH[4]}
   start_min=${BASH_REMATCH[5]}
   start_sec=${BASH_REMATCH[6]}
   start_hsec=$(echo ${start_min}\*60+${start_sec} | bc )
   start_hsec=$(echo 000${start_hsec} | tail -c5)
   start_dsec=$(echo ${start_hour}\*3600+${start_hsec} | bc )
   start_dsec=$(echo 0000${start_dsec} | tail -c6)
   start_oday=$(date -u -d "${start_year}-${start_month}-${start_day} 00:00:00 UTC" +%j)
else 
   tellerror "start time not in righ format" ; exit 1
fi
if [[ $endtime =~  ([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2}) ]] ; then
   end_year=${BASH_REMATCH[1]}
   end_month=${BASH_REMATCH[2]}
   end_day=${BASH_REMATCH[3]}
   end_hour=${BASH_REMATCH[4]}
   end_min=${BASH_REMATCH[5]}
   end_sec=${BASH_REMATCH[6]}
   end_hsec=$(echo ${end_min}\*60+${end_sec} | bc )
   end_hsec=$(echo 0000${end_hsec} | tail -c5)
   end_dsec=$(echo ${end_hour}\*3600+${end_hsec} | bc )
   end_dsec=$(echo 0000${end_dsec} | tail -c6)
   end_oday=$(date -u -d "$end_year-$end_month-$end_day 00:00:00 UTC" +%j)
else 
   tellerror "end time not in righ format" ; exit 1
fi
#echo hsec $start_hsec $end_hsec
#echo $starttime
#echo dsec $start_dsec $end_dsec
#echo oday $start_oday $end_oday
#exit

# Init error counter (global var used by function tellerror)
numerr=0

# Create data and scratch dirs. Enter scratch dir
if [ ! -d $D ] ; then
   mkdir -p $D   || { echo "Could not create data dir : $S " ; exit 1 ; }
fi
if [ ! -d $S ] ; then
   mkdir -p $S   || { echo "Could not create data dir : $S " ; exit 1 ; }
fi
cd       $S || { tellerror "no scratch dir $S" ;  exit 1 ;}



# --- fetch some things from blkdat.input 
cp $P/blkdat.input blkdat.input || tellerror "No blkdat.input file" 
export LBFLAG=`grep "'lbflag' =" blkdat.input | awk '{printf("%03d", $1)}'`
export EB=`grep "'iexpt ' =" blkdat.input | awk '{printf("%03d", $1)}'`
export PRIVER=`grep "'priver' =" blkdat.input | awk '{printf("%1d", $1)}'`
export NTRACR=`grep "'ntracr' =" blkdat.input | awk '{printf("%03d", $1)}'`
export YRFLAG=`grep "'yrflag' =" blkdat.input | awk '{printf("%1d", $1)}'`
export JERLV=`grep "'jerlv0' =" blkdat.input | awk '{printf("%1d", $1)}'`
export SSSRLX=`grep "'sssflg' =" blkdat.input | awk '{printf("%1d", $1)}'`
export SSTRLX=`grep "'sstflg' =" blkdat.input | awk '{printf("%1d", $1)}'`
export RLX=`grep "'relax ' =" blkdat.input | awk '{printf("%1d", $1)}'`
export TRCRLX=`grep "'trcrlx' =" blkdat.input | awk '{printf("%1d", $1)}'`
export THKDF4=`grep "'thkdf4' =" blkdat.input | awk '{printf("%f", $1)}'`
export VELDF4=`grep "'veldf4' =" blkdat.input | awk '{printf("%f", $1)}'`
export KAPREF=`grep "'kapref' =" blkdat.input | awk '{printf("%f", $1)}'`
export VSIGMA=`grep "'vsigma' =" blkdat.input | awk '{printf("%1d", $1)}'`
export FLXOFF=`grep "'flxoff' =" blkdat.input | awk '{printf("%1d", $1)}'`
export STDFLG=`grep "'stdflg' =" blkdat.input | awk '{printf("%1d", $1)}'`
export BNSTFQ=$(blkdat_get blkdat.input bnstfq)
export NESTFQ=$(blkdat_get blkdat.input nestfq)
export THKDF2=$(blkdat_get blkdat.input thkdf2)
export THFLAG=$(blkdat_get blkdat.input thflag)
export BACLIN=$(blkdat_get blkdat.input baclin)
export BATROP=$(blkdat_get blkdat.input batrop)
export CPLIFQ=$(blkdat_get blkdat.input cplifq)
export ICEFLG=$(blkdat_get blkdat.input iceflg)
export MOMTYP=$(blkdat_get blkdat.input momtyp)
export VISCO2=$(blkdat_get blkdat.input visco2)
export VELDF2=$(blkdat_get blkdat.input veldf2)
export IDM=$(blkdat_get blkdat.input idm)
export JDM=$(blkdat_get blkdat.input jdm)
export LWFLAG=`grep "'lwflag' =" blkdat.input | awk '{printf("%1d", $1)}'`

restarti=$(blkdat_get_string blkdat.input nmrsti "restart_in")

# Add period to restart file name if not present...
if [ $restarti != "restart_in" -a ${restarti: -1} != "."  ] ;then
   restarti="${restarti}."
fi


#
# --- Set up time limits
#
cmd="$BINDIR/hycom_limits.py $starttime $endtime $initstr"
echo "*Setting up HYCOM time limits : $cmd "
eval $cmd ||  tellerror "$cmd failed"
if [ $ICEFLG -eq 2 ] ; then
   cmd="$BINDIR/cice_limits.py $initstr $starttime $endtime $NMPI $P/ice_in"
   echo "*Setting up CICE  time limits : $cmd "
   eval $cmd ||  tellerror "$cmd failed"
fi

# --- fetch some things from the ice model
if [ $ICEFLG -eq 2 ] ; then
   #export icedt=`egrep "^[ ,]*dt[ ]*=" ../ice_in | sed "s/.*=//" | sed "s/,.*//"`
   export icedt=$($HYCOM_PYTHON_ROUTINES/namelist_extract.py ice_in setup_nml dt)
   export ice_restart_dir=$($HYCOM_PYTHON_ROUTINES/namelist_extract.py ice_in setup_nml restart_dir)
   export ice_restart_file=$($HYCOM_PYTHON_ROUTINES/namelist_extract.py ice_in setup_nml restart_file)
   export ice_restart_pointer_file=$($HYCOM_PYTHON_ROUTINES/namelist_extract.py ice_in setup_nml pointer_file)
fi


# --- check that iexpt from blkdat.input agrees with E from this script.
[ "$EB" == "$E" ] || tellerror " blkdat.input iexpt ${EB} different from this experiment ${E} set in EXPT.src"
echo "Fetched from blkdat.input:"
echo "--------------------------"
echo "EB     is $EB    "
echo "PRIVER is $PRIVER"
echo "YRFLAG is $YRFLAG"
echo "JERLV  is $JERLV "
echo "SSS    is $SSSRLX"
echo "SST    is $SSTRLX"
echo "BNSTFQ is $BNSTFQ"
echo "NESTFQ is $NESTFQ"
echo "FLXOFF is $FLXOFF"
echo "STDFLG is $STDFLG"
echo "--------------------"

# Check baroclinic time step 
if [ $YRFLAG -ge 1 ] ; then
   tmp2=21600.
else 
   tmp2=86400.
fi
tmp=$(echo $tmp2"%"$BACLIN | bc)
testbc=$(echo $tmp"=="0 | bc )
if [ $testbc -ne 1 ] ; then
   tellerror "$tmp2 seconds is not divisible with baclin=$BACLIN".
fi

# Check barotropic time step 
tmp=$(echo $BACLIN"/"$BATROP | bc -l)
testbt=$(echo $tmp"%2==0" | bc -l )
if [ $testbt -ne 1 ] ; then
   tellerror "($BACLIN / $BATROP ) %2 not zero ".
fi


# Check coupling time step
if [ $ICEFLG -eq 2 ] ; then
   # Check that icedt = cplifq*baclin
   test1=$(echo ${CPLIFQ}'<'0.0 | bc -l)

   # cplifq negative = number of time steps between coupling
   if [ $test1 -eq 1 ] ; then
      dtcpl=$(echo "-1.*"$CPLIFQ"*"$BACLIN | bc -l)
   # cplifq positive = fraction of days
   else 
      dtcpl=$(echo $CPLIFQ"*"86400. | bc -l)
   fi
   echo "Coupling time step=$dtcpl"

   # Check ice model dt against $dtcpl. dtcpl must be a multiple of icedt
   test2=$(echo "scale=0;"$dtcpl"%"$icedt | bc -l)
   test2=$(echo ${test2}'=='0.0 | bc -l)
   if [ $test2 -eq 0 ] ; then
      tellerror "coupling time step $dtcpl not an integer multiple of ice model time step $icedt"
   fi
fi


# TODO: Limitation for now. Note that newest hycom can initialize with yrflag=3 (realistic forcing)
if [ $YRFLAG -ne 3 ] ; then
   tellerror "must use yrflag=3 in blkdat.input"
   exit 1
fi

# Check momtyp. veldf2 and visco2 must be zero in this case
if [ $MOMTYP -eq 4 ] ; then
   test1=$(echo ${VISCO2}'!='0.0 | bc -l)
   test2=$(echo ${VELDF2}'!='0.0 | bc -l)
   if [ $test1 -ne 0 -o $test2 -ne 0 ] ; then
      tellerror "momtyp=4; visco2 and veldf2 must be zero"
   fi
fi



touch ports.input tracer.input
rm    ports.input tracer.input
[ -s  $P/tracer.input ] && cp $P/tracer.input tracer.input


# Get init flag from start time
if [ "$initstr" == "--init" ] ;then
   init=1
else 
   init=0
fi
tstart=$(cat limits | tr -s " " | sed "s/^[ ]*//" | cut -d " " -f1)
tstop=$(cat limits  | tr -s " " | sed "s/^[ ]*//" | cut -d " " -f2)
echo "Fetched from limits:"
echo "--------------------"
echo "init   is $init"
echo "tstart is $tstart"
echo "tstop  is $tstop"
echo "--------------------"

#C
#C --- turn on detailed debugging.
#C
#touch PIPE_DEBUG

#
# --- pget, pput "copy" files between scratch and permanent storage.
# --- Can both be cp if the permanent filesystem is mounted locally.
# --- KAL - we use pget=pput=cp
export pget=/bin/cp
export pput=/bin/cp
export plink='ln -sf'

echo "Initialization complete - now copying necessary files to scratch area"
echo

#
#
# --- Set up MPI partition for this segment
#
#
#echo "we go to copy the damn file",$NMPI
[  "$NMPI" == "" ] && NMPI=0  # Init if unset
if [ $NMPI != 0 ] ; then
  export NPATCH=`echo $NMPI | awk '{printf("%04d", $1)}'`
  /bin/rm -f patch.input
  echo /bin/cp $BASEDIR/topo/partit/depth_${R}_${T}.${NPATCH}  patch.input
  /bin/cp $BASEDIR/topo/partit/depth_${R}_${T}.${NPATCH}  patch.input || \
     tellerror "Could not get patch.input  ($BASEDIR/topo/partit/depth_${R}_${T}.${NPATCH})"
fi

#
# --- input files from file server.
#
echo "**Retrieving grid and topography files"
${pget} $BASEDIR/topo/regional.grid.a regional.grid.a || tellerror "no grid file regional.grid.a" 
${pget} $BASEDIR/topo/regional.grid.b regional.grid.b || tellerror "no grid file regional.grid.a" 
${pget} $BASEDIR/topo/depth_${R}_${T}.a regional.depth.a || tellerror "no topo file depth_${R}_${T}.a" 
${pget} $BASEDIR/topo/depth_${R}_${T}.b regional.depth.b || tellerror "no topo file depth_${R}_${T}.b" 
${pget} $BASEDIR/topo/kmt_${R}_${T}.nc cice_kmt.nc     || tellerror "no kmt file $BASEDIR/topo/kmt_${R}_${T}.nc "
${pget} $BASEDIR/topo/cice_grid.nc cice_grid.nc        || tellerror "no cice grid file $BASEDIR/topo/cice_grid.nc "


if [ "$SSTRLX" -eq 3 ] ; then
   [ -f  $CLMDIR/seatmp.a ] || tellerror "File $CLMDIR/seatmp.a does not exist"
   [ -f  $CLMDIR/seatmp.b ] || tellerror "File $CLMDIR/seatmp.b does not exist"
   ln -sf $CLMDIR/seatmp.a forcing.seatmp.a || tellwarn "Could not link $CLMDIR/seatmp.a"
   ln -sf $CLMDIR/seatmp.b forcing.seatmp.b || tellwarn "Could not link $CLMDIR/seatmp.b"
fi

# TODO: Not handled yet
if [ 0 -eq 1 ] ; then
   [ -f  $CLMDIR/surtmp.a ] || tellerror "File $CLMDIR/surtmp.a does not exist"
   [ -f  $CLMDIR/surtmp.b ] || tellerror "File $CLMDIR/surtmp.b does not exist"
   ln -sf $CLMDIR/surtmp.a forcing.surtmp.a || tellwarn "Could not link $CLMDIR/surtmp.a"
   ln -sf $CLMDIR/surtmp.b forcing.surtmp.b || tellwarn "Could not link $CLMDIR/surtmp.b"
fi

# ---
# --- Pre-prepared forcing option - for now just sets up links
# --- TODO: To verify, we have to go into file and check timings
# ---
echo "**Setting up pre-prepared synoptic forcing from force/synoptic/$E"
DIR=$BASEDIR/force/synoptic/$E/
echo $DIR

declare -a arr=("radflx" "shwflx" "vapmix" "airtmp" "precip" "mslprs" "wndewd" "wndnwd" "dewpt")
#for i in tauewd taunwd wndspd radflx shwflx vapmix \
#   airtmp precip uwind vwind clouds relhum slp ; do

#for i in radflx shwflx vapmix \
#   airtmp precip mslprs \
#   wndewd wndnwd ; do
for i in "${arr[@]}"; do
   echo "|--> $i"
   [ -f  $DIR/$i.a ] || tellerror "File $DIR/$i.a does not exist"
   [ -f  $DIR/$i.b ] || tellerror "File $DIR/$i.b does not exist"

   # Check range of file against start and stop times
   if [ -f  $DIR/$i.a -a  -f $DIR/$i.b ] 
   then
      [ ! -s  $DIR/${i}.a  ] && tellerror "$DIR/$i.a: File is empty"
      [ ! -s  $DIR/${i}.b  ] && tellerror "$DIR/$i.b: File is empty"
      ln -sf $DIR/${i}.a forcing.${i}.a ||  tellerror "Could not fetch $DIR/$i.a"
      ln -sf $DIR/${i}.b forcing.${i}.b ||  tellerror "Could not fetch $DIR/$i.b"

      frcstart=$(head -n  6 forcing.$i.b | tail -n1 | sed "s/.*=//" | sed "s/^[ ]*//" | cut -d " " -f 1)
      frcstop=$(tail -n 1 forcing.$i.b             | sed "s/.*=//" | sed "s/^[ ]*//" | cut -d " " -f 1)
      test1=$(echo ${tstart#-}'>='$frcstart | bc -l)
      test2=$(echo $tstop '<='$frcstop  | bc -l)
      [ $test1 -eq 1 ] || tellerror "File $S/forcing.$i.b: forcing starts after  model starts"
      [ $test2 -eq 1 ] || tellerror "File $S/forcing.$i.b: forcing stops  before model stops"
   fi
      

done


# MOSTAFA: END
#
# --- time-invarent heat flux offset
#
#TODO handle offsets
#setenv OF ""
#setenv OF "_413"
if [ "$CLMDIR" != "" ] ; then
   if [  "$OF" != "" ];  then
     touch  forcing.offlux.a  forcing.offlux.b
     if [ ! -s  forcing.offlux.a ] ; then
        ${pget} $CLMDIR/offlux${OF}.a forcing.offlux.a &
     fi
     if [ ! -s  forcing.offlux.b ] ; then
        ${pget} $CLMDIR/offlux${OF}.b forcing.offlux.b &
     fi 
   fi 
fi 

# MOSTAFA: BEGIN
# For time-invariant offlux
# copy flux off set files if flxoff=1
echo "FLXOFF =  $FLXOFF"
if [ $FLXOFF -eq 1 ] ; then
 echo "===================================================="
 echo " -------flux off set true: copy flux off set files-"
   cp ${D}/../../relax/${E}/offlux.a forcing.offlux.a || tellerror "Could not get offlux .a file"
   cp ${D}/../../relax/${E}/offlux.b forcing.offlux.b || tellerror "Could not get offlux .b file"
 echo "===================================================="
 else
    echo "fLxoff=F: No attempt to use flux offset correction" 
fi


# MOSTAFA: END

#
# --- river forcing
# --- KAL: rivers are experiment-dependent
#
if [ $PRIVER -eq 0 ] ; then
echo "**No river forcing. Set the priver to 1 to add river forcing"
fi
if [ $PRIVER -eq 1 ] ; then
  echo "**Setting up river forcing  from priver"
  cp $BASEDIR/force/rivers/$E/rivers.a forcing.rivers.a || tellerror "Could not get river .a file"
  cp $BASEDIR/force/rivers/$E/rivers.b forcing.rivers.b || tellerror "Could not get river .b file"
   if [ $NTRACR -ne 0 ] ; then
      echo "**Setting up bio river forcing"
      cp $BASEDIR/force/rivers/$E/ECO_no3.a rivers.ECO_no3.a || tellwarn "Could not get NO3 river .a file"
      cp $BASEDIR/force/rivers/$E/ECO_no3.b rivers.ECO_no3.b || tellwarn "Could not get NO3 river .b file"
      cp $BASEDIR/force/rivers/$E/ECO_sil.a rivers.ECO_sil.a || tellwarn "Could not get SIL river .a file"
      cp $BASEDIR/force/rivers/$E/ECO_sil.b rivers.ECO_sil.b || tellwarn "Could not get SIL river .b file"
      cp $BASEDIR/force/rivers/$E/ECO_pho.a rivers.ECO_pho.a || tellwarn "Could not get PHO river .a file"
      cp $BASEDIR/force/rivers/$E/ECO_pho.b rivers.ECO_pho.b || tellwarn "Could not get PHO river .b file"
   fi
fi


#
# --- kpar forcing
#
if [ $JERLV -eq 0 ] ; then
   echo "**Setting up kpar forcing"
   ln -sf $BASEDIR/force/seawifs/kpar.a forcing.kpar.a || tellerror "Could not get kpar.a file"
   ln -sf $BASEDIR/force/seawifs/kpar.b forcing.kpar.b || tellerror "Could not get kpar.b file"
fi


#
# --- relaxation
#
if [ $SSSRLX -eq 1 -o $SSTRLX -eq 1 -o $RLX -eq 1 -o $init -eq 1 ] ; then
   echo "**Setting up relaxation"
   for i in saln temp intf ; do
      j=$(echo $i | head -c3)
      [ ! -f  $BASEDIR/relax/${E}/relax_$j.a ] && tellerror "$BASEDIR/relax/${E}/relax_$j.a does not exist"
      [ ! -f  $BASEDIR/relax/${E}/relax_$j.b ] && tellerror "$BASEDIR/relax/${E}/relax_$j.b does not exist"
      ln -sf $BASEDIR/relax/${E}/relax_$j.a relax.$i.a  || tellerror "Could not get relax.$i.a"
      ln -sf $BASEDIR/relax/${E}/relax_$j.b relax.$i.b  || tellerror "Could not get relax.$i.b"
   done
fi
if [ $RLX -eq 1 ] ; then
   echo "**Setting up relaxation masks"
   [ ! -f  $BASEDIR/relax/${E}/relax_rmu.a ] && tellerror "$BASEDIR/relax/${E}/relax_rmu.a does not exist"
   [ ! -f  $BASEDIR/relax/${E}/relax_rmu.b ] && tellerror "$BASEDIR/relax/${E}/relax_rmu.b does not exist"
   ln -sf $BASEDIR/relax/${E}/relax_rmu.a relax.rmu.a  || tellerror "Could not get relax.rmu.a"
   ln -sf $BASEDIR/relax/${E}/relax_rmu.b relax.rmu.b  || tellerror "Could not get relax.rmu.b"
fi
#
# --- tracer relaxation
#
if [ $TRCRLX -ne 0 -o $NTRACR -eq -1 ] ; then
   echo "**Setting up tracer relaxation"
   for i in ECO_no3 ECO_pho ECO_sil ECO_oxy CO2_dic CO2_alk; do
      j=$(echo $i | head -c7)
      [ ! -f  $BASEDIR/relax/${E}/relax.$j.a ] && tellerror "$BASEDIR/relax/${E}/relax.$j.a does not exist"
      [ ! -f  $BASEDIR/relax/${E}/relax.$j.b ] && tellerror "$BASEDIR/relax/${E}/relax.$j.b does not exist"
      ln -sf $BASEDIR/relax/${E}/relax.$j.a relax.$i.a  || tellerror "Could not get relax.$i.a"
      ln -sf $BASEDIR/relax/${E}/relax.$j.b relax.$i.b  || tellerror "Could not get relax.$i.b"
   done
   echo "**Setting up tracer relaxation masks"
   [ ! -f  $BASEDIR/relax/${E}/relax_rmu.a ] && tellerror "$BASEDIR/relax/${E}/relax_rmutr.a does not exist"
   [ ! -f  $BASEDIR/relax/${E}/relax_rmu.b ] && tellerror "$BASEDIR/relax/${E}/relax_rmutr.b does not exist"
   ln -sf $BASEDIR/relax/${E}/relax_rmu.a relax.rmutr.a  || tellerror "Could not get relax.rmutr.a"
   ln -sf $BASEDIR/relax/${E}/relax_rmu.b relax.rmutr.b  || tellerror "Could not get relax.rmutr.b"

   [ ! -f  $INPUTDIR/co2_annmean_gl.txt ] && tellerror "$INPUTDIR/co2_annmean_gl.txt does not exist"
   ln -sf $INPUTDIR/co2_annmean_gl.txt co2_annmean_gl.txt || tellerror "Could not get co2_annmean_gl.txt"
fi
#
# - thermobaric reference state?
# 
echo "**Setting up various forcing files"
[ -f tbaric.a ] && rm tbaric.a
[ -f tbaric.b ] && rm tbaric.b
if [ ${KAPREF:0:1} == "-" ] ;then
   ${pget} $BASEDIR/topo/tbaric.a tbaric.a  || tellerror "Could not get tbaric.a"
   ${pget} $BASEDIR/topo/tbaric.b tbaric.b  || tellerror "Could not get tbaric.b"
fi

#
# Spatially varying isopycnal densities
#
if [ ${VSIGMA} -ne  0 ] ;then
   ${pget} $BASEDIR/relax/${E}/iso_sigma.a iso.sigma.a  || tellerror "Could not get iso.sigma.a"
   ${pget} $BASEDIR/relax/${E}/iso_sigma.b iso.sigma.b  || tellerror "Could not get iso.sigma.b"
fi


# Thickness diffusion
testthkdf2=$(echo $THKDF2'<'0.0 | bc -l)
testthkdf4=$(echo $THKDF4'<'0.0 | bc -l)
[ -f thkdf4.a ] && rm thkdf4.a
[ -f thkdf4.b ] && rm thkdf4.b
if [ ${testthkdf4} -eq 1 ] ; then 
   ${pget} ${D}/../../relax/${E}/thkdf4.a thkdf4.a  || tellerror "Could not get thkdf4.a"
   ${pget} ${D}/../../relax/${E}/thkdf4.b thkdf4.b  || tellerror "Could not get thkdf4.b"
fi

[ -f thkdf2.a ] && rm thkdf2.a
[ -f thkdf2.b ] && rm thkdf2.b
if [ ${testthkdf2} -eq 1 ] ; then 
   ${pget} ${D}/../../relax/${E}/thkdf2.a thkdf2.a  || tellerror "Could not get thkdf2.a"
   ${pget} ${D}/../../relax/${E}/thkdf2.b thkdf2.b  || tellerror "Could not get thkdf2.b"
fi
testveldf4=$(echo $VELDF4'<'0.0 | bc -l)
if [ ${testveldf4} -eq 1 ] ; then 
   ${pget} ${D}/../../relax/${E}/veldf4.a veldf4.a  || tellerror "Could not get veldf4.a"
   ${pget} ${D}/../../relax/${E}/veldf4.b veldf4.b  || tellerror "Could not get veldf4.b"
fi


# TODO Limited set of tests for now. 
# Link in nest dir if nesting activated
#tmp=$(echo $BNSTFQ'>'0.0 | bc -l)
#tmp2=$(echo $NESTFQ'>'0.0 | bc -l)
tmp=$(echo $BNSTFQ'!='0.0 | bc -l)
tmp2=$(echo $NESTFQ'!='0.0 | bc -l)
if [ $tmp -eq 1 -o $tmp2 -eq 1 ] ; then
   nestdir=$BASEDIR/nest/$E
   echo "Nesting input from $nestdir"
   ls nest
   if [ -d $nestdir ]  ; then
      [ -e nest ] && rm nest
      ln -s $nestdir nest
   else 
      tellerror "Nesting dir $nest does not exist"
   fi
fi

#export waveSDIR=/work/shared/nersc/msc/STOKES/Globww3/tmp
#export  waveSDIR=/work/shared/nersc/msc/STOKES/Globww3
echo "===================================================="
echo "STDFLG =  $STDFLG"
if [ $STDFLG -eq 1 ] ; then
 echo "===================================================="
 echo " -------Setting up Wave Stokes forcing----"
 for foo in  forcing.stokex.a forcing.stokex.b forcing.stokey.a forcing.stokey.b\
     forcing.transx.a forcing.transx.b forcing.transy.a forcing.transy.b\
     forcing.tauwx.a  forcing.tauwx.b forcing.tauwy.a  forcing.tauwy.b \
     forcing.twomx.a  forcing.twomx.b forcing.twomy.a  forcing.twomy.b \
     forcing.t01.a    forcing.t01.b     ; do

    echo "|--> linking to $waveSDIR/${foo}"
     ln -sf $waveSDIR/${foo} . ||  tellerror "Could not fetch stokes forcing"
 done
 echo "===================================================="
 else
    echo "STD=F: No attempt to link to the Stokes/Wave forcing in SCRATCH" 
fi
#read -t 5 
# Need ports.input file in these cases
if [ $tmp -eq 1 -a $LBFLAG -ne 2 -a $LBFLAG -ne 4 ] ; then
      tellerror "Must have lbflag = 2 or 4 when bnstfq <> 0.0 "
elif [ $tmp -eq 1 -a $LBFLAG -eq 1 ] ; then
   # Port flow - file must be present in experiment dir
   cp $P/ports.input . || tellerror "Could not get port file ${P}/ports.input for port flow"
elif [ $tmp -eq 1 -a $LBFLAG -eq 2 ] ; then
   # Nest flow - use file in  experiment dir if present. Otherwise look in nest dir
   if [ -f $P/ports.nest ] ; then
      echo "Using file $P/ports.nest for nesting: $P/ports.nest -> ./ports.input"
      cp $P/ports.nest ports.input       || tellerror "Could not get port file ${P}/ports.nest for nest flow"
   elif [ -f $nestdir/ports.nest ] ; then
      echo "Using file $nestdir/ports.nest for nesting: $P/ports.nest -> ./ports.input"
      cp $nestdir/ports.nest ports.input || tellerror "Could not get port file ${nestdir}/ports.nest for nest flow"
   else 
      tellerror "Could not get port file ports.nest in $P or  ${nestdir} for nest flow"
   fi
fi

# Need nest rmu in this case:
if [ $tmp2 -eq 1  ] ; then
   # Nest relaxation - use file in  experiment dir if present. Otherwise look in nest dir
#   if [ -f $P/rmu_nest.a -a -f $P/rmu_nest.a ] ; then
#      echo "Using file $P/rmu_nest.[ab] for nesting relaxation: $P/rmu_nest.[ab] -> ./rmu.[ab]"
#      cp $P/rmu_nest.a rmu.a       || tellerror "Could not get port file ${P}/rmu_nest.a for nest relax"
#      cp $P/rmu_nest.b rmu.b       || tellerror "Could not get port file ${P}/rmu_nest.b for nest relax"
#   elif [ -f $nestdir/rmu_nest.a -a -f $nestdir/rmu_nest.a ] ; then
#      echo "Using file $nestdir/rmu_nest.[ab] for nesting: $nestdir/rmu_nest.[ab] -> ./rmu.[ab]"
#      cp $nestdir/rmu_nest.a rmu.a       || tellerror "Could not get port file ${nestdir}/rmu_nest.a for nest relax"
#      cp $nestdir/rmu_nest.b rmu.b       || tellerror "Could not get port file ${nestdir}/rmu_nest.b for nest relax"
#   fi
  if [ -f $nestdir/rmu.a -a -f $nestdir/rmu.b ] ; then
      echo "Using file $nestdir/rmu.[ab] for nesting"
   else 
      tellerror "Could not find files $nestdir/rmu.[ab] for nest relaxation"
   fi
fi


## Retrieve nersc tidal data set if active
#if [ "$tideflag" == "T" ] ; then
#   ${pget} ${BASEDIR}/tides_nersc/$E/${tidechoice}obc_elev.dat . || tellerror "Could not get tidal data ${tidechoice}obc_elev.dat "
#fi
   



echo

#
# --- move old restart files to KEEP, typically from batch system rerun.
#
touch restart.dummy.b restart.dummy.a
for f  in restart* ; do
  /bin/mv $f KEEP/$f
done

#
# --- try to get HYCOM restart from various areas
#
if [ $init -eq 1 ] ; then
   echo "No restart needed"
else

   #HYCOM restart
   filename=${restarti}${start_year}_${start_oday}_${start_hour}_${start_hsec}
   echo $D/${filename}_mem001.a

   # Try to fetch restart from data dir $D
   if [ -f $D/${filename}.a -a -f $D/${filename}.b ] ; then
      echo "using HYCOM restart files ${filename}.[ab] from data dir $D"
      cp $D/${filename}.a .
      cp $D/${filename}.b .

   elif [ -f $D/${filename}_mem001.a -a -f $D/${filename}_mem001.b ]; then
      echo "using HYCOM restart files ${filename}_mem???.[ab] from data dir $D"
      for f in ${plink} $D/${filename}_mem*.? ; do
         ${plink} $f .
      done

   else
      tellerror "Could not find HYCOM restart file ${filename}.[ab] in $D"
   fi

   #CICE restart
   if [ $ICEFLG -eq 2 ] ; then
      filenameice="${ice_restart_dir}/${ice_restart_file}.${start_year}-${start_month}-${start_day}-${start_dsec}"

      # Try to fetch restart from data dir $D
      if [ -f $D/${filenameice}.nc ] ; then
         echo "using CICE restart file ${filenameice}.nc from data dir $D"
         cp $D/${filenameice}.nc ${filenameice}.nc
         echo ${filenameice}.nc > ${ice_restart_pointer_file}

      elif [ -f $D/${filenameice}_mem001.nc ]; then
         echo "using CICE restart file ${filenameice}_mem???.nc from data dir $D"
         for f in $D/${filenameice}_mem*.nc ; do
            ${plink} $f cice/.
         done
         echo ${filenameice}_mem000.nc > ${ice_restart_pointer_file}

      else
         tellerror "Could not find CICE restart file ${filenameice} in $D"
      fi
   fi


fi



#
# --- model executable. One executable to rule them all (if it wasnt for CICE, that is...)
#
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
else
   echo "SIGVER = $SIGVER"
   echo "So far only 7 term eq of state is supported (SIGVER=1 or 2) is supported"
   exit 1
fi
TERMS2=$(echo 0$TERMS | tail -c3)
echo "SIGVER      = $SIGVER .There are $TERMS terms in equation of state"

# Set up rel path and stmt fnc
compdir=$(source_dir $V $TERMS $THFLAG)
compdir=$P/build/${compdir}
if [ $ICEFLG -eq 2 ] ; then
  echo "Retrieving  hycom_cice from $compdir"
  /bin/cp $compdir/hycom_cice  . || tellerror "Could not get hycom_cice executable at "
elif [ $ICEFLG -eq 0 ] ; then
  echo "Retrieving  hycom_oasis from $compdir"
  /bin/cp $compdir/hycom_oasis  . || tellerror "Could not get hycom_oasis executable at "
fi


#
# --- summary printout
#
[ ! -d old ] && mkdir -p old
touch   summary_out
/bin/mv summary_out old/summary_old


#
# --- heat transport output
#
touch   flxdp_out.a flxdp_out.b
/bin/mv flxdp_out.a old/flxdp_old.a
/bin/mv flxdp_out.b old/flxdp_old.b
#
touch   ovrtn_out
/bin/mv ovrtn_out old/ovrtn_old


#
# --- clean up old archive files, typically from batch system rerun.
#
[ ! -d KEEP ] && mkdir KEEP
touch archv.dummy.b archv.dummy.a
for f  in arch* ; do
  /bin/mv $f KEEP/
done


#
# --- let all file copies complete.
#
wait
#
# --- zero file length means no rivers.
#
if [ ! -s forcing.rivers.a ] ; then
   /bin/rm forcing.rivers.[ab]
   echo "problem with rivers"
   ls -l *rivers*
fi


#C
#C --- Nesting input archive files for next segment.
#C
#if (-e ./nest) then
#  cd ./nest
#  touch archv_${NB}.tar
#  if (-z archv_${NB}.tar) then
#    ${pget} ${D}/nest/archv_${NB}.tar archv_${NB}.tar &
#  endif
#  cd ..
#endif
#C
#chmod ug+x hycom
#/bin/ls -laFq
#C
#if (-e ./nest) then
#  ls -laFq nest
#endif


if [ $numerr -eq 0 ] ; then
   echo "No fatal errors. Ok to start model set up in $S"
else
   echo "Some fatal errors occured. See above"
fi

# Tell where stuff ended up
exit $numerr # Fails if any fatal errors occured


