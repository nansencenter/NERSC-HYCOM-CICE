#!/bin/bash

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
echo "Stop  time is $stoptime"

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
export YRFLAG=`grep "'yrflag' =" blkdat.input | awk '{printf("%1d", $1)}'`
export JERLV=`grep "'jerlv0' =" blkdat.input | awk '{printf("%1d", $1)}'`
export SSSRLX=`grep "'sssflg' =" blkdat.input | awk '{printf("%1d", $1)}'`
export SSTRLX=`grep "'sstflg' =" blkdat.input | awk '{printf("%1d", $1)}'`
export RLX=`grep "'relax ' =" blkdat.input | awk '{printf("%1d", $1)}'`
export THKDF4=`grep "'thkdf4' =" blkdat.input | awk '{printf("%f", $1)}'`
export KAPREF=`grep "'kapref' =" blkdat.input | awk '{printf("%f", $1)}'`
export VSIGMA=`grep "'vsigma' =" blkdat.input | awk '{printf("%1d", $1)}'`
export BNSTFQ=$(blkdat_get blkdat.input bnstfq)
export NESTFQ=$(blkdat_get blkdat.input nestfq)
export THKDF2=$(blkdat_get blkdat.input thkdf2)
export THFLAG=$(blkdat_get blkdat.input thflag)
export BACLIN=$(blkdat_get blkdat.input baclin)
export BATROP=$(blkdat_get blkdat.input batrop)
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
tmp=$(echo $BATROP"/"$BACLIN | bc -l)
testbt=$(echo $tmp"%2==0" | bc -l )
if [ $testbt -ne 1 ] ; then
   tellerror "($BACLIN / $BATROP ) %2 not zero ".
fi

# TODO: Limitation for now. Note that newest hycom can initialize with yrflag=3 (realistic forcing)
if [ $YRFLAG -ne 3 ] ; then
   tellerror "must use yrflag=3 in blkdat.input"
   exit 1
fi

#
# --- Set up time limits
#
echo "*Setting up time limits"
cmd="$BASEDIR/../python/hycom_limits.py $starttime $endtime $initstr"
eval $cmd ||  tellerror "$cmd failed"
cmd="$BASEDIR/../python/cice_limits.py $initstr $starttime $endtime $NMPI $P/ice_in"
eval $cmd ||  tellerror "$cmd failed"


touch ports.input tracer.input
rm    ports.input tracer.input
[ -s  $P/tracer.input ] && cp $P/tracer.input tracer.input



##
## --- Fetch some things from infile.in
##
## TODO iday may be ordinal day > days in year
#rungen=$(head -n 2 infile.in  | tail -n1 | cut -c1-3)
#refyear=$(head -n 3 infile.in  | tail -n1 | cut -c1-4)
#refyear=$(printf "%4.4d" $refyear)
#iday=$(head -n 4 infile.in  | tail -n1 | cut -c1-4)
#iday=`echo 00$iday | tail -4c`
#ihour=$(head -n 4 infile.in  | tail -n1 | cut -c7-8)
#ihour=`echo 0$ihour | tail -3c`
#frc=$(head -n 6 infile.in  | tail -n1 | cut -c1-5 | tr -d " ")
#clm=$(head -n 6 infile.in  | tail -n1 | cut -c7-11 | tr -d " ") 
#nestoflag=$(head -n 11 infile.in | tail -n1 | awk '{printf("%1s", $1)}')
#nestiflag=$(head -n 12 infile.in | tail -n1 | awk '{printf("%1s", $1)}')
#tideflag=$(head -n 13 infile.in | tail -n1 | cut -c1 )
#tidechoice=$(head -n 13 infile.in | tail -n1 | cut -c3-5 )
#gpflag=$(head -n 14 infile.in | tail -n1 | awk '{printf("%1s", $1)}')
#echo "Fetched from infile.in:"
#echo "-----------------------"
#echo "rungen  is              : $rungen"
#echo "refyear is              : $refyear"
#echo "iday    is              : $iday"
#echo "ihour   is              : $ihour"
#echo "Forcing option is       : ${frc}"
#echo "Climatology option is   : ${clm}"
#echo "Outer nest flag is      : $nestoflag"
#echo "Inner nest flag is      : $nestiflag"
#echo "Tide flag and choice is : $tideflag $tidechoice"
#echo "Gidpoint flag is        : $gpflag "
#echo


##
## --- Fetch some things from limits
##
#tstart=$(cat limits | tr " " -s | sed "s/^[ ]*//" | cut -d " " -f1)
#tstop=$(cat limits | tr " " -s | sed "s/^[ ]*//" | cut -d " " -f2)

## Get init flag from start time
#init=0
#[ `echo $tstart | awk '{print ($1 < 0.0)}'` == 1  ] && init=1
#echo "Fetched from limits:"
#echo "--------------------"
#echo "init is $init"
#echo


# Get init flag from start time
if [ $initstr == "--init" ] ;then
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
echo "**Retrieving grid ant topography files"
${pget} $BASEDIR/topo/regional.grid.a regional.grid.a || tellerror "no grid file regional.grid.a" 
${pget} $BASEDIR/topo/regional.grid.b regional.grid.b || tellerror "no grid file regional.grid.a" 
${pget} $BASEDIR/topo/depth_${R}_${T}.a regional.depth.a || tellerror "no topo file depth_${R}_${T}.a" 
${pget} $BASEDIR/topo/depth_${R}_${T}.b regional.depth.b || tellerror "no topo file depth_${R}_${T}.b" 
${pget} $BASEDIR/topo/kmt_${R}_${T}.nc cice_kmt.nc     || tellerror "no kmt file $BASEDIR/topo/kmt_${R}_${T}.nc "
${pget} $BASEDIR/topo/cice_grid.nc cice_grid.nc        || tellerror "no cice grid file $BASEDIR/topo/cice_grid.nc "


##
## --- Check forcing and climatology option in infile.in
##
#if [ "$clm" == "era40" ] ; then
#   CLMDIR=$BASEDIR/force/nersc_era40/$E/
#elif [ "$clm" == "old" ] ; then
#   CLMDIR=$BASEDIR/force/nersc_old/$E/
#elif [ "$clm" == "ncepr" ] ; then
#   CLMDIR=$BASEDIR/force/nersc_ncepr/$E/
#elif [ "$clm" != "prep" ] ; then
#   tellerror "Unknown climate option $clm"
#   CLMDIR=""
#fi

   

##
##
## --- Climatology atmospheric forcing - always copied
##
##
#if [ "$CLMDIR" != "" -a "$clm" != "prep" ] ; then
#   echo "**Retrieving climatology forcing files"
#   for i in tauewd taunwd wndspd radflx shwflx vapmix \
#      airtmp precip uwind vwind clouds relhum slp ; do
#      # TODO - doesnt handle offsets
#      #${pget} $CLMDIR/$i.a      forcing.$i.a || tellerror "Could not fetch $CLMDIR/$i.a"
#      [ -f  $CLMDIR/$i.a ] || tellerror "File $CLMDIR/$i.a does not exist"
#      [ -f  $CLMDIR/$i.b ] || tellerror "File $CLMDIR/$i.b does not exist"
#      ln -sf $CLMDIR/$i.a      forcing.$i.a || tellerror "Could not link $CLMDIR/$i.a"
#      ln -sf $CLMDIR/$i.b      forcing.$i.b || tellerror "Could not link $CLMDIR/$i.b"
#   done
#fi


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

##
## --- Synoptic forcing option - for now just sets up links
##
#if [ "$frc" != "month" -a  "$clm" != "prep" ] ; then
#   echo "**Setting up INLINE_FORCING synoptic forcing"
#
#   # pathvar is name of variable
#   pathvar=""
#   if [ "${frc}" == "era40" ] ; then
#      pathvar="ERA40_PATH"
#   elif [ "${frc}" == "era-i" ] ; then
#      pathvar="ERAI_PATH"
#   elif [ "${frc}" == "ncepr" ] ; then
#      pathvar="NCEP_PATH"
#   elif [ "${frc}" == "ecmwf" ] ; then
#      pathvar="ECMWF_PATH"
#      pathvar2="./Ecmwfr"
#   elif [ "${frc}" == "metno" ] ; then
#      pathvar="METNO_PATH"
#      pathvar2="./Met.no"
#   elif [ "${frc}" == "ecnc" ] ; then
#      pathvar="ECNC_PATH"
#      pathvar2="./Ecmwf.nc"
#   elif [ "${frc}" == "ncepr" ] ; then
#      pathvar="NCEP_PATH"
#   else
#      tellerror "No method for forcing $frc (did you set it up?)"
#   fi
#
#   if [ "$pathvar" != "" ] ; then
#      pathval="${!pathvar}" # pathval is value of variable named pathvar
#      echo "Forcing is $frc, environment variable $pathvar is $pathval"
#      if [ "${pathval}" == "" -o ! -d ${pathval:=tomtomtom} ] ; then  #tomtomtom is substituted if pathval is empty
#         tellerror "$pathvar not set or not a directory"
#      else # Could be present
#         if [ "$pathvar2" != ""  ] ; then
#            ln -s ${pathval} $pathvar2 || tellwarn "Could not link $pathvar ${pathval} as $pathvar2 (it might be present already)"
#         else
#            ln -s ${pathval} . || tellwarn "Could not link $pathvar ${pathval} (it might be present)"
#         fi
#      fi
#   else
#      tellerror "No pathvar $pathvar set up for forcing $frc. Set it in REGION.src "
#   fi
##else
##   tellerror "Forcing variable frc undefined"
#fi
#
#

##
## --- Pre-prepared forcing option - for now just sets up links
##
#if [ "$clm" == "prep" ] ; then
#   echo "**Setting up pre-prepared synoptic forcing from force/synoptic/$E"
#
#   DIR=$BASEDIR/force/synoptic/$E/
#
#   for i in tauewd taunwd wndspd radflx shwflx vapmix \
#      airtmp precip uwind vwind clouds relhum slp ; do
#
#      echo "|--> $i"
#
#      [ -f  $DIR/$i.a ] || tellerror "File $DIR/$i.a does not exist"
#      [ -f  $DIR/$i.b ] || tellerror "File $DIR/$i.b does not exist"
#      ln -sf $DIR/${i}.a forcing.${i}.a ||  tellerror "Could not fetch $DIR/$i.a"
#      ln -sf $DIR/${i}.b forcing.${i}.b ||  tellerror "Could not fetch $DIR/$i.b"
#   done
#fi

# ---
# --- Pre-prepared forcing option - for now just sets up links
# --- TODO: To verify, we have to go into file and check timings
# ---
echo "**Setting up pre-prepared synoptic forcing from force/synoptic/$E"
DIR=$BASEDIR/force/synoptic/$E/
#for i in tauewd taunwd wndspd radflx shwflx vapmix \
#   airtmp precip uwind vwind clouds relhum slp ; do
for i in radflx shwflx vapmix \
   airtmp precip mslprs \
   wndewd wndnwd ; do
   echo "|--> $i"
   [ -f  $DIR/$i.a ] || tellerror "File $DIR/$i.a does not exist"
   [ -f  $DIR/$i.b ] || tellerror "File $DIR/$i.b does not exist"

   # Check range of file against start and stop times
   if [ -f  $DIR/$i.a -a  -f $DIR/$i.b ] 
   then
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


#
# --- river forcing
# --- KAL: rivers are experiment-dependent
#
if [ $PRIVER -eq 1 ] ; then
   echo "**Setting up river forcing"
   cp $BASEDIR/force/rivers/$E/rivers.a forcing.rivers.a || tellerror "Could not get river .a file"
   cp $BASEDIR/force/rivers/$E/rivers.b forcing.rivers.b || tellerror "Could not get river .b file"
fi


#
# --- kpar forcing
#
if [ $JERLV -eq 0 ] ; then
   echo "**Setting up kpar forcing"
   cp $BASEDIR/force/seawifs/kpar.a forcing.kpar.a || tellerror "Could not get kpar.a file"
   cp $BASEDIR/force/seawifs/kpar.b forcing.kpar.b || tellerror "Could not get kpar.b file"
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

# Need ports.input file in these cases
if [ $tmp -eq 1 -a $LBFLAG -eq 1 ] ; then
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
      tellerror "Could not get port file port.nest in $P or  ${nestdir} for nest flow"
   fi
fi

# Need nest rmu in this case:
if [ $tmp2 -eq 1  ] ; then
   # Nest relaxation - use file in  experiment dir if present. Otherwise look in nest dir
   if [ -f $P/rmu_nest.a -a -f $P/rmu_nest.a ] ; then
      echo "Using file $P/rmu_nest.[ab] for nesting relaxation: $P/rmu_nest.[ab] -> ./rmu.[ab]"
      cp $P/rmu_nest.a rmu.a       || tellerror "Could not get port file ${P}/rmu_nest.a for nest relax"
      cp $P/rmu_nest.b rmu.b       || tellerror "Could not get port file ${P}/rmu_nest.b for nest relax"
   elif [ -f $nestdir/rmu_nest.a -a -f $nestdir/rmu_nest.a ] ; then
      echo "Using file $nestdir/rmu_nest.[ab] for nesting: $nestdir/rmu_nest.[ab] -> ./rmu.[ab]"
      cp $nestdir/rmu_nest.a rmu.a       || tellerror "Could not get port file ${nestdir}/rmu_nest.a for nest relax"
      cp $nestdir/rmu_nest.b rmu.b       || tellerror "Could not get port file ${nestdir}/rmu_nest.b for nest relax"
   else 
      tellerror "Could not get file rmu.[ab] in $P or ${nestdir} for nest relaxation"
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
# --- try to get restart input from various areas
#
#if [ $YRFLAG -ne 3 -a $init -eq 1 ] ; then
if [ $init -eq 1 ] ; then
   echo "No restart needed"
else
   filename=restart${refyear}_${iday}_${ihour}

   # Special case if CURVIINT flag set
   if [ -f $P/CURVIINT ] ; then
      if [ -f $BASEDIR/curviint/$E/${filename}.a -a $BASEDIR/curviint/$E/${filename}.b ] ; then
         echo "using CURVIINT restart files"
         for i in $BASEDIR/curviint/$E/${filename}* ; do
            echo "Using restart file $i"
            nname=$(basename $i)
            cp $i $nname
         done
         rm $P/CURVIINT
      else
         tellerror "CURVIINT set but could not find CURVIINT restart file ${filename}.[ab]"
      fi

   # Try to fetch restart from data dir $D
   elif [ -f $D/${filename}.a -o -f $D/${filename}_mem001.a ] ; then
      echo "using restart files ${filename}.[ab] from data dir $D"
      ln -sf $D/${filename}* .

   else
      tellerror "Could not find restart file ${filename}.[ab] (for init make sur yrflag < 3)"
   fi

fi




#
# --- model executable. One executable to rule them all!
#
if [ $SIGVER -eq 1 ] ; then
   TERMS=7
   MYTHFLAG=0
elif [ $SIGVER -eq 2 ] ; then
   TERMS=7
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
/bin/cp $compdir/hycom_cice  . || tellerror "Could not get hycom_cice executable"



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
   echo "No fatal errors."
#   echo "$logfile"  
#   echo "$logfile.err"
else
   echo "Some fatal errors occured. See above"
#   echo "$logfile"  
#   echo "$logfile.err"
#   echo
#   echo "Model now set up to run in $S"
fi

# Tell where stuff ended up

exit $numerr # Fails if any fatal errors occured



