abort()
{
   errcode=$?
    echo >&2 '
***************
*** ABORTED ***
***************
'

    echo "error $errorcode"
    echo "the command executing at the time of the error was"
    echo "$BASH_COMMAND"
    echo "on line ${BASH_LINENO[0]}"


    echo "An error occurred. Exiting..." >&2
    exit 1
}



array_contains () {
    local seeking=$1; shift
    local in=1
    for element; do
        if [[ $element == $seeking ]]; then
            in=0
            break
        fi
    done
    return $in
}




function blkdat_get {
   # Function retruns specified parameters from blkdat.input. Useful if you must test values
   if [ $# -ne 2 ]  ;then
      echo "Need blkdat.input file and parameter " 1>&2
      return 1
   fi

   blkdat=$1
   par=$2

   if [ ! -f $blkdat ] ; then
      echo "Unable to find $blkdat in $(pwd)" 1>&2
      return 1
   fi


   integers=("iexpt" "priver" "yrflag" "jerlv0" "sssflg" "sstflg" "relax" "vsigma" "idm" "jdm" "kdm" "nhybrd" "nsigma" "lbflag" "thflag" "iceflg" "momtyp")
   floats=("thkdf4","kapref" "sigma" "dp00" "dp00x" "dp00f" "ds00" "ds00x" "ds00f" "dp0k" "ds0k" "skmap" "nestfq" "bnstfq" "thkdf2" "baclin" "batrop" "cplifq" "visco2" "veldf2")

   param=$(printf %-6s $par)
   if    array_contains "$par" ${integers[@]}  ; then
      #value=$(grep "'$param' =" $blkdat | awk '{printf("%d", $1)}')  
      value=$(grep "'$param'" $blkdat | awk '{printf("%d", $1)}')    # Some files do not have the "=" sign
      if [ "$par" == "iexpt" ] ; then
         value=$(echo 00$value | tail -c4)
      fi
   elif  array_contains "$par" ${floats[@]}  ; then
      value=$(grep "'$param' =" $blkdat | awk '{printf("%f ", $1)}')
   else
      echo "$par is unknown" 1>&2
      return 1
   fi

   printf "blkdat_get: %-7s = $value\n" $par 1>&2
   echo $value
   return 0

}


function blkdat_get_string {
   # Function retruns specified parameters from blkdat.input. Useful if you must test values
   if [ $# -ne 3 ]  ;then
      echo "Need blkdat.input file and parameter " 1>&2
      return 1
   fi

   blkdat=$1
   par=$2
   default=$3

   if [ ! -f $blkdat ] ; then
      echo "Unable to find $blkdat in $(pwd)" 1>&2
      return 1
   fi

   strings=("nmrsti" "nmrsto" "nmarcv" "nmarcm" "nmarcs" )

   param=$(printf %-6s $par)
   if    array_contains "$par" ${strings[@]}  ; then
      value=$(egrep "'$param'" $blkdat)
      #echo "test of $param ${value}" 1>&2

      if [ ${#value} -eq 0 ] ; then
         value=$default
      else
         #echo "test2 of $param ${#value}" 1>&2
         IFS="'"
         set $value
         value=$4
         #echo "test3 of $param $value" 1>&2
         #echo "test4 of $param $4" 1>&2
      fi
   else 
      value=$default
   fi
   #echo "test5 of $param $value" 1>&2

   #echo length of $param ${#value} 1>&2
   if [ ${#value} -eq 0 ] ; then
      value=$default
   fi
   #echo "test6 of $param $value" 1>&2


   printf "blkdat_get_string: %-7s = $value\n" $par 1>&2
   echo "$value"
   return 0

}




function blkdat_pipe {
   # Function returns lines with specified parameters from blkdat.input. Useful if building a new blkdat-like file
   if [ $# -ne 2 ]  ;then
      echo "Need blkdat.input file and parameter " 1>&2
      return 1
   fi
   if [ ! -f $blkdat ] ; then
      echo "Unable to find $blkdat in $(pwd)" 1>&2
      return 1
   fi
   blkdat=$1
   par=$2
   param=$(printf %-6s $par)
   #value=$(grep "'$param' =" $blkdat)
   value=$(grep "'$param'" $blkdat)
   echo "$value"
   return 0
}


function topo_file {
   if [ $# -ne 2 ]  ;then
      echo "Need region name and topo version" 2>&1
      return 1
   fi
   R=$1
   topov=$(echo 0$2 | tail -c3)
   #echo $topov $R 1>&2
   echo "depth_${R}_${topov}"
   return 0
}


function kmt_file {
   if [ $# -ne 2 ]  ;then
      echo "Need region name and topo version" 2>&1
      return 1
   fi
   R=$1
   topov=$(echo 0$2 | tail -c3)
   #echo $topov $R 1>&2
   echo "kmt_${R}_${topov}.nc"
   return 0
}

function archv_property {
   # Get properties from archive .b  file 
   if [ $# -ne 2 ]  ;then
      echo "Need file and what to look for" 2>&1
      return 1
   fi
   
   levels=()
   if [ "$2" == "kdm" ] ; then 
      while read fieldname tmp nstep dtime k dens range ; do 
         #echo $k
         levels+=($k)
      done <<<"$( tail -n +11  $1)"
      #echo levels=${levels[@]}
      #echo $levels
      IFS=$'\n' sortedlevels=($(sort -n <<<"${levels[*]}"))
      unset IFS
      #echo sortedlevels=${sortedlevels[@]}
      printf %s "${sortedlevels[@]: -1}"
   else 
      echo "No method for $2 ... " 1>&2
   fi

}

function copy_grid_files {
   if [ $# -ne 1 ]  ;then
      echo "Need target location and what to look for" 2>&1
      return 1
   fi
   TARGETDIR=$1
   [ ! -d $TARGETDIR ] && mkdir -p $TARGETDIR
   cd $TARGETDIR
   touch regional.grid.a regional.grid.b
   rm    regional.grid.a regional.grid.b
   cp ${BASEDIR}/topo/regional.grid.b regional.grid.b     || { echo "Could not get regional.grid.b file " ; exit 1 ; }
   cp ${BASEDIR}/topo/regional.grid.a regional.grid.a     || { echo "Could not get regional.grid.a file " ; exit 1 ; }
}



function copy_topo_files {
   if [ $# -ne 1 ]  ;then
      echo "Need target location and what to look for" 2>&1
      return 1
   fi
   TARGETDIR=$1
   [ ! -d $TARGETDIR ] && mkdir -p $TARGETDIR
   cd $TARGETDIR
   touch regional.depth.a regional.depth.b
   rm    regional.depth.a regional.depth.b
   cp ${BASEDIR}/topo/depth_${R}_${T}.a regional.depth.a  || { echo "Could not get regional.depth.a file " ; exit 1 ; }
   cp ${BASEDIR}/topo/depth_${R}_${T}.b regional.depth.b  || { echo "Could not get regional.depth.b file " ; exit 1 ; }
}


function copy_setup_files {
   if [ $# -ne 1 ]  ;then
      echo "common_functions.sh - copy_setup_files : Need target location and what to look for" 2>&2
      return 1
   fi
   TARGETDIR=$1
   [ ! -d $TARGETDIR ] && mkdir -p $TARGETDIR
   cd $TARGETDIR
   touch infile.in blkdat.input
   rm    infile.in blkdat.input
   copy_grid_files $S
   copy_topo_files $S
   cp ${BASEDIR}/expt_${X}/blkdat.input blkdat.input      || { echo "Could not get file blkdat.input " ; exit 1 ; }
   #cp ${BASEDIR}/expt_${X}/infile.in infile.in            || { echo "Could not get file infile.in " ; exit 1 ; }
}


function gmap_file {
   if [ $# -ne 1 ]  ;then
      echo "$(basename $0) - $FUNCNAME : Need region name as input" 1>&2
      return 1
   fi
   echo "${1}.gmap"
}


function source_dir {
   if [ $# -ne 3 ]  ;then
      echo "$(basename $0) - $FUNCNAME : Need version, number of terms and thflag as input" 1>&2
      return 1
   fi
   V=$1
   TERMS=$2
   THFLAG=$3
   TERMS2=$(echo 0$TERMS | tail -c3)
   # Set up rel path
   bd=src_${V}ZA-${TERMS2}Tsig${THFLAG}-i-sm-sse_relo_mpi/
   echo $bd
   return 0
}

#function parse_isotime {
#   if [[ $1 =~  ([0-9]){4}-([0-9]){2}-([0-9]){2}T([0-9]){2}:([0-9]){2}:([0-9]){2} ]] ; then
#      yr=${BASH_REMATCH[1]}
#      mo=${BASH_REMATCH[2]}
#      dm=${BASH_REMATCH[3]}
#      hr=${BASH_REMATCH[4]}
#      min=${BASH_REMATCH[5]}
#      sec=${BASH_REMATCH[6]}
#      echo "$yr $mo $dm $hr $min $sec"
#      return 0
#   else 
#      echo "end time not in righ format" ; 
#      return 1
#   fi
#}


tellerror () {
  echo "[FATAL  ] $1" 
  let numerr=$numerr+1; 
  #echo "[FATAL  ] $1" >> $logfile
}

tellwarn () {
  echo "[WARNING] $1" 
  #echo "[WARNING] $1" >> $logfile
}
