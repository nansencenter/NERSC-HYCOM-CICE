#!/bin/bash

usage="

   This routine must run in a experiment directory.

   The script uses an existing mapping to map bathymetry used by this experiment to the 
   bathymetry used by the experiment in the directory specified

   Usage:
      $(basename $0) new_expt

   Example:
      $(basename $0) /work/$USER/hycom/TP4a0.12/expt_03.1
"
options=$(getopt -o m:  -- "$@")
[ $? -lt 4 ] || {
    echo "Incorrect options provided"
    exit 1
}
maxinc=50
eval set -- "$options"
while true; do
    case "$1" in
    -m)
       shift;
       maxinc=$1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done


export EDIR=$(pwd)
export BINDIR=$(cd $(dirname $0) && pwd)
export BASEDIR=$(cd .. && pwd)
source $BINDIR/common_functions.sh

# Experiment  needs experiment number
if [ $# -ne 1 ] ; then
   echo -e "$usage"
   exit
fi

# Explore provided input path to get target region 
newexptpath=$1
newregionpath=$(dirname $newexptpath)
thisexpt=$1

source ${newregionpath}/REGION.src || { echo "Could not source ${newregionpath}/REGION.src" ; exit 1 ; }
source ${newexptpath}/EXPT.src || { echo "Could not source ${newexptpath}/EXPT.src" ; exit 1 ; }
NR=$R
NX=$X
NE=$E
NT=$T
echo "New  region name    :$NR"
echo "New  experiment     :$NX"
echo "New  experiment topo:$NT"
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${EDIR}/EXPT.src || { echo "Could not source ${EDIR}/EXPT.src" ; exit 1 ; }
echo "This region name    :$R"
echo "This experiment     :$X"
echo "This experiment topo:$T"

target_dir=$BASEDIR/subregion/${E}/${NR}_${NE}/
[[ ! -d ${target_dir} ]] && mkdir -p ${target_dir}
export SCRATCH=${target_dir}/SCRATCH
[ ! -d $SCRATCH ] && mkdir -p $SCRATCH


# Get regional.grid of expt region
source_grid=regional.grid
tmp=$BASEDIR/topo/regional.grid
cp ${tmp}.a $SCRATCH/${source_grid}.a || { echo "Could not get ${tmp}.a file " ; exit 1 ; }
cp ${tmp}.b $SCRATCH/${source_grid}.b || { echo "Could not get ${tmp}.b file " ; exit 1 ; }


# Get regional.grid of target region
target_grid=$NR.regional.grid
tmp=$newregionpath/topo/regional.grid
cp ${tmp}.a $SCRATCH/${target_grid}.a || { echo "Could not get ${tmp}.a file " ; exit 1 ; }
cp ${tmp}.b $SCRATCH/${target_grid}.b || { echo "Could not get ${tmp}.b file " ; exit 1 ; }

# Clear topo file of target
target_topo=tmp
touch $SCRATCH/${target_topo}.a && rm $SCRATCH/${target_topo}.a 
touch $SCRATCH/${target_topo}.b && rm $SCRATCH/${target_topo}.b

## Get topo file of source
source_topo=$(topo_file $R $T)
tmp=$BASEDIR/topo/$source_topo
cp $tmp.a $SCRATCH/${source_topo}.a || { echo "Could not get $tmp.a " ; exit 1 ; }
cp $tmp.b $SCRATCH/${source_topo}.b || { echo "Could not get $tmp.b " ; exit 1 ; }

# Get gmap file
target_gmap=$NR.gmap
tmp=$BASEDIR/subregion/${target_gmap}
if [ ! -f $tmp.a -o ! -f  $BASEDIR/$tmp.b ] ; then
   cp $tmp.a $SCRATCH/${target_gmap}.a
   cp $tmp.b $SCRATCH/${target_gmap}.b
else 
   echo "Could not get gmap files $tmp.a and  $tmp.b " 
   exit 1
fi

# Get dimensions from target region
cd $SCRATCH || { echo "Could not cd to $SCRATCH" ; exit 1 ; }
echo blkdat_get ${target_grid}.b 
target_idm=$(blkdat_get ${target_grid}.b idm)
target_jdm=$(blkdat_get ${target_grid}.b jdm)



prog=${HYCOM_ALL}/subregion/src/isuba_topog 
logfile=$SCRATCH/isuba_topog.log
touch $logfile && rm $logfile

#cat <<EOF
#${target_gmap}.a
#${source_topo}.a
#${target_topo}.a
#test
#  ${target_idm}    'idm   '                                                               
#  ${target_jdm}    'jdm   '    
#EOF


echo
echo
echo "Runing bathymetry interpolation routine $prog in directory $SCRATCH"
echo "Log can be found in $logfile"
$prog >> $logfile <<EOF
${target_gmap}.a
${source_topo}.a
${target_topo}.a
test
  ${target_idm}    'idm   '                                                               
  ${target_jdm}    'jdm   '    
EOF

cp ${target_topo}.a  ${target_dir}/${source_topo}.a
cp ${target_topo}.b  ${target_dir}/${source_topo}.b

echo "New topography file for region in $NR in ${target_dir}/${source_topo}.[ab]"
echo "Note that topo filename still uses 'old' region and topo version but is interpolated to the new region $NR"

echo "Normal exit"
exit 0
