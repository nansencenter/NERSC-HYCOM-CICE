#!/bin/bash
export STARTDIR=$(pwd)
export BASEDIR=$(cd $(dirname $0)/.. && pwd)
source $BASEDIR/bin/common_functions.sh
export SCRATCH=$BASEDIR/subregion/SCRATCH
[ ! -d $SCRATCH ] && mkdir $SCRATCH

usage="
   Example:
      $(basename $0) expt  new_region_path new_region_expt 

   Example:
      $(basename $0) 90.9 /work/$USER/hycom/TP4a0.12/ 03.1
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


# Experiment  needs experiment number
if [ $# -ne 3 ] ; then
   echo -e "$usage"
   exit
fi

# Explore provided input path to get target region 
thisexpt=$1
newregionpath=$2
newregionexpt=$3
source ${newregionpath}/REGION.src || { echo "Could not source ${newregionpath}/REGION.src" ; exit 1 ; }
source ${newregionpath}/expt_${newregionexpt}/EXPT.src || { echo "Could not source ${newregionpath}/expt_${newregionexpt}/EXPT.src" ; exit 1 ; }
NR=$R
NX=$X
NE=$E
NT=$T
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$thisexpt/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$thisexpt/EXPT.src" ; exit 1 ; }
echo "This region name    :$R"
echo "This experiment     :$X"
echo "This experiment topo:$T"
echo "New  region name    :$NR"
echo "New  experiment     :$NX"
echo "New  experiment topo:$NT"

# Get regional.grid of target region
target_grid=$NR.regional.grid
tmp=$newregionpath/topo/regional.grid
cp ${tmp}.a $SCRATCH/${target_grid}.a || { echo "Could not get ${tmp}.a file " ; exit 1 ; }
cp ${tmp}.b $SCRATCH/${target_grid}.b || { echo "Could not get ${tmp}.b file " ; exit 1 ; }

# Get topo file of target
#target_topo=$(topo_file $NR $NT)
#tmp=$newregionpath/topo/$target_topo
#cp $tmp.a $SCRATCH/${target_topo}.a || { echo "Could not get $tmp.a " ; exit 1 ; }
#cp $tmp.b $SCRATCH/${target_topo}.b || { echo "Could not get $tmp.b " ; exit 1 ; }
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



logfile=$SCRATCH/isuba_topog.log
touch $logfile && rm $logfile
echo "Running isuba_topog. Log can be found in $logfile"

cat <<EOF
${target_gmap}.a
${source_topo}.a
${target_topo}.a
test
  ${target_idm}    'idm   '                                                               
  ${target_jdm}    'jdm   '    
EOF

/home/nersc/knutali/Models/hycom/HYCOM_ALL_2.2.72/ALL/subregion/src/isuba_topog >> $logfile <<EOF
${target_gmap}.a
${source_topo}.a
${target_topo}.a
test
  ${target_idm}    'idm   '                                                               
  ${target_jdm}    'jdm   '    
EOF


echo "Normal exit"
exit 0
