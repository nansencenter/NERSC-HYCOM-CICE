#!/bin/bash
export STARTDIR=$(pwd)
export BASEDIR=$(cd $(dirname $0)/.. && pwd)
source $BASEDIR/bin/common_functions.sh
export SCRATCH=$BASEDIR/subregion/SCRATCH
[ ! -d $SCRATCH ] && mkdir $SCRATCH

iscan=15
usage="
   Example:
      $(basename $0) [-s iscan] expt  new_region_path new_region_expt archive1 archive2 ....

   Example:
      $(basename $0) 90.9 /work/$USER/hycom/TP4a0.12/ 03.1 archive1 
      $(basename $0) -s $iscan 90.9 /work/$USER/hycom/TP4a0.12/ 03.1 archive1 

   Optional argument 'iscan' has default value of 15
"
options=$(getopt -o s:  -- "$@")
[ $? -lt 4 ] || {
    echo "Incorrect options provided"
    exit 1
}
eval set -- "$options"
while true; do
    case "$1" in
    -s)
       shift;
       iscan=$1
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done


# Experiment  needs experiment number
if [ $# -lt 4 ] ; then
   echo -e "$usage"
   exit
fi

# Explore provided input path to get target region 
thisexpt=$1
newregionpath=$2
newregionexpt=$3
shift 3
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
target_topo=$(topo_file $NR $NT)
tmp=$newregionpath/topo/$target_topo
cp $tmp.a $SCRATCH/${target_topo}.a || { echo "Could not get $tmp.a " ; exit 1 ; }
cp $tmp.b $SCRATCH/${target_topo}.b || { echo "Could not get $tmp.b " ; exit 1 ; }

# Get topo file of source
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
#echo blkdat_get ${target_grid}.b 
target_idm=$(blkdat_get ${target_grid}.b idm)
target_jdm=$(blkdat_get ${target_grid}.b jdm)

target_dir=../${R}_${E}/${NR}_${NE}/
mkdir -p $target_dir || { echo "Could not create ${trget_dir}" ; exit 1 ; }
target_dir=$(cd ${target_dir} && pwd)/
echo "Target dir = ${target_dir}"

logfile=$SCRATCH/isubaregion.log
touch $logfile && rm $logfile
for source_archv in $@ ; do
   echo "***"

   my_source_archv=$(echo $STARTDIR/$source_archv |  sed "s/\.[ab]$//")
   echo "Processing $my_source_archv"

   if [ ! -f ${my_source_archv}.a -o ! -f ${my_source_archv}.b ] ; then
      echo "Source file ${my_source_archv}.[ab] does not exist"
      exit 1
   fi
   target_archv=$(basename $my_source_archv)
   touch ${target_archv}.a && rm ${target_archv}.a
   touch ${target_archv}.b && rm ${target_archv}.b


#c --- 'flnm_reg'  = target sub-region grid       filename
#c --- 'flnm_map'  = target sub-region grid map   filename
#c --- 'flnm_top'  = target bathymetry filename, or 'NONE'
#c --- 'flnm_tin'  = input  bathymetry filename, or 'NONE'
#c --- 'flnm_in'   = input  archive    filename
#c --- 'flnm_out'  = output archive    filename
#c --- 'cline_out' = output title line (replaces preambl(5))
#c --- 'idm   ' = longitudinal array size
#c --- 'jdm   ' = latitudinal  array size
#c --- 'smooth' = smooth interface depths (0=F,1=T)
#c --- 'iscan ' = How many grid points to search for landfill process
echo "Running isubaregion. Log can be found in $logfile"
prog=${HYCOM_ALL}/subregion/src/isubaregion
#cat << EOF # For diag
${prog} >> $logfile <<EOF
${target_grid}.a
${target_gmap}.a
${target_topo}.a
${source_topo}.a
${my_source_archv}.a
${target_archv}.a
test
  ${target_idm}    'idm   '                                                               
  ${target_jdm}    'jdm   '    
    0    'iceflg'    
    0    'smooth'    
   $iscan   'iscan '    
EOF


   if [ ! -f ${target_archv}.a -o ! -f ${target_archv}.b ] ; then
      echo "Error : New file ${target_archv}.[ab] does not exist"
      echo "This is probably due to errors during extrapolation. Check log file $logfile"
      exit 1
   else 
      echo "Moving ${target_archv}.[ab] to $BASEDIR/subregion/${target_dir}"
      mv ${target_archv}.* ${target_dir}
   fi
   

done


echo "Normal exit"
exit 0
