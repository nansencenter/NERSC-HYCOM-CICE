#!/bin/bash
export BASEDIR=$(cd $(dirname $0)/.. && pwd)
export SCRATCH=$BASEDIR/subregion/SCRATCH
[ ! -d $SCRATCH ] && mkdir $SCRATCH

usage="
   This script will set up mapping from this hycom region to the region in the input dir.
   The input dir contains a hycom \"tree\" with forcing, experimetns, etc etc.
   The script uses isuba_gmapi, which uses a conjugate gradient method to find the
   mapping. From source code:

   c     create an array index map to a diferent-grid subregion from a 
   c     full region hycom file.
   c
   c     subregion grid is arbitrary, except that any part of this grid
   c     that is outside the input grid must be land.
   c
   c     Alan J. Wallcraft,  NRL,  January 2009.
   c

   Usage  :
      $(basename $0) [-m n] [path to hycom region]

   
   Example:
      $(basename $0) /work/$USER/hycom/TP4a0.12/

   Optional argument m n denotes maxinc parameter sent to isuba_gmapi. From source code:
   c     'maxinc' is the user provided maximum input array index jump 
   c     on the target grid, i.e. abs(x_out(i,j-1)-x_out(i,j)) and
   c     abs(y_out(i,j-1)-y_out(i,j)) must both be less than maxinc (after
   c     allowing for periodic and arctic boundaries) for all i and j.
   c     The total run time is proportional to maxinc**2, but if maxinc
   c     is set too small the results will be inaccurate.
"
options=$(getopt -o m:  -- "$@")
[ $? -eq 0 ] || {
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
if [ $# -ne 1 ] ; then
   echo -e "$usage"
   exit
fi


# Explore provided input path to get target region 
newregionpath=$1
source ${newregionpath}/REGION.src || { echo "Could not source ${newregionpath}/REGION.src" ; exit 1 ; }
NR=$R
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
echo "This region name :$R"
echo "New  region name :$NR"

#TODO: Get regional.grid of source region
cp $BASEDIR/topo/regional.grid.a $SCRATCH || { echo "Could not get regional.grid.a file " ; exit 1 ; }
cp $BASEDIR/topo/regional.grid.b $SCRATCH || { echo "Could not get regional.grid.b file " ; exit 1 ; }

#TODO: Get regional.grid of destination region
echo "New  region files from  :$newregionpath/topo/"
newgrid=${NR}.regional.grid
cp $newregionpath/topo/regional.grid.a $SCRATCH/$newgrid.a || { echo "Could not get $newgrid.a file " ; exit 1 ; }
cp $newregionpath/topo/regional.grid.b $SCRATCH/$newgrid.b || { echo "Could not get $newgrid.b file " ; exit 1 ; }

# Get dimensions from target region
cd $SCRATCH || { echo "Could not cd to $SCRATCH" ; exit 1 ; }
idm=$(blkdat_get regional.grid.b idm)
jdm=$(blkdat_get regional.grid.b jdm)
echo "New idm          :$idm"
echo "New jdm          :$jdm"


newmap=${NR}.gmap
touch ${newmap}.a && rm ${newmap}.a
touch ${newmap}.b && rm ${newmap}.b
prog=/home/nersc/knutali/Models/hycom/HYCOM_ALL_2.2.72/ALL/subregion/src_2.2.12/isuba_gmapi
$prog << EOF 
${newgrid}.a
${newmap}.a
Mapping from region $R to region $NR
  $idm   'idm   '
  $jdm   'jdm   '
  50    'maxinc'
EOF

# Not exhaustive
pwd
if [ -s ${newmap}.a ] ; then
   cp ${newmap}.a ..
else 
   { echo "Error - File $newmap.a is empty " ; exit 1 ; }
fi
if [ -s ${newmap}.b ] ; then
   cp ${newmap}.b ..
else 
   { echo "Error - File $newmap.b is empty " ; exit 1 ; }
fi

echo "Normal exit. ${newmap}.[ab] placed in $(cd .. && pwd)"

