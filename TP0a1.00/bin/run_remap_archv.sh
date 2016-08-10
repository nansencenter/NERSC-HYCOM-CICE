# Remapping of archive file to new vertical grid
usage="
   Example:
      $(basename $0) expt  archive1 archive2 ....


   Example:
      $(basename $0) 90.9  archive1 archive2 ....
"
options=$(getopt -o m:  -- "$@")
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
if [ $# -lt 2 ] ; then
    echo "Incorrect options provided"
    echo "$usage"
    exit 1
fi

export STARTDIR=$(pwd)
export BASEDIR=$(cd $(dirname $0)/.. && pwd)
source $BASEDIR/bin/common_functions.sh
export SCRATCH=$BASEDIR/nest/SCRATCH
[ ! -d $SCRATCH ] && mkdir $SCRATCH


thisexpt=$1
shift 
source ${BASEDIR}/REGION.src || { echo "Could not source ${BASEDIR}/REGION.src" ; exit 1 ; }
source ${BASEDIR}/expt_$thisexpt/EXPT.src || { echo "Could not source ${BASEDIR}/expt_$thisexpt/EXPT.src" ; exit 1 ; }
echo "This region name    :$R"
echo "This experiment     :$X"
echo "This experiment topo:$T"
TARGETDIR=$BASEDIR/nest/$E/
[ ! -d ${TARGETDIR} ] && mkdir -p ${TARGETDIR}

# Get regional.grid 
tmp=$BASEDIR/topo/regional.grid
cp ${tmp}.a $SCRATCH/ || { echo "Could not get ${tmp}.a file " ; exit 1 ; }
cp ${tmp}.b $SCRATCH/ || { echo "Could not get ${tmp}.b file " ; exit 1 ; }

# Get topo file of target
target_topo=$(topo_file $R $T)
tmp=$BASEDIR/topo/$target_topo
cp $tmp.a $SCRATCH/regional.depth.a || { echo "Could not get $tmp.a " ; exit 1 ; }
cp $tmp.b $SCRATCH/regional.depth.b || { echo "Could not get $tmp.b " ; exit 1 ; }

# Copy blkdat.input to SCRATCH
cp $BASEDIR/expt_${X}/blkdat.input $SCRATCH || { echo "could not get $BASEDIR/expt_${X}/blkdat.input" ; exit 1 ; }


# Get dimensions from target region
cd $SCRATCH || { echo "Could not cd to $SCRATCH" ; exit 1 ; }
echo blkdat_get ${target_grid}.b 
target_idm=$(blkdat_get regional.grid.b idm)
target_jdm=$(blkdat_get regional.grid.b jdm)
echo "New idm          :$target_idm"
echo "New jdm          :$target_jdm"


for source_archv in $@ ; do

   my_source_archv=$(echo $STARTDIR/$source_archv |  sed "s/\.[ab]$//")
   echo "Processing $my_source_archv"

   #kdm old from archive file (largest "k")
   kdmold=$(archv_property ${my_source_archv}.b kdm)
   echo "Old kdm          :$kdmold"

   target_archv=$(basename ${my_source_archv})
   [ -f ${target_archv}.a ] && rm ${target_archv}.a
   [ -f ${target_archv}.b ] && rm ${target_archv}.b

## Get from target blkdat.input. NB: some of these can be empty at this stage
#iexpt=$(blkdat_get blkdat.input iexpt)
#yrflag=$(blkdat_get blkdat.input yrflag)
#idm=$(blkdat_get blkdat.input idm)
#jdm=$(blkdat_get blkdat.input jdm)
#nhybrd=$(blkdat_get blkdat.input nhybrd)
#nsigma=$(blkdat_get blkdat.input nsigma)
#dp00x=$(blkdat_get blkdat.input dp00x)
#dp00f=$(blkdat_get blkdat.input dp00f)
#ds00=$(blkdat_get blkdat.input ds00)
#ds00x=$(blkdat_get blkdat.input ds00x)
#ds00f=$(blkdat_get blkdat.input ds00f)
#dp0k=$(blkdat_get blkdat.input dp0k)
#ds0k=$(blkdat_get blkdat.input dp00)
#skmap=$(blkdat_get blkdat.input skmap)
#sigma=$(blkdat_get blkdat.input sigma)

   dp00=$(blkdat_get blkdat.input dp00)
   kdmnew=$(blkdat_get blkdat.input kdm)



#c --- 'thbase' = new reference density (sigma units)

prog=${HYCOM_ALL}/archive/src/remap_archv
logfile="$SCRATCH/remap_archv_${target_archv}.log"
echo "Running $prog - logfile is $logfile"
if [ -n "$dp00" ] ; then
#c --- 'flnm_i' = name of original archive file
#c --- 'flnm_o' = name of target   archive file
#c --- 'iexpt ' = experiment number x10  (000=from archive file)
#c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
#c --- 'idm   ' = longitudinal array size
#c --- 'jdm   ' = latitudinal  array size
#c --- 'kdmold' = original number of layers
#c --- 'kdmnew' = target   number of layers
#c --- 'nhybrd' = new number of hybrid levels (0=all isopycnal)
#c --- 'nsigma' = new number of sigma  levels (nhybrd-nsigma z-levels)
#c --- 'dp00'   = new deep    z-level spacing minimum thickness (m)
#c --- 'dp00x'  = new deep    z-level spacing maximum thickness (m)
#c --- 'dp00f'  = new deep    z-level spacing stretching factor (1.0=const.z)
#c --- 'ds00'   = new shallow z-level spacing minimum thickness (m)
#c --- 'ds00x'  = new shallow z-level spacing maximum thickness (m)
#c --- 'ds00f'  = new shallow z-level spacing stretching factor (1.0=const.z)
$prog<<EOF > $logfile 2>&1
${my_source_archv}.a
${target_archv}.a
$(blkdat_pipe blkdat.input iexpt)
$(blkdat_pipe blkdat.input yrflag)
$(blkdat_pipe blkdat.input idm)
$(blkdat_pipe blkdat.input jdm)
  $kdmold      'kdmold' = original   number of layers
  $kdmnew      'kdmnew' = target   number of layers
$(blkdat_pipe blkdat.input nhybrd)
$(blkdat_pipe blkdat.input nsigma)
$(blkdat_pipe blkdat.input dp00)
$(blkdat_pipe blkdat.input dp00x)
$(blkdat_pipe blkdat.input dp00f)
$(blkdat_pipe blkdat.input ds00)
$(blkdat_pipe blkdat.input ds00x)
$(blkdat_pipe blkdat.input ds00f)
$(blkdat_pipe blkdat.input thbase)
$(blkdat_pipe blkdat.input sigma)
EOF
else
#c --- 'flnm_i' = name of original archive file
#c --- 'flnm_o' = name of target   archive file
#c --- 'iexpt ' = experiment number x10  (000=from archive file)
#c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
#c --- 'idm   ' = longitudinal array size
#c --- 'jdm   ' = latitudinal  array size
#c --- 'kdmold' = original number of layers
#c --- 'kdmnew' = target   number of layers
#c --- 'nhybrd' = new number of hybrid levels (0=all isopycnal)
#c --- 'nsigma' = new number of sigma  levels (nhybrd-nsigma z-levels)
#c --- 'dp0k'   = new deep    z-level spacing minimum thickness (m)
#c --- 'ds0k '  = new deep    z-level spacing maximum thickness (m)
echo
$prog<<EOF > $logfile 2>&1
${my_source_archv}.a
${target_archv}.a
$(blkdat_pipe blkdat.input iexpt)
$(blkdat_pipe blkdat.input yrflag)
$(blkdat_pipe blkdat.input idm)
$(blkdat_pipe blkdat.input jdm)
  $kdmold      'kdmold' = original   number of layers
  $kdmnew      'kdmnew' = target   number of layers
$(blkdat_pipe blkdat.input nhybrd)
$(blkdat_pipe blkdat.input nsigma)
$(blkdat_pipe blkdat.input dp0k)
$(blkdat_pipe blkdat.input ds0k)
$(blkdat_pipe blkdat.input thbase)
$(blkdat_pipe blkdat.input sigma)
EOF
fi

if [ -f ${target_archv}.a -a -f ${target_archv}.b ] ; then
   echo "Found ${target_archv}.[ab] - moving to $TARGETDIR"
   cp ${target_archv}.* ${TARGETDIR}
else 
   echo "${target_archv}.[ab] not found - aborting. See $logfile for diag"
   exit 1
fi

echo
done 

echo "Normal exit"
exit 0
