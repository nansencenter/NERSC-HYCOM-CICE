#!/bin/bash

myclim="woa2018" # Climatology to use
#myclim="phc" # Climatology to use

Icore=23
Jcore=23

# Must be in expt dir to run this script
if [ -f EXPT.src ] ; then
   export BASEDIR=$(cd .. && pwd)
else
   echo "Could not find EXPT.src. This script must be run in expt dir"
   exit 1
fi
EDIR=$(pwd)/                           # Location of this script
BASEDIR=$(cd $(dirname $0)/.. && pwd)/ # Location of basedir
source $BASEDIR/REGION.src
source $EDIR/EXPT.src

KSIGMA=$(egrep "'thflag'"  $EDIR/blkdat.input  | sed "s/.thflag.*$//" | tr -d "[:blank:]")
export NTRACR=`grep "'ntracr' =" $EDIR/blkdat.input | awk '{printf("%03d", $1)}'`

# Create z-level relaxation files
cd $EDIR
echo "z climatology"
z_generic.sh $myclim > $EDIR/log/ref_z_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_z_relax.out"
echo ".."

# Create relaxation files on hybrid coordinates
cd $EDIR
echo "hybrid climatology"
relaxi.sh $myclim   > $EDIR/log/ref_hybrid_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_hybrid_relax.out"
echo ".."

if [ $NTRACR -ne 0 ] ; then
  # Create biology z-relaxation files from WOA2013
  cd $EDIR
  echo "bio relax climatology"
  z_woa2013_bio.sh $KSIGMA > $EDIR/log/ref_bio_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_bio_relax.out"
  echo ".."

  # Create biology hybrid-relaxation files from WOA2013
  cd $EDIR
  echo "bio hybrid relax climatology"
  relax_sil.sh ${X} woa2013 > $EDIR/log/ref_sil_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_sil_relax.out"

  relax_pho.sh ${X} woa2013 > $EDIR/log/ref_pho_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_pho_relax.out"

  relax_nit.sh ${X} woa2013 > $EDIR/log/ref_no3_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_no3_relax.out"

  relax_oxy.sh ${X} woa2013 > $EDIR/log/ref_oxy_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_oxy_relax.out"

  # Create CO2 z-relaxation files from GLODAPV2
  cd $EDIR
  echo "co2 relax climatology"
  z_glodap_co2.sh $KSIGMA > $EDIR/log/ref_co2_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_co2_relax.out"
  echo ".."

  # Create CO2 hybrid-relaxation files from GLODAPV2
  cd $EDIR
  echo "co2 hybrid relax climatology"
  relax_alk.sh ${X} woa2013 > $EDIR/log/ref_alk_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_alk_relax.out"

  relax_dic.sh ${X} woa2013 > $EDIR/log/ref_dic_relax.out 2>&1
  res=$?
  [ $res -eq 0 ] && echo "Success"
  [ $res -ne 0 ] && echo "Failure... Log in $EDIR/log/ref_pho_relax.out"

  echo ".."
fi

# Create relaxation mask
echo "relaxation mask"
cat <<EOF | relax_rmu.sh ${X}  > $EDIR/log/ref_rmu_mask.out 2>&1
F
F
F
T
T
20
20
EOF
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."
echo ".."

# Create simple river forcing
cd $EDIR
echo "river forcing, if biology active, may take some time"
if [ $NTRACR -ne 0 ] ; then
   river_nersc.sh 100 300 $INPUTDIR/rivers_ahype-ehype_clim_rev2.dat $INPUTDIR/biorivers.dat > $EDIR/log/ref_river_nersc.out 2>&1
else
   river_nersc.sh 100 300 $INPUTDIR/rivers_ahype-ehype_clim_rev2.dat > $EDIR/log/ref_river_nersc.out 2>&1
fi

res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_river_nersc.out"
echo ".."

# Create kpar file
echo "kpar forcing"
seawifs_mon_kpar.sh > $EDIR/log/ref_seawifs.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_seawifs.out"
echo ".."

# Create tiling. 
echo "grid tiling"
tile_grid.sh -${Icore} -${Jcore} ${T} > $EDIR/log/ref_tiling.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure...  Log in  $EDIR/log/ref_tiling.out" 
echo ".."

echo "If things went fine, you can now generate test forcing like this: "
echo "    atmo_synoptic.sh erai 2001-01-01T00:00:00  2001-01-05T00:00:00"
echo "Then edit the job script pbsjob.sh to read and make sure expt_preprocess.sh is called like this"
echo "    expt_preprocess.sh 2001-01-01T00:00:00 2001-01-05T00:00:00 --init "
