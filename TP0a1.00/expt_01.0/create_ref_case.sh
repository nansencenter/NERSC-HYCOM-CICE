HERE=$(pwd)
X=01.0

echo "Logs can be found in $HERE/log"

# Create hycom_all exacutables
cd $HERE
cd ../../hycom_ALL/hycom_2.2.72_ALL/
echo "compiling hycom_all"
csh Make_all.com ARCH=amd64 > $HERE/log/ref_hycom_all.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

# Create hycom exacutable
cd $HERE
cd ../../hycom/HYCOM_2.2.98_ESMF5/RELO/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi/
source set_env.sh
echo "compiling hycom"
csh Make_cice.csh > $HERE/log/ref_hycom.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

# Create z-level relaxation files
cd $HERE
echo "z climatology"
../bin/z_generic.sh 01.0 woa2013 > $HERE/log/ref_z_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

# Create relaxation files on hybrid coordinates
cd $HERE
echo "hybrid climatology"
../bin/relaxi.sh 01.0 woa2013 2 > $HERE/log/ref_hybrid_relax.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

# Create relaxation mask
echo "relaxation mask"
cat <<EOF | ../bin/relax_rmu.sh 01.0  > $HERE/log/ref_rmu_mask.out 2>&1
F
F
F
T
F
20
20
EOF
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

# Create simple river forcing
cd $HERE
echo "river forcing"
../bin/river_nersc.sh 01.0 100 300 > $HERE/log/ref_river_nersc.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

# Create kpar file
echo "kpar forcing"
../bin/seawifs_mon_kpar.sh > $HERE/log/ref_seawifs.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

# Create tiling 
echo "grid tiling"
../bin/tile_grid.sh -2 -2 9.5 01 > $HERE/log/ref_tiling.out 2>&1
res=$?
[ $res -eq 0 ] && echo "Success"
[ $res -ne 0 ] && echo "Failure..."

echo "If things went fine, you can now generate test forcing like this: "
echo "    ../bin/atmo_synoptic_new.sh 01.0 erai 2013-01-02T00:00:00  2013-01-05T00:00:00"
echo "Then edit the job script pbsjob.sh to read and make sure expt_preprocess.sh is called like this"
echo "    ../bin/preprocess_expt.sh 2015-01-02T00:00:00 2015-01-05T00:00:00 --init    ... "

