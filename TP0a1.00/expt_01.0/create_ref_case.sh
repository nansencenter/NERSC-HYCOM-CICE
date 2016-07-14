HERE=$(pwd)
X=01.0

echo "Logs can be found in $HERE/log"

# Create hycom_all exacutables
cd $HERE
cd ../../hycom_ALL/hycom_2.2.72_ALL/
echo "compiling hycom_all"
csh Make_all.com ARCH=amd64 > $HERE/log/ref_hycom_all.out 2>&1

# Create hycom exacutable
cd $HERE
cd ../../hycom/HYCOM_2.2.98_ESMF5/RELO/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi/
source set_env.sh
echo "compiling hycom"
csh Make_cice.csh > $HERE/log/ref_hycom.out 2>&1

# Create z-level relaxation files
cd $HERE
echo "z climatology"
../bin/z_generic.sh 01.0 woa2013 > $HERE/log/ref_z_relax.out 2>&1

# Create relaxation files on hybrid coordinates
cd $HERE
echo "hybrid climatology"
../bin/relaxi.sh 01.0 woa2013 2 > $HERE/log/ref_hybrid_relax.out 2>&1

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

# Create simple river forcing
cd $HERE
echo "river forcing"
../bin/river_nersc.sh 01.0 100 300 > $HERE/log/ref_river_nersc.out 2>&1

# Create kpar file
echo "kpar forcing"
../bin/seawifs_mon_kpar.sh > $HERE/log/ref_seawifs.out 2>&1

