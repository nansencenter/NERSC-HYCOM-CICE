# setup.sh: get regional.grid.[ab] files and find pivot points

# data name (should match $FA, $FB and $PP variables below)
dname="wamnsea"


# name of outputs from script 
# (these will be renamed depending on the model)
# (NB: these should be the same as 'filenamea' and 'filenameb' in ../p_wamnsea.F90):
FA="pivots_"${dname}".a"
FB="pivots_"${dname}".b"

# location of ModelInput directory
# (has the regional.grid.a files for all the models)
# (among other things)
M_INP="/work/shared/nersc/msc/ModelInput"

# Model dependent stuff:
MODEL="BS1a0.045"
topo_dir=${M_INP}/Barents_Hyc2.2.12/topo/
# MODEL="FR1a0.03"
# topo_dir=${M_INP}/FramStrait_Hyc2.2.12/${MODEL}/topo
# MODEL="TP4a0.12"
# topo_dir=${M_INP}/TOPAZ4/${MODEL}/topo


ln -sf ${topo_dir}/regional.grid.a .
ln -sf ${topo_dir}/regional.grid.b .
ln -sf ${topo_dir}/regional.depth.a .
ln -sf ${topo_dir}/regional.depth.b .
ln -sf /work/shared/nersc/msc/WAMNSEA/wam_nsea.an.20120228.nc .

# name of script to find pivot points:
PP="../"${dname}"-piv-hyc2.2.12"
echo " "
echo "Getting pivot points..."
echo ${PP}
echo " "
${PP}

# rename output files;
# echo mv ${FA} ${MODEL}_${FA}
# echo mv ${FB} ${MODEL}_${FB}
# mv ${FA} ${MODEL}_${FA}
# mv ${FB} ${MODEL}_${FB}
