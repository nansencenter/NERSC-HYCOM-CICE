# setup.sh: get regional.grid.[ab] files and find pivot points
wmsc_dir="/work/shared/nersc/msc"
M_INP="/work/shared/nersc/msc/ModelInput"

# data source name
# (should match $FA, $FB and $PP file names below)
dname=wamnsea
#where to put the pivot files
outdir=.
# outdir=$wmsc_dir/WAVES_INPUT/WAMNSEA/pivots

# name of script to find pivot points
PP="../"${dname}"-piv-hyc2.2.12"

# Model-dependent stuff
# (Do all models at once, or adjust the sequence below (`seq ...`))
MOD[1]="TP4a0.12"
tdir[1]=${M_INP}/TOPAZ4/${MOD[1]}/topo
MOD[2]="BS1a0.045"
tdir[2]=${M_INP}/Barents_Hyc2.2.12/topo/
MOD[3]="FR1a0.03"
tdir[3]=${M_INP}/FramStrait_Hyc2.2.12/${MOD[3]}/topo

for num in `seq 3 3`
do
   MODEL=${MOD[$num]}
   topo_dir=${tdir[$num]}
   echo " "
   echo "*****************************************************************************"
   echo $MODEL
   echo $topo_dir
   echo "*****************************************************************************"
   echo " "

   ln -sf ${topo_dir}/regional.grid.a .
   ln -sf ${topo_dir}/regional.grid.b .
   ln -sf ${topo_dir}/regional.depth.a .
   ln -sf ${topo_dir}/regional.depth.b .
   
   # run $PP to get pivot points;
   echo " "
   echo ${PP}
   echo " "
   ${PP}
   
   # name of outputs from script:
   # (NB: these should be the same as 'filenamea' and 'filenameb' in ../p_getpivots.F90)
   FA="pivots_"${dname}".a"
   FB="pivots_"${dname}".b"
   
   # rename output files;
   echo " "
   echo mv ${FA} ${outdir}/${MODEL}_${FA}
   echo mv ${FB} ${outdir}/${MODEL}_${FB}
   echo " "
   mv ${FA} ${outdir}/${MODEL}_${FA}
   mv ${FB} ${outdir}/${MODEL}_${FB}
done

#clean up:
rm regional.*
