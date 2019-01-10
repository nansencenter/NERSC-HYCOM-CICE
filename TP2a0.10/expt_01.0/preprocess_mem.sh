#!/bin/bash -l

export SLURM_SUBMIT_DIR=$(pwd)
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }

# ------------------- Fetch Environment ------------------------------
# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }


for ik in `seq 1 1` ; do
  if [ $ik -eq 0 ] ; then
    START="2001-12-25T00:00:00"
    END="2002-01-01T00:00:00"
    INITFLG="--init"
    cp -f ice_in.0 ice_in
  else
    INITFLG=""
    cp -f ice_in.1 ice_in
    if [ $ik -eq 1 ] ; then
      START="2007-07-10T00:00:00"
      END="2009-01-01T00:00:00"
    else 
      START="2010-01-02T00:00:00"
      END="2016-01-02T00:00:00"
    fi
  fi
  echo "Start time in pbsjob.sh: $START"
  echo "End   time in pbsjob.sh: $END"

  # Generate atmospheric forcing :
  atmo_synoptic.sh erai+all $START $END 

  # Transfer data files to scratch - must be in "expt_XXX" dir for this script
  expt_preprocess.sh $START $END $INITFLG        ||  { echo "Preprocess had fatal errors "; exit 1; }

done


# prepare for ensemble run ...
echo date " prepare for ensemble run ... "
Idir=$(pwd)

echo "Idir=${Idir} $SLURM_SUBMIT_DIR"

iPert=1
Pervars="forcing.airtmp forcing.mslprs forcing.precip forcing.wndewd forcing.wndnwd"

if [ ${iPert} -eq 1 -a ! -s ./force_perturb-2.2 ]; then
   PerDir='/cluster/home/xiejp/NERSC-HYCOM-CICE/hycom/MSCPROGS/src/Perturb_Forcing-2.2.98'
   ln -sf ${PerDir}/force_perturb-2.2 .
fi

if [ "$1" -gt "-1" -a "$2" -gt "-1" ]; then
  mem1=$1
  mem2=$2
  (( Nens=mem2-mem1+1 ))
  echo 'mems: ' $1 '~' $2 ': ' ${Nens}
else
  exit $? 
fi

echo "Cloning the ensemble subdirectory with perturbed forcing'"
for imem in `seq $mem1 $mem2`; do
   cd ${Idir}       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }
   # setup the model information files
   Exdir=mem`echo 00${imem} | tail -4c`
   [ ! -r ${Exdir} ] && mkdir ${Exdir}
   [ -r ${Exdir}/hycom.stop ] && rm -rf ${Exdir}/hycom.stop
   [ -r ${Exdir}/SCRATCH ] && rm -rf ${Exdir}/SCRATCH
   [ ! -r ${Exdir}/SCRATCH ] && mkdir ${Exdir}/SCRATCH

   # link the neccesary model file
   cd ${Exdir}/SCRATCH
   [ ! -r ./cice ] && mkdir cice 
   ln -sf ${S}/* .
   ln -sf ${S}/cice/* ./cice/.

   [ -r summary_out ] && rm summary_out
   [ -r ./cice/ice.restart_file ] && rm ./cice/ice.restart_file
   rm restart.*.?
   for f in ${S}/restart.*_${Exdir}.? ; do
     Nf=${#f}
     ((Nf0 = Nf - 33 ))
#     echo $f $Nf ${f:Nf0:24}${f:Nf-2:2} 
     ln -sf $f ${f:Nf0:24}${f:Nf-2:2} 
   done
   rm ./cice/iced.*.nc
   for f in ${S}/cice/iced.*_${Exdir}.nc ; do
      Nf=${#f}
      ((Nf0 = Nf -31 ))
    #  echo $f $Nf $Nf0 ${f:Nf0:21}.nc
      ln -sf $f ./cice/${f:Nf0:21}.nc
      echo "./cice/${f:Nf0:21}.nc" > ./cice/ice.restart_file
   done

  # create the perturbeted forcing
   if [ ${iPert} -eq 1 ]; then
     if [ ! -s ./force_perturb-2.2 ]; then
       ln -sf ${Idir}/force_perturb-2.2 .
       cp ${Idir}/infile2.in_init infile2.in
     fi
     ./force_perturb-2.2 era-i era40
     echo " imem=${imem}"
     echo "replace ${Pervars} ... "
     for ifile in ${Pervars}; do
       if [ -r tst.${ifile}.a -a -r tst.${ifile}.b ]; then
         rm ${ifile}.[ab]
         ln -sf tst.${ifile}.a ${ifile}.a
         ln -sf tst.${ifile}.b ${ifile}.b
       fi
     done
   fi
done
   








exit $?

