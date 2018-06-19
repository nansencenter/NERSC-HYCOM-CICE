#!/bin/bash -l

#SBATCH --account=nn9481k

#SBATCH -J  TP2a010

#SBATCH -N 32   # number of nodes
##SBATCH -n 30   # number of cores in each node 
## Each compute node has 24 cores (See more details in section Hardware on Sisu User Guide).

## System message output file
#SBATCH -o log/hycom.%J.out

## System error message file
#SBATCH -e log/hycom.%J.err

## How long job takes, wallclock time hh:mm:ss
#SBATCH -t 02:00:00

#SBATCH  --mail-type=END
#SBATCH --mail-user=jiping.xie@nersc.no

module restore system
module load NCL/6.4.0-intel-2017a
module load FFTW/3.3.6-intel-2017a

set -x

export NMPI=130
export SLURM_SUBMIT_DIR=$(pwd)
export Nbatch=8 # Note: this number is designed by the whole cores and 
                # and it is better to consistent with Ncore for each job run.
export Ncor=4 

MONITORINTERVAL=10 # time interval for periodic checks on job status


# Enter directory from where the job was submitted
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }

# ------------------- Fetch Environment ------------------------------
# -------- these are needed in preprocess scripts---------------------
echo "SLURM_JOBID    = $SLURM_JOBID     "
echo "SLURM_SUBMIT_DIR= $SLURM_SUBMIT_DIR "
[ -z "$NOMP" ] && NOMP=0

# Enter directory from where the job was submitted
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_O_WORKDIR  "; exit 1; }

iPert=1
if [ ${iPert} -eq 1 ]; then
   PerDir='/nird/home/xiejp/NERSC-HYCOM-CICE/hycom/MSCPROGS/src/Perturb_Forcing-2.2.98'
   PerDir='/cluster/work/users/xiejp/TOPAZ/TP2a0.10/expt_02.0'
   Pervars="forcing.wndewd forcing.wndnwd"
   Pervars="forcing.wndewd forcing.wndnwd forcing.airtmp forcing.mslprs" 
   #forcing.mslprs forcing.precip "
   #Pervars="forcing.airtmp forcing.mslprs forcing.precip forcing.wndewd forcing.wndnwd"
fi

if [ "$1" -gt "-1" -a "$2" -gt "-1" ]; then
  echo 'mems: ' $1 '~' $2
  mem1=$1
  mem2=$2
  (( Nens=mem2-mem1+1 ))
  echo 'mems: ' $1 '~' $2 ': ' ${Nens}
else
  exit $? 
fi

# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

# step 1:  cloning the ensemble subdirectory with perturbed forcing
for imem in `seq $mem1 $mem2`; do
   cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }
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
   [ -r summary_out ] && rm summary_out
   [ -r ./cice/ice.restart_file ] && rm ./cice/ice.restart_file
   cp ${S}/cice/* ./cice/.

   # create the perturbeted forcing
   if [ ${iPert} -eq 1 ]; then
     if [ ! -s ./force_perturb-2.2 ]; then
       ln -sf ${PerDir}/force_perturb-2.2 .
       cp ${PerDir}/infile2.in_init infile2.in
     fi
     ./force_perturb-2.2 era-i era40
     for ifile in ${Pervars}; do
       if [ -r ${ifile}.a -a -r ${ifile}.b ]; then
         rm ${ifile}.[ab]
         ln -sf tst.${ifile}.a ${ifile}.a
         ln -sf tst.${ifile}.b ${ifile}.b
       fi
     done
   fi
done
   
# step 2: submit a job
if [ ${Nens} <= ${Nbatch} ]; then
  #  a few ensemble <=Nbatch 
  for imem in `seq $mem1 $mem2`; do
     cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }
     Exdir=mem`echo 00${imem} | tail -4c`
     echo "Start time in $?: $START " $imem
     echo "End   time in $?: $END " $imem

     cd ${Exdir}/SCRATCH
     if [ ${imem} ! -eq ${mem2} ]; then
       srun --mpi=pmi2 -N $Ncor -n $NMPI --cpu_bind=cores ./hycom_cice  &
     else
       srun --mpi=pmi2 -N $Ncor -n $NMPI --cpu_bind=cores ./hycom_cice  
     fi
   done
else
  #  lot of members (>Nbatch)
  (( Ncycl = (${Nens} -1) / ${Nbatch} +1 )) 
  for icycl in `seq 1 ${Ncycl}`; do
    if [ ${icycl} -eq ${Ncycl} ]; then
       for jmem in `seq 1 ${Nbatch}`; do
         cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }
         (( ii = $mem1 + ${jmem} -1 + (${icycl} - 1) * ${Nbatch} ))
         if [ $ii -lt ${mem2} ]; then
           Exdir=mem`echo 00${ii} | tail -4c`
           cd ${Exdir}/SCRATCH
           srun --mpi=pmi2 -N $Ncor -n $NMPI --cpu_bind=cores ./hycom_cice  &
         else
           Exdir=mem`echo 00${ii} | tail -4c`
           cd ${Exdir}/SCRATCH
           sleep 30
           srun --mpi=pmi2 -N $Ncor -n $NMPI --cpu_bind=cores ./hycom_cice  
           break
         fi
       done
    else
       (( i1 = $mem1 + (${icycl} - 1) * ${Nbatch} ))
       (( i2 = $mem1 + ${Nbatch} -1 + (${icycl} - 1) * ${Nbatch} ))
       for ii in `seq $i1 $i2`; do
         cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }
         Exdir=mem`echo 00${ii} | tail -4c`
         cd ${Exdir}/SCRATCH
         echo $(pwd)
         srun --mpi=pmi2 -N $Ncor -n $NMPI --cpu_bind=cores ./hycom_cice  &
       done
       wait
   #    cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }
   #    finished=0
   #    while (( ! finished )) ; do
   #      i0=0 
   #      finished=0
   #      for ii in `seq $i1 $i2`; do
   #        cmem=mem`echo 00$ii|tail -4c`
   #        if [ -s ${cmem}/SCRATCH/summary_out ]; then
   #          cline=$(grep stop ${cmem}/SCRATCH/summary_out)
   #          if [ ${#cline} -gt 4 ]; then
   #            let i0=i0+1
   #          fi
   #        fi
   #      done
   #      if [ $i0 -eq ${Nbatch} ]; then
   #        finished=1
   #      else
   #        sleep ${MONITORINTERVAL} 
   #      fi
   #    done
    fi
  done
fi

# postprocess before to averge in the ensemble 
for imem in `seq $mem1 $mem2`; do
 cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_SUBMIT_DIR  "; exit 1; }
  Exdir=mem`echo 00${imem} | tail -4c`
  cd ${Exdir}/SCRATCH
  # Must be done under mem???/SCRATCH to run this script
  [ ! -r ../data ] && mkdir ../data
  [ ! -r ../data/cice ] && mkdir ../data/cice

  # restart. and archm
  Fhycom="restart. archm. archv."
  for ifile in ${Fhycom} ; do
  for i in $(ls -- ${ifile}*.[ab]) ; do
    if [ -L "${i}" ]; then
      rm $i
    else
      mv $i ../data/.
    fi
  done
  done

  # sea ice fields
  Fcice="iced. iceh."
  for ifile in ${Fcice}; do
  for i in $(ls -- cice/${ifile}*.nc); do
    if [ -L "${i}" ]; then
      rm $i
    else
      mv $i ../data/cice/.
    fi
  done
  done

  # Copy some files useful for analysis
  for i in $(ls regional.grid.* regional.depth.* blkdat.input ice_in cice_*.nc) ; do
     cp $i ../data/.
  done
  #
  # --- HYCOM error stop is implied by the absence of a normal stop.
  #
  if  [ `tail -1 summary_out | grep -c "^normal stop"` == 0 ] ; then
    echo "BADRUN"  > ../hycom.stop
  else
    echo "GOODRUN" > ../hycom.stop
  fi

done


exit $?
