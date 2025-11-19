#!/bin/bash

#SBATCH -J HycAVE
#SBATCH --output=HYCAVE_%J.out
#SBATCH --error=HYCAVE_%J.err
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=24:00:00  # Adjust the time limit as needed
#SBATCH --account=nn9481k
#SBATCH --qos=preproc
#SBATCH --mem-per-cpu=7GB
ml purge
ml load CMake/3.23.1-GCCcore-11.3.0
ml load ESMF/8.3.0-iomkl-2022a
ml load FFTW/3.3.10-GCC-11.3.0
ml load UDUNITS/2.2.28-GCCcore-11.3.0
ml load Python/3.10.4-GCCcore-11.3.0
ml load GSL/2.7-intel-compilers-2022.1.0

export basedir=$(pwd)

region=$1
experiment=$2
year1=$3
year2=$4
workdir=$5

# Function to process a subset of files
process_files() {

    #workdir=$WORK"/" # sigma2 computers assign this as /clustee/work/users/$USER

    local region=${1}
    local experiment=${2}
    local year=${3}
    local month=${4}
    local basedir=${5}
    local workdir=${6}

    if [ "$region" = "TP5" ] ; then
       reg="TP5a0.06"
    elif [ "$region" = "TP2" ] ; then
       reg="TP2a0.10"
    elif [ "$region" = "TP0" ] ; then
       reg="TP0a1.00"
    elif [ "$region" = "NAT" ] ; then
       reg="NATa1.00"
    fi

    expt="expt_"$(echo $experiment | cut -c1-2)"."$(echo $experiment | cut -c3)
    folder=$workdir$reg"/"$expt"/data/"
    trunk=$folder$year$month

    echo $trunk $year $month

    mkdir -p $trunk
    if [ "$region" = "TP5" ]  || [ "$region" = "TP2" ] || [ "$region" = "TP0" ] || [ "$region" = "NAT" ] ; then


       rm -rf $trunk"/hycave" $trunk"/hyc2proj" 
#       MSPROGSbinPATH="/cluster/home/cagyum/NERSC-HYCOM-CICE/hycom/MSCPROGS/bin/"
#       cp $MSPROGSbinPATH"hycave" $trunk"/hycave" || tellerror "Could not get hycave"
#       cp $MSPROGSbinPATH"hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
#       ln -s $MSPROGSbinPATH"/hycave" $trunk"/hycave" || tellerror "Could not get hycave"
#       ln -s $MSPROGSbinPATH"/hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
       ln -s $folder"/hycave" $trunk"/hycave" || tellerror "Could not get hycave"
       ln -s $folder"/hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
#       ln -s "/cluster/work/users/cagyum/TP2a0.10/expt_03.8/data/hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
       execute=" ./hycave archv "
    fi

    rm -f $trunk"/regional.grid.a" $trunk"/regional.grid.b" 
    rm -f $trunk"/regional.depth.a" $trunk"/regional.depth.b" 
    rm -f $trunk"/proj.in" $trunk"/extract.archm" $trunk"/depthlevels.in"  
    rm -f $trunk"/archfiles.txt"
    rm -f $trunk"/grid.info"

    ln -s $folder"/regional.grid.a" $trunk"/regional.grid.a"
    ln -s $folder"/regional.grid.b" $trunk"/regional.grid.b"
    ln -s $folder"/regional.depth.a" $trunk"/regional.depth.a"
    ln -s $folder"/regional.depth.b" $trunk"/regional.depth.b"
    ln -s $folder"/grid.info" $trunk"/grid.info"
    ln -s $folder"/proj.in" $trunk"/proj.in" # projection fits WOA grid
    ln -s $folder"/extract.archm" $trunk"/extract.archm"
    ln -s $folder"/depthlevels.in" $trunk"/depthlevels.in"    

    d=$year"-01-01"
    next_year=$((year + 1))

    append=" "
    while [ "$d" != $next_year"-01-01" ]; do

        n=$(date -d "$d" '+%j')
        m=$(date -d "$d" '+%m')
        y=$(date -d "$d" '+%Y')

        if [ "$region" = "TP5" ] || [ "$region" = "TP2" ] || [ "$region" = "TP0" ]; then
            if [ "$m" == "$month" ]; then
                printf -v nn "%03d" "$((10#$n))"
                append=$append"../archm."$y"_"$nn"_12.a "
            fi
        fi
        d=$(date -I -d "$d + 1 day")
    done
    echo $execute$append >> $trunk"/archfiles.txt"

    cd $trunk"/"
    $execute$append

    

    mv $trunk'/AVE.a' $trunk"/archm."$year"_001_12.a" # necessary to have the name in archm format
    mv $trunk'/AVE.b' $trunk"/archm."$year"_001_12.b"

    ./hyc2proj "archm."$year"_001_12.b"
    mv $trunk"/archm_"$year"0101_12.nc" $folder'/AVE.WOA.'$month'.'$year'.'$year'.nc'


    mv $trunk"/archm."$year"_001_12.a" $folder'/AVE.'$month'.'$year'.'$year'.a'
    mv $trunk"/archm."$year"_001_12.b" $folder'/AVE.'$month'.'$year'.'$year'.b'

    cd $basedir

}


# Export the function so it can be used by parallel
export -f process_files

for year in $(seq $year1 $year2); do
    for month in $(seq -w 1 12); do
#        for experiment2 in 030 033 034; do

	#echo $year $month
#        srun -n1 -c1 bash -c "process_files $region $experiment2 $year $month $basedir" &
        srun -n1 -c1 bash -c "process_files $region $experiment $year $month $basedir $workdir" &
#        bash -c "process_files $region $experiment $year $month $basedir" &
        done
    done
#done

wait
