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
process_daily_files() {

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
       ln -s $folder"/hycave" $trunk"/hycave" || tellerror "Could not get hycave"
       ln -s $folder"/hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
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
    ln -s $folder"/proj.in" $trunk"/proj.in" 
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
    mv $trunk"/archm_"$year"0101_12.nc" $folder'/AVE.'$month'.'$year'.'$year'.nc'


    mv $trunk"/archm."$year"_001_12.a" $folder'/AVE.'$month'.'$year'.'$year'.a'
    mv $trunk"/archm."$year"_001_12.b" $folder'/AVE.'$month'.'$year'.'$year'.b'

    cd $basedir

}

# Function to process a subset of files
process_monthly_files() {

    local region=${1}
    local experiment=${2}
    local year1=${3}
    local year2=${4}
    local month=${5}
    local basedir=${6}
    local workdir=${7}

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
    trunk=$folder"AVE"$month

    echo $trunk $month

    mkdir -p $trunk

    if [ "$region" = "TP5" ]  || [ "$region" = "TP2" ] || [ "$region" = "TP0" ] || [ "$region" = "NAT" ] ; then


       rm -rf $trunk"/hycave" $trunk"/hyc2proj" 
       ln -s $folder"/hycave" $trunk"/hycave" || tellerror "Could not get hycave"
       ln -s $folder"/hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
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
    ln -s $folder"/proj.in" $trunk"/proj.in" 
    ln -s $folder"/extract.archm" $trunk"/extract.archm"
    ln -s $folder"/depthlevels.in" $trunk"/depthlevels.in" 

    y=$year1

    append=" "

    while (( y <= year2 )); do
        cp -fu $folder"/AVE."$month"."$y"."$y".a" $trunk"/archm."$y"_001_12.a"
        cp -fu $folder"/AVE."$month"."$y"."$y".b" $trunk"/archm."$y"_001_12.b" 
        append=$append"archm."$y"_001_12.a "
        (( y++ ))
    done
    echo $execute$append >> $trunk"/archfiles.txt"

    cd $trunk"/"
    $execute$append

    mv $trunk'/AVE.a' $trunk"/archm.2200_001_12.a" # necessary to have the name in archm format
    mv $trunk'/AVE.b' $trunk"/archm.2200_001_12.b"

    ./hyc2proj "archm.2200_001_12.b"
    mv $trunk"/archm_22000101_12.nc" $folder'/CLIM.'$month'.'$year1'.'$year2'.nc'    


    mv $trunk"/archm.2200_001_12.a" $folder'/CLIM.'$month'.'$year1'.'$year2'.a'
    mv $trunk"/archm.2200_001_12.b" $folder'/CLIM.'$month'.'$year1'.'$year2'.b'

    rm $trunk"/archm.*"

    cd $basedir

}


# Function to process a subset of files
process_yearly_files() {

    local region=${1}
    local experiment=${2}
    local year=${3}
    local basedir=${4}
    local workdir=${5}

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
    trunk=$folder"AVE"$year

    echo $trunk $year

    mkdir -p $trunk

    if [ "$region" = "TP5" ]  || [ "$region" = "TP2" ] || [ "$region" = "TP0" ] || [ "$region" = "NAT" ] ; then

       rm -rf $trunk"/hycave" $trunk"/hyc2proj" 
       ln -s $folder"/hycave" $trunk"/hycave" || tellerror "Could not get hycave"
       ln -s $folder"/hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
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
    ln -s $folder"/proj.in" $trunk"/proj.in" 
    ln -s $folder"/extract.archm" $trunk"/extract.archm"
    ln -s $folder"/depthlevels.in" $trunk"/depthlevels.in" 

    append=" "

    for month in $(seq -w 1 12); do
        d=$year"-"$month"-15"
        n=$(date -d "$d" '+%j')
        printf -v nn "%03d" "$((10#$n))"

        cp -fu $folder"/AVE."$month"."$year"."$year".a" $trunk"/archm."$year"_"$nn"_12.a"
        cp -fu $folder"/AVE."$month"."$year"."$year".b" $trunk"/archm."$year"_"$nn"_12.b" 
        append=$append"archm."$year"_"$nn"_12.a "
    done

    echo $execute$append >> $trunk"/archfiles.txt"

    cd $trunk"/"
    $execute$append

    mv $trunk'/AVE.a' $trunk"/archm.2200_001_12.a" # necessary to have the name in archm format
    mv $trunk'/AVE.b' $trunk"/archm.2200_001_12.b"

    ./hyc2proj "archm.2200_001_12.b"
    mv $trunk"/archm_22000101_12.nc" $folder'/AVE.'$year'.'$year'.nc'    


    mv $trunk"/archm.2200_001_12.a" $folder'/AVE.'$year'.'$year'.a'
    mv $trunk"/archm.2200_001_12.b" $folder'/AVE.'$year'.'$year'.b'

    rm $trunk"/archm.*"

    cd $basedir

}

# Function to process a subset of files
process_model_climatology() {

    local region=${1}
    local experiment=${2}
    local year1=${3}
    local year2=${4}
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
    trunk=$folder"AVEALL"

    echo $trunk $year1 $year2

    mkdir -p $trunk

    if [ "$region" = "TP5" ]  || [ "$region" = "TP2" ] || [ "$region" = "TP0" ] || [ "$region" = "NAT" ] ; then

       rm -rf $trunk"/hycave" $trunk"/hyc2proj" 
       ln -s $folder"/hycave" $trunk"/hycave" || tellerror "Could not get hycave"
       ln -s $folder"/hyc2proj" $trunk"/hyc2proj" || tellerror "Could not get hycave"
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
    ln -s $folder"/proj.in" $trunk"/proj.in" 
    ln -s $folder"/extract.archm" $trunk"/extract.archm"
    ln -s $folder"/depthlevels.in" $trunk"/depthlevels.in" 

    append=" "

    for year in $(seq $year1 $year2); do

        cp -fu $folder"/AVE."$year"."$year".a" $trunk"/archm."$year"_001_12.a"
        cp -fu $folder"/AVE."$year"."$year".b" $trunk"/archm."$year"_001_12.b" 
        append=$append"archm."$year"_001_12.a "        

    done

    echo $execute$append >> $trunk"/archfiles.txt"

    cd $trunk"/"
    $execute$append

    mv $trunk'/AVE.a' $trunk"/archm.2200_001_12.a" # necessary to have the name in archm format
    mv $trunk'/AVE.b' $trunk"/archm.2200_001_12.b"

    ./hyc2proj "archm.2200_001_12.b"
    mv $trunk"/archm_22000101_12.nc" $folder'/CLIM.'$year1'.'$year2'.nc'    


    mv $trunk"/archm.2200_001_12.a" $folder'/CLIM.'$year1'.'$year2'.a'
    mv $trunk"/archm.2200_001_12.b" $folder'/CLIM.'$year1'.'$year2'.b'

    rm $trunk"/archm.*"

    cd $basedir


}


# Export the function so it can be used by parallel
export -f process_daily_files
export -f process_monthly_files
export -f process_yearly_files
export -f process_model_climatology

for year in $(seq $year1 $year2); do
    for month in $(seq -w 1 12); do

        srun --exclusive -n1 -c1 bash -c "process_daily_files $region $experiment $year $month $basedir $workdir" &

    done
done

wait


for month in $(seq -w 1 12); do

        srun --exclusive -n1 -c1 bash -c "process_monthly_files $region $experiment $year1 $year2 $month $basedir $workdir" &

done

wait


for year in $(seq $year1 $year2); do

        srun --exclusive -n1 -c1 bash -c "process_yearly_files $region $experiment $year $basedir $workdir" &

done

wait

srun --exclusive -n1 -c1 bash -c "process_model_climatology $region $experiment $year1 $year2 $basedir $workdir" &

wait


# cd $folder
# for month in $(seq -w 1 12); do
#    d="2200-"$month"-15"
#    n=$(date -d "$d" '+%j')
#    printf -v nn "%03d" "$((10#$n))"
#    cp $folder'/CLIM.'$month'.'$year1'.'$year2'.a'  $folder"/archm.2200_"$nn"_12.a"
#    #cp $folder'/CLIM.'$month'.'$year'.'$year'.a' $folder
# done

# srun --exclusive -n1 -c1 "$folder/hycave" archv $folder"/archm.2200_*_12.a" 
# mv $folder"AVE.a" $folder'/CLIM.00.'$year1'.'$year2'.a' 
# mv $folder"AVE.b" $folder'/CLIM.00.'$year1'.'$year2'.b'
# rm $folder"/archm.2200_*_12.a"
# wait

