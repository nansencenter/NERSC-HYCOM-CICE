# Environment for compiling/running HYCOM-CICE on betzy
ml purge
ml load ESMF/8.3.0-iomkl-2022a
ml load UDUNITS/2.2.28-GCCcore-11.3.0
ml load Python/3.10.4-GCCcore-11.3.0
ml load GSL/2.7-intel-compilers-2022.1.0
ml load FFTW/3.3.10-GCC-11.3.0
ml load CMake/3.23.1-GCCcore-11.3.0

ulimit -s 2000000
