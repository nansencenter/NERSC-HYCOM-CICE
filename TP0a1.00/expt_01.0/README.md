[toc]

This is the experiment directory. It will contain necessary input files as well 
as the hycom source and executable, and the output files. Below is a description of the
subdirectories

    ├── build     # where the model is compiled
    ├── data      # where model output is stored
    ├── log       # logs ...
    └── SCRATCH   # where the model runs

# HYCOM and CICE input files

The HYCOM input file is 

    blkdat.input
    
    Specifications about heat flux options in blkdat:
    c --- 'flxflg' = thermal forcing flag (0=none,3=net_flux,1-2,4-6=sst-based)
    c --- (=1 MICOM bulk parameterization)
    c --- (=2 Kara bulk parameterization)
    c --- (=4 COARE bulk parameterization, approx.)
    c --- (=5 L&Y bulk parameterization)
    c --- (=6 COARE bulk parameterization, approx., better pressure)

    flxflg=3 is using the net flux (say, from ERA-Interim) but traditionally we rather 
    use the formulas that use the model variables (SST, ice etc). 

    flxflg=5 is the real CORE2 described by Large and Yeager 2009 (L&Y): 
    https://rda.ucar.edu/datasets/ds260.2/docs/OSGC-000-000-003-157.pdf

    flxflg 4 and 6 are two versions of COARE (Fairall et al.  2003) 
    
    flxflg=6 takes sea level pressure as input to the drag coefficient and uses wind-ocean 
    speed  http://journals.ametsoc.org/doi/abs/10.1175/1520-0442(2003)016%3C0571%3ABPOASF%3E2.0.CO%3B2 

    empflg = 4,5,6 are all sst-based because they use the model SST to calculate evaporation. 
    Evaporation is not big in the Arctic, but still we should use consistent flags for empflg 
    and flxflg (and wndflg as well!) for the sake of consistency (except that flxflg=6 implies 
    wndflg=4, its the same).

The CICE input file is
 
    ice_in


# Environment setup (EXPT.src)

First of all, for convenience,  source the region file in the parent directory ($REGION.src). This will set up a few
variables, but most importantly, it will add NHCROOT/bin to your path. That way you dont have to
type the full path to commands.

Scratch/data directories are set in EXPT.src. Scratch directories are where the
models run, the expt_preprocess.sh script (see below) will copy necessary files to that location
before running the model. Data directories are where model output is stored. The
expt_postprocess.sh (see below) script will move data files from the scratch directory to the
data directory. It is possible to modify these by setting these environment variables
in EXPT.src:

    export D=$P/data     # Where data ends up
    export S=$P/SCRATCH  # Scratch area 

by default these are set relative to the experiment directory, given by environment variable $P.

Note that by sourcing the REGION.src file (in the region directory, above this one)
and the EXPT.src file (in the experiment directory), the shell wil set the
variables S and D which point to the scratch and data directories respectively.
type 

    cd $S  

or

    cd $D

to go to the respective directories.

Most of the shell scripts in $NHCROOT/bin/ source the REGION.src and EXPT.src files before 
proceeding.

# HYCOM environment settings in EXPT.src

Some HYCOM settings are controlled by EXPT.src. Below is a list

    export SIGVER=2      # Version of equation of state (this is 7-term sigma 2). Must not cause conflict with thflag in blkdat.input
    export NMPI=4        # Number of MPI threads to use
    X="01.0"             # Experiment number, it should match experiment directory name (expt_01.0)
    E="010"              # E is X without "."
    T="01"               # Topography version
    export V=2.2.98      # hycom version              

The SIGVER variable is used by the compile script to get the correct equation of state (it basically modifies stmt_fns.h in the
hycom build directory). The $T variable is used to choose the bathymetries to use (you can have different bathymetries for different experiments).

# Compiling HYCOM-CICE

HYCOM-CICE compilation is complicated, since we have to compile two models, find the correct configuration files, etc etc. A script exists that is meant to take care of all this, and can be found in $NHCROOT/bin/compile_model.sh. It can be run without arguments, to get a brief explanation of its arguments. Here is some typical use cases

    #compile for hexagon using portland compilers (NB: module PrgEnv-pgi must be loaded!)
    $NHCROOT/bin/compile_model.sh pgi
    
    #compile for generic linux with gnu fortran and openmpi mpi  libraries.
    $NHCROOT/bin/compile_model.sh -m openmpi pgi

To configure for other architectures you will need to create macro files suitable for your system. For a more
in-depth description of the compile script, see  [TODO](TODO)

The compiled executable will be placed inside the build/src... directory as hycom_cice
    

# Preprocessing (script expt_preprocess.sh in $NHCROOT/bin)

This script will fetch the files needed to run the model and place them in the 
SCRATCH directory. * This script can be run interactively, this is useful as you can catch
most errors before starting a job *

In most cases you will be alerted if there is a missing files
or if there are inconcistencies, but it will not catch every error. The expt_preprocess.sh is called like this from the
experiment directory:

   expt_preprocess.sh 2015-01-01T00:00:00 2015-01-10T00:00:00 

or 

   expt_preprocess.sh 2015-01-01T00:00:00 2015-01-10T00:00:00  --init

The optional init flag tells that the script is to set up hycom from a initial state using climatology.
*NB: In new hycom versions, model initialization works even if hycom yrflag==3.*

You will have to set the iniflg in blkdat.input appropriately if you wish to initialize
the hycom model. Depending on the value of iniflg, you may also have to create relaxation 
fields first.


# Running jobs

Before jobs are run the script expt_preprocess.sh must be run, this can be set in
the job script. You can also run this script interactively from the experiment
directory, to make sure you have all data files you need. It should cover most
of the data files we use, and handle most (but not all) inconsistencies. If there
is a inconsistency not caught by the script, it will eventually be caught by hycom.

If you add new data files for the model, you will have to modify expt_preprocess.sh
so that the script copies them to the scratch directory. Otherwise you must copy
them by hand to the scratch directory.

expt_postprocess.sh copies files from the scratch directory to your data directory.
If you find some files are not copied from the  scratch directory you will have
to modify expt_postprocess.sh. Also, if a job runs out of time, expt_postprocess.sh will
not be run. In that case you can run the expt_postprocess.sh script interactively to
retrieve the files to the data directory.

Job scripts will have to be modified from machine to machine, but the
expt_preprocess.sh and expt_postprocess.sh scripts should be able to run on all machines.

# create_ref_case.sh

This script sets up basic input files .... TODO

# Submitting jobs

A example job script for the pbs system is given in pbsjob.sh. It is set up for hexagon,
but with some modifications it should work on most cray XT systems.

# Example EXPT.src file 


    # Example EXPT.src file
    #!/bin/bash
    #
    # Sets environment for this experiment, also makes the scratch (S) and
    # data (D) directories if not already present

    #
    #
    # --- R is region name.
    # --- V is source code version number.
    # --- T is topography number.
    # --- K is number of layers.
    # --- E is expt number.
    # --- P is primary path.
    # --- D is permanent directory.
    # --- S is scratch   directory, must not be the permanent directory.
    #
    # hycom executable will be retrieved from Build_V${V}_X${X}. ex: Build_V2.2.12_X01.0
    #
    mydir=$(cd $(dirname ${BASH_SOURCE}) && pwd)
    unset -v X E T V K P D S
    X="01.0"                # X based on dir name (expt_03.1)
    E="010"                 # E is X without "."
    T="01"                                                           # Topography version
    export V=2.2.98                                                  # hycom version              
    #export K=`grep "'kdm   ' =" blk* | awk '{printf("%03d", $1)}'`   # get kdm from blkdat
    export K=`grep "'kdm   ' =" $mydir/blkdat.input | awk '{printf("%03d", $1)}'`   # get kdm from blkdat
    export P=$mydir                                                  #  ---""---
    export D=$P/data                                                 # Where data ends up
    export S=$P/SCRATCH                  # Scratch area 

    export SIGVER=2   # Version of equation of state (this is 7-term sigma 2). Must not cause conflict with thflag in blkdat.input
    export NMPI=4

    # Consistency check. Ensures expt dir ends in expt_X
    #echo $tmp
    tmp=$(basename $P)
    if [ "$tmp" != "expt_${X}" ] ;then
       echo "Error: Mismatch between path of experiment $P and assumed name expt_${X}"
       exit 1
    fi

