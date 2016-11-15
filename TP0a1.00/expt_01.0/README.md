This is the experiment directory. It will contain the necessary input files as well 
as the hycom source and executable, and the output files. There are directories for logs, 
scratch (where the model runs) and data (where output data is placed).

# Environment

Scratch/data directories are set in EXPT.src. Scratch directories are where the
models run, the expt_preprocess.sh script will copy necessary files to that location
before running the model. Data directories are where model output is stored. The
expt_postprocess.sh script will move data files from the scratch directory to the
data directory. It is possible to modify these by setting these environment variables
in EXPT.src:

    export D=$P/data     # Where data ends up
    export S=$P/SCRATCH  # Scratch area 

by default these are set relative to the experiment directory, given by environment variable $P.

Note that by sourcing the REGION.src file (in the top-level region directory)
and the EXPT.src file (in the experiment directory), the shell wil set the
variables S and D which point to the scratch and data directories respectively.
type 

    cd $S  

or

    cd $D

to go to the respective directories.

Most of the scripts in ../../bin/ source the REGION.src and EXPT.src files before 
proceeding.

# HYCOM environment settings

Some HYCOM environment settings are controlled by EXPT.src. This is the setting for equation of state, and the number of MPI
tasks to use. These are set by setting these variables in EXPT.src:

    export SIGVER=2   # Version of equation of state (this is 7-term sigma 2). Must not cause conflict with thflag in blkdat.input
    export NMPI=4     # Number of MPI threads to use

The SIGVER variable is used by the compile script to get the correct equation of state (it basically modifies stmt_fns.h in the
hycom build directory).

# Compiling HYCOM-CICE

HYCOM-CICE compilation is complicated, since we have to compile two models, find the correct configuration files, etc etc. A script exists that is meant to take care of all this, it is called "compile_model.sh


# expt_preprocess.sh

This script will fetch the files needed to run the model and place them in the 
SCRATCH directory. * This script can be run interactively, this is useful as you can catch
most errors before starting a job *

In most cases you will be alerted if there is a missing files
or if there are inconcistencies, but it will not catch every error. The expt_preprocess.sh is called like this:

   ../../bin/expt_preprocess.sh 2015-01-01T00:00:00 2015-01-10T00:00:00 

or 

   ../../bin/expt_preprocess.sh 2015-01-01T00:00:00 2015-01-10T00:00:00  --init

The optional init flag tells that the script is to set up hycom from a initial state.
*NB: In new hycom versions, model initialization works even if hycom yrflag==3.*

You will have to set the iniflg in blkdat.input appropriately if you wish to initialize
the model. Depending on the value of iniflg, you may also have to create relaxation 
fields first.


# Running jobs

Before jobs are run the scripts expt_preprocess.sh must be run, this can be set in
the job queue script. You can also run this script from the experiment
directory, to make sure you have all data files you need. It should cover most
of the data files we use, and handle most (but not all) inconsistencies. If there
is a inconsistency not caught by the script, it will eventually be caught by hycom.

If you add new data files for the model, you will have to modify preprocess.sh
so that the script copies them to the scratch directory. Otherwise you must copy
them by hand to the scratch directory.

postprocess.sh copies files from the scratch directory to your data directory.
If you find some files are not copied from the  scratch directory you will have
to modify postprocess.sh. Also, if a job runs out of time, expt_postprocess.sh will
not be run. In that case you can run the postprocess.sh script interactively to
retrieve the files to the data directory.

Job scripts will have to be modified from machine to machine, but the
expt_preprocess.sh /expt_postprocess.sh scripts should be able to run on all machines.

# create_ref_case.sh

This script sets up basic input files .... TODO

# Submitting jobs

A example job script for the pbs system is given in pbsjob.sh. It is set up for hexagon,
but with some modifications it should work on most cray XT systems.

