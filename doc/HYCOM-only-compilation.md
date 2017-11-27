

The compilation of a standalone HYCOM from the NERSC-HYCOM-CICE model is mainly carried out using a modified compile script in "$NHCROOT/bin", so-called compile_model.sh. This script has been modified to account for standalone compilation of HYCOM. Here, I will explain briefly all required procedures and processes have been applied to make this compilation up and running. Note the general description can be found in "$NHCROOT/HYCOM-CICE-compilation.md".

The compile script has to be run in the experiment
dir, and a local directory "build" will be created that contains a copy of the
source code (copied from $NHCROOT/hycom/RELO/).  Its inside this directory that
the compilation of the code takes place and final executable will be located.


#  Invocation of compile_model.sh

      Example:
         $(basename $0) [ -u ]  [ -m mpi_library ]  compiler

      arguments   :
         compiler : compiler to use for. Currently supported: ifort gfortran pgi

      optional arguments :
         -u              : update code in build dir from $sourcedir
         -m mpi_library  : on some machines you need to specify what mpi library to use


# compile_models.sh, step by step 

## Initial stuff...
An environment variable called MACROID needs to be specified according to the following naming role:

MACROID=$ARCH.$SITE.$compiler

where ARCH denotes architecture (e.g. Linux, Darwin,…), SITE is deduced from hostname using the line senate command (e.g. Hexagon, Sisu, Fram,…), and compiler denotes the name of compiler (e.g. fortran, intel, pig,….). You may need to modify to include your specific site within the "compile_model.sh"

## Changes in $NHCROOT/hycom/RELO/config

As next step, we include a configuration file called $MACROID_hycom (e.g. Linux.sisu.intel_hycom) in the location of $NHCROOT/hycom/RELO/config. This file is quite similar to the $MACROID_cice (where export TYPE=cice), but we have removed the flags related to the ESMF which is introduced as (for example for Sisu HPC):

INCLUDES      = -I$(MPI_HOME)/include  -I${ESMF_MOD_DIR} -I${ESMF_DIR}/include -I./CICE/rundir/compile


CPPFLAGS      = -DIA32 -DREAL8 -DMPI -DSERIAL_IO -DNAN2003 -DTIMER -DRELO  -DUSE_ESMF -DUSE_ESMF_5 -DNERSC_HYCOM_CICE
EXTRALIBS     = -L${ESMF_LIB_DIR}/ -lesmf -lnetcdf
  

Above command lines are simply modified as

INCLUDES      = -I$(MPI_HOME)/include

CPPFLAGS      = -DIA32 -DREAL8 -DMPI -DSERIAL_IO -DNAN2003 -DTIMER -DRELO -DNERSC_HYCOM_CICE



* You need to take of changes according to your HPC and compiler.

## Changes in $NHCROOT/hycom/RELO/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi

In this folder, we replace the old "Makefile" by the new one which modified to include proper flags for compilation of standalone HYCOM.

Another script to compile model, so-called Make-hycom.csh, is copied in this folder. This script is called by "compile_model.sh" when "iceflg" is set to zero in "blkdat.input".



## Changes in experiment folder

Since the model has been meant for a coupled system, there are few inconsistency in two source codes in $NHCROOT/hycom/RELO/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi (i.e. hycom.F and mod_hycom.F). To avoid any mofdification of source codes at this stage, I created a folder called "mysource" including modified hycom.F and mod_hycom.F which are copied to the "$EXPT_PATH/build/src_2.2.98ZA-07Tsig0-i-sm-sse_relo_mpi" during compile time only when "iceflg" is zero in "blkdat.input". Furthermore, you need to have ice_in in experiment folder (this is get read in initialisation which later on be untangled from the standalone compilation).

# Compile HYCOM

The final step is to compile the HYCOM source code depending on the status of "iceflg" flag in "blkdat.input". If "iceflg" is non-zero, the model compiles and links both HYCOM and CICE models to generate of hycom-cice executable. If "iceflg" is zero the hycom-cice executable is generated from only compiling HYCOM model.

