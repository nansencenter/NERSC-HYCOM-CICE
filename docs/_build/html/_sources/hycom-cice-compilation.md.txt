

The compilation of the nersc-hycom-cice model is handled by the compile script,
simply called compile_model.sh, located in $NHCROOT/bin. This script does a lot
of things "behind the carpet", so here is a quick description of what it
actually does.


The compile script has to be run in the experiment
dir, and a local directory "build" will be created that contains a copy of the
source code (copied from $NHCROOT/hycom/RELO/).  Its inside this directory that
the compilation of the code takes place and final executable will be located.


#  Invocation of compile_model.sh

Running the compile_model.sh script without arguments will provide you with a
brief description of how to run it. Here is the short version:

      Example:
         $(basename $0) [ -u ]  [ -m mpi_library ]  compiler

      arguments   :
         compiler : compiler to use for. Currently supported: ifort gfortran pgi

      optional arguments :
         -u              : update code in build dir from $sourcedir
         -m mpi_library  : on some machines you need to specify what mpi library to use

The only necessary argument is the compiler to use. Accepted values are pgi,
gfortran, intel. The compile script needs this to choose the correct macro
files when compiling the model.

The optional argument "-u" will make sure that you copy the most recent
code from the $NHCROOT/hycom/RELO/ directory into the local build directory. If
you run without the "-u" option, you will make no changes to the source code
inside the build directory. (This is useful if you make local changes and want
to keep them before possibly copying these changes to the "main" directory at
$NHCROOT/hycom/RELO/).

The last optional argument "-m" determines what openmpi framework to use. On
hexagon (and meny other big HPC installations), the openmpi framework is set by loading modules, so this argument is
not needed in this case. On a mopre "generic" system, such as your local linux
installation, you will have to specify what mpi framework to use (typical values
are mpich, mpich2, openmpi, lam, intempi, etc etc)....  The compile script needs
this  information mainly to choose the correct location of  ESMF fortran modules when compiling the model.

# compile_models.sh, step by step 

## Initial stuff...
Initially, the compile scripts set somelocal variables, such as ARCH
(architecture) and SITE /(sitename) using the linux uname command.

Using these variables, a MACROID variable is set. This variable is used for
getting the correct Makefile macros when compiling the HYCOM and COCE model
code.

## ESMF setup

Next, the compile script tries to get the location of the ESMF libraries and the
ESMF include files. Three variables will be set (these are used in the makefile
macros). The variables are 

| Variable | Description |
|----------|-------------|
| ESMF_DIR | Top-level directory of the ESMF installation |
| ESMF_LIB_DIR | Directory of the ESMF library |
| ESMF_MOD_DIR | Directory of the ESMF fortran modules |

Three cases arise:

* If you have set ESMF_DIR, ESMF_MOD_DIR and ESMF_LIB_DIR before callign the
  compile script, then these will be used as they are.

* Hardcoded cases (only  hexagon for now) set ESMF_MOD_DIR and ESMF_LIB_DIR
  based on ESMF_DIR

* Generic cases set ESMF_MOD_DIR and ESMF_LIB_DIR based on ESMF_DIR and the mpi
  library you use.


## Source from EXPT.src

From the experiment source file (EXPT.src), the script will get equation of
state to use, based on the variable SIGVER. It also sets the number of terms
used by the equation of state.

## copy code from $NHCROOT/hycom/RELO/

If the build dir does not exist, then an rsync command is run to get the HYCOM
and CICE code from  $NHCROOT/hycom/RELO/. If the directory does exist, no code
will be copied unless tyhe "-u" option is provided to the script.

## hycom_feature_flag

It is possible to specify hycom preprocessor flags inside this file to activate or deactivate
features in the hycom code. If the file exists, it will be read and passed on to
the C preprocessor command.

# Set up eq of state

The equation of state will be set up next, and linked to the correct statement
function  based on the choice of SIGVER. This include file is needed when
compiling hycom.

# Compile CICE

The compilation of CICE then follows, using variables set up earlier in the
script. Note that the grid size needs to be passed to CICE at this stage.

# Compile HYCOM

The final step is to compile the HYCOM source code. This also compiles the main
hycom-cice executable, and links the CICE and HYCOM object files to create the
hycom_cice executable.
