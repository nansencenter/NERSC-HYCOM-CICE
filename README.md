[TOC]

# Description

This repository contains the source code necessary to run the coupled HYCOM-CICE model developed at NERSC. In the repository you will find the HYCOM and CICE source codes, utility routines, and a directory containing a simple model setup (Directory TP0a1.00). Below is a short description of how to set up and compile the code

*NB: Throughout this documentation there will be references to $NHCROOT, Experiment directory, Region directory, etc
etc. A Brief explanation of these concepts can be found in [doc/Naming.md](doc/Naming.md).*

# Prerequisites

The following tools are necessary 

* A working fortran compiler.

* A working installation of the Eearth System Modelling Framework: https://www.earthsystemcog.org/projects/esmf/download/. The code has been tested and verified to work with ESMF v 5.2.0.rp3. Some info on installing ESMF can be found [here](doc/ESMF.md)

* A working python 2.7 installation and several python modules, [installation info can be found here](doc/python.md)

* Input data for the model [TODO](TODO)

# Retrieving HYCOM-CICE
Clone code from main repository 

`
git clone https://github.com/nansencenter/NERSC-HYCOM-CICE
`

If you get errors about server certificates, see [here](../..//overview#markdown-header-server-certificates)

# HYCOM-CICE directory structure
The following structure is relative to the checked out code ($NHCROOT)

    ├── bin                 # Location of binaries and python routines
    ├── cice                # Location of CICE code
    ├── doc                 # Documentation in markdown format
    ├── hycom               # Location of hycom code and utilities
    │   ├── hycom_ALL       # Location of setup/diag routines 
    │   ├── MSCPROGS        # Location of setup/diag routines  developed at NERSC
    │   └── RELO            # Location of hycom source code
    ├── input               # Location of some input files 
    └── TP0a1.00            # Location of "Reference experiment"


# Compiling HYCOM-CICE support programs

* Compile MSCPROGS routine, following the instructions given in [./hycom/MSCPROGS/](./hycom/MSCPROGS/)

* Compile HYCOM_ALL programs, following the instructions given in [./hycom/hycom_ALL/](./hycom/hycom_ALL/)

# Compiling HYCOM-CICE

Before you attempt this, make sure that the HYCOM-CICE support programs are compiled, that the MSCPROGS programs are compiled, and that you have the other prerequisites described above. The compilation is handled by the script compile_model.sh in NHCROOT/bin, and is described in [./doc/HYCOM-CICE-compilation.md](./doc/HYCOM-CICE-compilation.md)
