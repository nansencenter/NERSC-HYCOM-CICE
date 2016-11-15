[TOC]

# Description

This repository contains the source code necessary to run the coupled HYCOM-CICE model developed at NERSC. In the repository you will find the HYCOM and CICE source codes, utility routines, and a directory containing a simple model setup (Directory TP0a1.00). Below is a short description of how to set up and compile the code

# Prerequisites

The following tools are necessary 

* A working fortran compiler.

* To run the HYCOM-CICE coupled code, you will need to have a working installation of the Eearth System Modelling Framework: https://www.earthsystemcog.org/projects/esmf/download/. The code has been tested and verified to work with ESMF v 5.2.0.rp3.

* A working python 2.7 installation, [more info can be found here](doc/python.md)


# Retrieving HYCOM-CICE
Clone code from main repository (TODO: Fix when moved/use markdown)

`
git clone https://bitbucket.org/knutal/nersc-hycom-cice
`

If you get errors about server certificates, see [here](../..//overview#markdown-header-server-certificates)

# Compiling HYCOM-CICE support programs

* Compile MSCPROGS routine, following the instructions given in [hycom/MSCPROGS/](hycom/MSCPROGS/)

* Compile HYCOM_ALL programs, following the instructions given in [hycom/HYCOM_ALL/](hycom/HYCOM_ALL/)

# Compiling HYCOM-CICE

Before yout attempt this, make sure that the HYCOM-CICE support programs are compiled, and that the MSCPROGS programs are compiled.
