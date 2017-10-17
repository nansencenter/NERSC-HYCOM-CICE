# Retrieving and installing ESMF

To run the HYCOM-CICE coupled code, you will need to have a working installation
of the Eearth System Modelling Framework:
https://www.earthsystemcog.org/projects/esmf/download/. The code has been tested
and verified to work with ESMF v 5.2.0.rp3. and ESMF v 6.3.0.

You could try to install his yourself, instructions are available in the downloaded distribution. 
However, it is probably recommended that you let your local IT support to do it
to make sure that MPI etc is properly set up.

Once installed, hycomc-cice can be compiled - [more on hycom compilation here](./HYCOM-CICE-compilation.md).

# Site-specific details

## Hexagon.bccs.uib.no

On hexagon, ESMF version 6.3.0rp1 has been tested for the pgi compilers. In
order to use it, issue the command

module load esmf/6.3.0rp1-pgi

In order to use this library, hycom-cice will have to be compiled with the pgi
compiler. 

## Sisu.csc.fi
## Generic linux machine with fortran

On debian, hycom was successfully compiled and tested  using the gfortran compiler and the
openmpi MPI library. The gfortran compiler and openmpi should be installed
first, then the ESMF library (which should use openmpi).

