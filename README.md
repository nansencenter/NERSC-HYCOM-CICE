[TOC]

# Description

# Prerequisites

The following tools are necessary in many parts of the code.

A working fortran compiler.

To run the HYCOM-CICE coupled code, you will need to have a working installation of the Eearth System Modelling Framework: https://www.earthsystemcog.org/projects/esmf/download/. The code has been tested and verified to work with ESMF v 5.2.0.rp3.

A working python 2.7 installation, with the following packages

* [numpy](https://pypi.python.org/pypi/numpy)

* [scipy](https://pypi.python.org/pypi/scipy)

* [pyproj](https://pypi.python.org/pypi/pyproj)

* [f90nml](https://pypi.python.org/pypi/f90nml)

* [cfunits](https://pypi.python.org/pypi/cfunits)

* [basemap](https://pypi.python.org/pypi/basemap)

* [matplotlib](https://pypi.python.org/pypi/matplotlib)
* [netCDF4](https://pypi.python.org/pypi/netCDF4)

Additional python packages may be required, but these are usually installed in most python distributions. In addition, these packages are required (developed by Knut and available on github):

* [gridxsec](https://github.com/knutalnersc/gridxsec), a tool for creating cross-sections on 2D grids 

* [abfile](https://github.com/knutalnersc/abfile), a tool for reading HYCOM .[ab] files 

* [modelgrid](https://github.com/knutalnersc/modelgrid), a tool for creating modelgrids using Bentsen et al conformal mapping or standard projections 

* [modeltools](https://github.com/knutalnersc/modeltools), a collection of various tools .. 

# Retrieval and compilation of code

## Retrieving NERSC-HYCOM-CICE
Clone code from main repository (TODO: Fix when moved/use markdown)
`
git clone https://bitbucket.org/knutal/nersc-hycom-cice
`

If you get errors about server certificates, see [here](../..//overview#markdown-header-server-certificates)


## Retrieving and installing python modules

Next you will need to install the python packages necessary to use the  python routines. With some luck, most of the necessary packages are already installed.  A quick way to check is to start ipython and write import <name_of_package>. Example:

    import numpy
    ---------------------------------------------------------------------------
    ImportError                               Traceback (most recent call last)
    <ipython-input-2-5a8446063dec> in <module>()
    ----> 1 import numpy
    ImportError: No module named numpy

If you find some packages are missing, either install them yourself, or have your friendly IT guy install them for you. The following packages can be installed by root in the python distribution, and are available in python installation tools (pip). Some of them are also available as binary packages in linux (RPM, debian package etc etc):
[numpy](https://pypi.python.org/pypi/numpy)  
[scipy](https://pypi.python.org/pypi/scipy) 
[pyproj](https://pypi.python.org/pypi/pyproj) 
[f90nml](https://pypi.python.org/pypi/f90nml) 
[cfunits](https://pypi.python.org/pypi/cfunits) 
[basemap](https://pypi.python.org/pypi/basemap) 
[matplotlib](https://pypi.python.org/pypi/matplotlib) 
[netCDF4](https://pypi.python.org/pypi/netCDF4). It is also possible to install these packages outside the main python distribution tree, but the process is a bit more elaborate.

The github packages are not yet available as installers, and can be installed as follows:

    cd [location_of_my_python_modules]
    git clone https://github.com/knutalnersc/gridxsec
    git clone https://github.com/knutalnersc/abfile
    git clone https://github.com/knutalnersc/modelgrid
    git clone https://github.com/knutalnersc/modeltools

Replace [location_of_my_python_modules] with the actual path. If you get errors about server certificates, see [here](../..//overview#markdown-header-server-certificates). 

In order to use these locally installed packages you will have to set PYTHONPATH in your shell setup file (.profile / .bash_profile ):

    export PYTHONPATH=$PYTHONPATH:[location_of_my_python_modules]

Again, replace [location_of_my_python_modules] with the actual path

# Retrieving and installing ESMF

To run the HYCOM-CICE coupled code, you will need to have a working installation of the Eearth System Modelling Framework: https://www.earthsystemcog.org/projects/esmf/download/. The code has been tested and verified to work with ESMF v 5.2.0.rp3.

You could try to install his yourself, instructions are available in the downloaded distribution. However, it is probably recommended that you let your local IT support to do it.


# Site-specific details

## Hexagon.bccs.uib.no

On hexagon, most tools are installed already. Make sure these commands are run first:

    module load python/2.7.9-dso
    export PYTHONPATH=$PYTHONPATH:/home/nersc/knutali/opt/python/lib/python2.7/site-packages/

Apart from that, make sure you use the pgi compiler and libraries in the cray PrgEnv.


# This and that...

## Server certificates
If you get an error like "server certificate verification failed", you will need to install certificates on the machine where you want to run the model(or contact IT support). More [here...](https://en.wikipedia.org/wiki/Certificate_authority). If certificate installation fails, you can try this as a last resort before issuing the git clone commands:
`
export GIT_SSL_NO_VERIFY=true
`