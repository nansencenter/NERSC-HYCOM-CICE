# Description

# Prerequisites

The following tools are necessary in many parts of the code

# Prerequisites
A Fortran compiler

To run the HYCOM-CICE coupled code, you will need to have a working installation of the Eearth System Modelling Framework.

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

Clone code from main repository (TODO: Fix when moved/use markdown)
`
git clone https://bitbucket.org/knutal/nersc-hycom-cice
`

If you get an error like "server certificate verification failed", you will need to install certificates on the machine where you want to run the model(or contact IT support). More [here...](https://en.wikipedia.org/wiki/Certificate_authority). If certificate installation fails, you can try this as a last resort:
`
GIT_SSL_NO_VERIFY=true git clone https://bitbucket.org/knutal/nersc-hycom-cice
`

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
    cd where_i_want_to_install_python_modules


# Site-specific details