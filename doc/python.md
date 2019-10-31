[TOC]

# Prerequisites

You will need a  working python 2.7 installation, with the following packages.

* [numpy](https://pypi.python.org/pypi/numpy)

* [scipy](https://pypi.python.org/pypi/scipy)

* [pyproj](https://pypi.python.org/pypi/pyproj)

* [f90nml](https://pypi.python.org/pypi/f90nml)

* [cfunits](https://pypi.python.org/pypi/cfunits)

* [basemap](https://pypi.python.org/pypi/basemap)

* [matplotlib](https://pypi.python.org/pypi/matplotlib)

* [netCDF4](https://pypi.python.org/pypi/netCDF4)

These depend on other non-python packages being installed, such as udunits and netcdf4. Most of these packages are usually installed on a linux system. If they are missing, you can ask an IT guy to install them on the system  or install them yourself.

In addition, these packages are required (developed by Knut and available on github):

* [gridxsec](https://github.com/knutalnersc/gridxsec), a tool for creating cross-sections on 2D grids 

* [abfile](https://github.com/knutalnersc/abfile), a tool for reading HYCOM .[ab] files 

* [modelgrid](https://github.com/knutalnersc/modelgrid), a tool for creating modelgrids using Bentsen et al conformal mapping or standard projections 

* [modeltools](https://github.com/knutalnersc/modeltools), a collection of various tools .. 

* [modeltools](https://github.com/MostafaBakhoda/modeltools), a collection of various tools, few changes done by Mostafa ..

# Checking for missing python modules

When running the HYCOM-CICE python routines you will usually get an error message if a module is missing. But you can also check more directly.

A quick way to check is to start python (or ipython) and write import <name_of_module>. The following shows the error message if a module is missing

    import numpy
    ---------------------------------------------------------------------------
    ImportError                               Traceback (most recent call last)
    <ipython-input-2-5a8446063dec> in <module>()
    ----> 1 import numpy
    ImportError: No module named numpy

If you find some packages are missing have an IT guy install them on the system you are working on, or install them yourself (some tips below). 

# Installation of missing packages

There are many ways of installing python packages, you can install binary packages (.rpm, .deb etc) for your distribution, you can download source code and compile and install locally, or use python installation tool "pip".
If you need to install a package locally, try to use pip, as it is usually the easiest. The following will install f90nml under your local user account:

    knutal@debian:~$ pip install --user f90nml
    Downloading/unpacking f90nml
      Downloading f90nml-0.19.tar.gz
      Running setup.py (path:/tmp/pip-build-2D2rgF/f90nml/setup.py) egg_info for package f90nml
    
    Installing collected packages: f90nml
      Running setup.py install for f90nml
    
    Successfully installed f90nml
    Cleaning up...
    knutal@debian:~$ ls -altr

Note that some machines (ex hexagon) have issues  with SSL certificates. If you get errors when running the command above (or if "pip search f90nml" returns errors about SSL certificates), you can try this command in stead:

    pip install --user  --index-url=http://pypi.python.org/simple/ --trusted-host pypi.python.org f90nml

These files will be installed in your home directory under $HOME/.local/lib/python2.7/site-packages/, where python will automatically find them.





# Installing the python modules from github

The github packages can also be installed using the pip tool.  Use the following command to install them:

    pip install --user git+http://github.com/knutalnersc/abfile
    pip install --user git+http://github.com/knutalnersc/modeltools
    pip install --user git+http://github.com/knutalnersc/modelgrid
    pip install --user git+http://github.com/knutalnersc/gridxsec
    pip install --user git+http://github.com/MostafaBakhoda/modeltools

We have currently moved all python tools from github to our NERSC-HYCOM-CICE repository (pythonlibs) and replace "path-to-NERSC-HYCOM-CICE":

   pip install --user file:path-to-NERSC-HYCOM-CICE/pythonlibs/modeltools
   pip install --user file:path-to-NERSC-HYCOM-CICE/pythonlibs/modelgrid
   pip install --user file:path-to-NERSC-HYCOM-CICE/pythonlibs/gridxsec 
   pip install --user file:path-to-NERSC-HYCOM-CICE/pythonlibs/abfile 

To upgrade the libraries, just call the same install commands with the extra option: "--upgrade"

These files will be installed in your home directory under $HOME/.local/lib/python2.7/site-packages/, where python will automatically find them.
