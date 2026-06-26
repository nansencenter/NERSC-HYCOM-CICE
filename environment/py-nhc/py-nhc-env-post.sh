#!/usr/bin/bash

# NB: Set this path!
export PYLIBSPATH=/path/to/NERSC-HYCOM-CICE/pythonlibs

pip install $PYLIBSPATH/modeltools
pip install $PYLIBSPATH/modelgrid
pip install $PYLIBSPATH/gridxsec
pip install $PYLIBSPATH/abfile
