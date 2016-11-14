Prep:

1. netcdf4 [not working - need correct modules]
A. cp makefile_nc4 makefile
B. source modules_nc4.src

2. netcdf3 [default modules work]
A. cp makefile_nc3 makefile

Compilation:
1. Edit p_getpivots.F90: location of example data files, data grid size etc, names of lon/lat...
2. Edit name of executable in Makefile
3. Compile with "make"
!. Used to need to set modules with "source modules.src"
   due to netcdf4 format (never got this working),
   but now convert files to netcdf3 at download time

Running and testing:
cd RUN_HERE

Running:
1. Edit get_pivots.sh
2. ./get_pivots.sh

Testing:
1. Edit test_pivots.m or test_pivots.py
2. Run preferred script [m_proj not working currently for matlab]

To use in HYCOM:
1. Set link in preprocess_waves.sh
2. Make a CPP flag
3. Edit mod_readwaves.F to make read_pivots_[data] and read_waves_[data] subroutines
