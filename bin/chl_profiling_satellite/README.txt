CAGLAR: Nov22
how to use:

_FORTRAN_subroutine_

The code uses a fortran loop for a higly efficient loop for mapping satellite to model grid.
To use the fortran subroutine, notice that at the top, we imported machine specific library (e.g. _BETZY_the_mapping_loop_with_ice)
You need to import your own machine library. To create it (inside the bin/chl_profiling_satellite/ folder) specify the MACHINE as your own

>>  f2py -c -m _MACHINE_the_mapping_loop_with_ice satellite_map.f90

This should create .so file  with a long name, rename it to '_MACHINE_the_mapping_loop_with_ice.so' format 

You can then add this machine specific library to the code. Search in the code 'gethostname' and copy/modify the if statement following the others. 
____________________

The code works with mandatory inputs. 
You need to specify these:

domain, experiment, year and day

Also consider the following options:

--nerscdir=/cluster/home/cagyum/NERSC-HYCOM-CICE
--satdir=/cluster/work/users/cagyum/OCCCI/L3/copernicus_globcolor_daily_rep/chl/

--opendap_rep # this overwrites user satdir and gets the Copernicus L3 reprocessed data through opendap
              # you need to have a Copernicus account.
              # the code will use your username and password.
              # you can provide your user and password as a system variable 
              # e.g. in bash: export copernicus_user=XXX  && export copernicus_pass=XXX  

--opendap_rep # this overwrites user satdir and gets the Copernicus L3 NRT data through opendap
              # you need to have a Copernicus account.
              # the code will use your username and password.
              # you can provide your user and password as a system variable 
              # e.g. in bash: export copernicus_user=XXX  && export copernicus_pass=XXX

--debug       # save a netcdf file containing satellite, old chl and new chl.



# The script has been tested with the following inside the experiment folder with a link to bin folder in the upper directory:

python ../bin/chl_profiling_satellite/chl_profile.py /cluster/work/users/cagyum/TP5a0.06 020 $y $nn --nerscdir=/cluster/home/cagyum/NERSC-HYCOM-CICE --opendap_nrt --debug
