#### Experiment folder for HYCOM-OASIS-NeXtSIM
This directory is an example of how to run HYCOM coupled to NeXtSIM.
The run rquires a compiled version of NeXtSIM and a compiled version of hycom.

Files:\
**OASIS_input/ice.fabm.nc**\
OASIS restart file for the fields sent by the ice model when fabm is activated.\
\
**OASIS_input/ice.nc**\
OASIS restart file for the fields sent by the ice model when running with only physics.\
\
**OASIS_input/ocean.nc**\
OASIS restart file for the fields sent by the ocean model.\
\
**OASIS_input/namcouple.fabm**\
OASIS coupling specifications when fabm is activated.
\
**OASIS_input/namcouple**\
OASIS coupling specifications when running with only physics.\
\

**blkdat.input**: this must be put in the experiment folder where hycom is compiled.  \
To compile run coupled to nextsim, the following must be setin blkdat:\
```
   0   'iceflg' = sea ice model flag (0=none,1=energy loan,2=coupled/esmf)\
   7   'wndflg' = wind stress input flag (0=none,1=u/v-grid,2,3=p-grid,4,5=wnd10m,7=nextsim)\
   4   'ustflg' = ustar forcing     flag        (3=input,1,2=wndspd,4=stress)\
   3   'flxflg' = thermal forcing   flag (0=none,3=net_flux,1-2,4-6=sst-based)\
   3   'empflg' = E-P     forcing   flag (0=none,3=net_E-P, 1-2,4-6=sst-based_E)\
```
   
**cpl_run.cfg**\
Configuration file the Nextsim. \

**slurm_cr.sh**
Slurm submission script for the coupled run.\

The run is submitted with this command.\
>sbatch slurm_cr.sh cpl_run.cfg
