This directory contains routines for processing and postprocessing model output.
The shell scripts (.sh) mainly utilize programs in the HYCOM_ALL and MSCPROGS,
and some ofd the python routines in this directory. Meny of the shell scripts
require you to run them from an experiment directory, because it needs to source
the REGION.src and EXPT.src files, copy files to the correct location etc etc.

Most of the python scripts can be called from the command line, but some need
the hycom setup files (grid and bathymetry) to be present.

# Common functions 



|source file     | purpose|
|-------- | -------------|
|common_functions.sh  | functions used by many of the experiment processing scripts|


# Tools in the "useful" category

|source file     | purpose|
|-------- | -------------|
|hycom_date.py| Can be used to convert to/from actual date/time to time used by hycom and ordinal day used by hycom|




# Experiment setup, configuration, maintenance, backup 


|executable     | purpose|
|-------- | -------------|
|region_new.sh  | Sets up a skeleton directory of a region directory, based on another region directory |
|region_backup.sh | Backs up files for a region (grid and topo), either to local path or to norstore |
|expt_new.sh  | Creates a new experiment directory , based on an existing experiment directory |
|expt_config_backup.sh  | Backs up a experiment configuration (basic input files, blkdat, ice-in etc) to a local path or to norstore |
|expt_transfer_norstore_scratch.sh | Transfers data from an experiment directory to scratch on norstore (actual storage on norstore must be done after this |
|cleanSCRATCH.sh | removes files in every directory named SCRATCH in the experiment or region directory (run wo args will not delete anything, just describe what is going to be deleted)|
|tar_in_chunks.py | Create tarfiles from files. The tarfiles have a fixed size (specified on input). Meant for use on norstore |
|hycom_limits.py | Calculates HYCOM time limits from actual input dates. mainly used by expt_preprocess.sh when setting up model |
|expt_postprocess.sh  | Transfers data from the experiment scratch directory to its data directory, for safe(r) keeping |
|expt_preprocess.sh  | Trransfers necessary input files to scratch in preparation of a model simulation |
|compile_model.sh | Script for compiling the HYCOM-CICE model. |
|tile_grid.sh | Sets up the MPI partitioning of the hycom model. |
|cice_limits.py | Create a CICE input file based on existing file and some arguments. Used to set up correct time steps etc |
|hycom_sigma_setup.py| Can be used to set up sigma coordinates as a set of linear profiles stitched together|
|hycomsteps.py| Can be used to calculate suitable time steps for hycom, based on some hycom restrictions|


# Plotting

|executable     | purpose|
|-------- | -------------|
|hycom_plot_archive.py  | 2D plot  of fields from archives. Has many options to fine-tune plots |
|hycom_plot_field.py  | 2D plot  of fields from any .a file. Requires grid size as input |
|hycom_plot_section.py | section plot  of fields from archives. Has many options to fine-tune plots |
|cice_icevolume.py| Calculates ice volume from a set of CICE netcdf files. Input is actually file pattern strings interpreted by python glob module. If you dont know what this means, its probably a bit difficult to use ...|



# Generation of grids


|executable     | purpose|
|-------- | -------------|
|hycom_grid.py       | Create bathymetry based on conformal mapping or on proj4 projection. Can do  the same as old confmap routines in MSCPROGS |



# Generation of bathymetries

|executable     | purpose|
|-------- | -------------|
|hycom_bathy_merge.py       | Merges two bathymetries on the same grid. |
|hycom_bathy_consistency.py | Does a consistency check of a hycom bathymetri (removes isolated basins, single cell islands, etc)|
|hycom_bathy.py             | Generate hycom bathymetry for a predefined hycom grid |
|cice_kmt.py | Simple script for creating CICE land mask from hycom bathymetry |
| cice_kmt.sh | Basically a wrapper around coce_kmt.py |
|hycom_bathy_modify.py| Can be used to modify bathymetry in points or in blocks based on input from stdin|


# Generation of relaxation fields

|executable     | purpose|
|-------- | -------------|
|z_generic.sh  | Generates a climatology on fixed z levels for the current experiment domain |
|relaxi.sh  | Generates a climatology in hybrid coordinates forr the current experiment. Needs z_generic to be run first |
|relax_rmu.sh  | Generates sidewall relaxation weights |
|hycom_woa13_zfiles.py  | Generates .d files from WOA2013 data. These can be used by routine z_generic.sh |
|hycom_woa13.py  | Generates gridded z-level  files from WOA2013 data in ab format. Output similar to z_generic.sh (this routine can replace z_generic.sh) |
| soda_to_hycom.py | Starting poinr for reading SODA data and interpolating them to hycom for nesting. Not finished!| 

# Generation of various forcing input


|executable     | purpose|
|-------- | -------------|
|calc_thkdf4.sh| Generates  hardcoded  diffusion fields for hycom   (needed if relevant flag in blkdat.input is -1). NB: hardcoded to increase  diffusiuon for Brazil/Gluf Stream. NOT generic ...|
|hycom_kapref.py| Can be used to create spatially varying thermobaric reference density. The thermobaric parameterization tends to be unstable, so maybe you shouldnt use thios .... |
|seawifs_mon_kpar.sh| Can be used to set up Jerlov water type based on a seawifs climatology provided for hycom.|
|hycom_vsigma.py| Statring point for generating spatially varying target densities (vsigma in blkdat.input) NOT finished !|




# Generation of atmospheric forcing

|executable     | purpose|
|-------- | -------------|
|hycom_atmfor.py | generates atmospheric forcing, based on xml file input |
|atmo_synoptic.sh | Wrapper around hycom_atmfor.py |

# Generation of river forcing

|executable     | purpose|
|-------- | -------------|
| river_nersc.sh | Use simple nersc river routine |
| river_trip.sh | River fluxes based on TRIP and ERAI or ERA40 climatology. Can now use pre-existing riverflow in TRIP_PATH |
| river_trip_realtime.sh | River fluxes based on TRIP, can generate realtime river forcing based on ERA40 or ERAI as well as river climatology|



# Nesting/interpolation tools

|executable     | purpose|
|-------- | -------------|
|isuba_gmapi.sh        | Create a horizontal mapping from one region to another |
|isuba_topog.sh        | Use mapping (see above) to interpolate topo file from one region to another |
|isubaregion_archv.sh  | Use mapping (see above) to interpolate archive from one region to another |
|remap_archv.sh        | Do vertical remapping of archive files |
|hycom_topo_ports.py   | Generate a ports.input file  and rmu file suitable for use when nesting |
| nest_setup_ports.sh  | Generates a ports.input file and rmu file suitable for nesting for an experiment (wrapper around hycom_topo_ports.py)|
| nest_check_ports.sh  |Checks a port setup for errors |
| nemo_mesh_to_hycom.py | Converts NEMO grid to hycom file types (regional.grid and regional.depths)|
| nemo_to_hycom.py  | Converts NEMO input files  to hycom archv files  (also calls nemo_mesh_to_hycom)|


# Tide tools 

|executable     | purpose|
|-------- | -------------|
|hycom_add_tides.py        | Used to tidal data to a existing hycom archv file or create archive files suitable for barotropic nesting only |
|tidal_forcing.py       | Generates tidal data from FES2004 or FES2014, then runs hycom_add_tides.py to generate archive files with tidal data at boundaries |
|tidal_forcing.sh       | Wrapper around tidal_forcing.py - this is usually the routine you want to run |


# Diagnostic tools 


|executable     | purpose|
|-------- | -------------|
|hycom_vcoord_plot.py          | Create plots showing the layout of the vertical coordinate specified by blkdat.input|
|woa2013_sigma_profile.py      | Creates a sigma profile from woa2013 data. Tries to set up a vertical coordinate matching a blkdat.input setup |
|woa2013_density_distribution.py | Calculates density distribution ++ for a specified region. Based on WOA2013 data |
|hycom_sigma_compare.py| Can be used to show the difference in hybrid layers between setups. Input is different blkdat.input files |


# Other tools

|executable     | purpose|
|-------- | -------------|
|hycom_montg1.py| Can be used to calculate montg1 for archive files, based on 3D profile in archive and other sources|
|namelist_extract.py| Routine for extracting info from FORTRAN namelist. Mainly used to process CICE namelists|


# The rest

|path     | purpose|
|-------- | -------------|
|Notebooks| Contains some notebooks that can be used together with ipython|
|scripts_specific| some scripts that are not generic ...|
|not_implemented| Some scripts that are not fully implemented ...|
