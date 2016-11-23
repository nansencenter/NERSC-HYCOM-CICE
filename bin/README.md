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



# Experiment setup, configuration, maintenance, backup 

|executable     | purpose|
|-------- | -------------|
|region_new.sh  | Sets up a skeleton directory of a region directory, based on another region directory |
|region_backup.sh | Backs up files for a region (grid and topo), either to local path or to norstore |
|expt_new.sh  | Creates a new experiment directory , based on an existing experiment directory |
|expt_postprocess.sh  | Transfers data from the experiment scratch directory to its data directory, for safe(r) keeping |
|expt_preprocess.sh  | Trransfers necessary input files to scratch in preparation of a model simulation |
|expt_config_backup.sh  | Backs up a experiment configuration (basic input files, blkdat, ice-in etc) to a local path or to norstore |
|expt_transfer_norstore_scratch.sh | Transfers data from an experiment directory to scratch on norstore (actual storage on norstore must be done after this |
|cleanSCRATCH.sh | removes files in every directory named SCRATCH in the experiment or region directory (run wo args will not delete anything, just describe what is going to be deleted)|
|tar_in_chunks.py | Create tarfiles from files. The tarfiles have a fixed size (specified on input). Meant for use on norstore |

# Plotting

|executable     | purpose|
|-------- | -------------|
|hycom_plot_archive.py  | 2D plot  of fields from archives. Has many options to fine-tune plots |
|hycom_plot_field.py  | 2D plot  of fields from any .a file. Requires grid size as input |
|hycom_plot_section.py | section plot  of fields from archives. Has many options to fine-tune plots |




# Generation of input data


# Tools


# Nesting/interpolation tools
|isuba_gmapi.sh        | Create a horizontal mapping from one region to another |
|isuba_topog.sh        | Use mapping (see above) to interpolate topo file from one region to another |
|isubaregion_archv.sh  | Use mapping (see above) to interpolate archive from one region to another |
|remap_archv.sh        | Do vertical remapping of archive files |


Notebooks
atmo_nersc_clim.sh
atmo_nersc_synoptic_oldversions.sh
atmo_synoptic.sh
calc_thkdf4.sh
cice_icevolume.py
cice_kmt.py
cice_kmt.sh
cice_limits.py
cleanEXP.sh
compile_model.sh
hycom_atmfor.py
hycom_bathy.py
hycom_bathy_consistency.py
hycom_bathy_modify.py
hycom_date.py
hycom_grid.py
hycom_kapref.py
hycom_limits.py
hycom_sigma_compare.py
hycom_sigma_setup.py
hycom_topo_merge.py
hycom_topo_ports.py
hycom_vcoord_plot.py
hycom_vsigma.py
hycom_woa13.py
hycom_woa13_zfiles.py
hycomsteps.py
namelist_extract.py
nemo_mesh_to_hycom.py
nemo_to_hycom.py
nest_check_ports.sh
nest_setup_ports.sh
old
relax_rmu.sh
relaxi.sh
river_nersc.sh
river_nersc_bio.sh
river_trip.sh
scripts_specific
seawifs_mon_kpar.sh
setlimits.py
soda_to_hycom.py
tile_grid.sh
woa2013_density_distribution.py
woa2013_sigma_profile.py
z_generic.sh
