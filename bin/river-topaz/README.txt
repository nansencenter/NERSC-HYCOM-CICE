Welcome to the TOPAZ_runoff program. This program is designed to compute river forcing for TOPAZ using GloFAS data. 
Here is a quick rundown of the folders you will find here:

- "scripts" contains the main scripts to compute the river forcing files. If you're here to compute the river forcing file from GloFAS data. Go to scripts/GloFAS/README.txt

- "GloFAS_correction_process" contains scripts and instructions to compute the GloFAS correction file.

- "data" contains sample data necessary to compute the river forcing. Most of the files are there. Head to data/README.txt to see where you can find the file you need.

- "outputs" contains examples of output files from the program to edit the estuary position

Finally, you may want to have a look at the technical documentation named 'river_forcing_documentation.pdf"
#  First developed by Adrien Acchiardi at Augest 2025 
#---------------------------------------


# After November 2025, the Python routine has been merged into the GitHub of HYCOM_CICE 
#-------------------------
1) The key Python routine is bin/river_synoptic.sh. It can be performed in the routine of expt_preprocess.sh with online or offline by setting River_online.
   
2) In river_synoptic.sh, the output directory for the river files is automatically created as riverh.ab under the forcing/river/$expnum/ when using River_online=1. In some cases, if the time period is more than one year, it will take a long time. In order to save computing resources, performing the command in your experiment directory, but don't forget to switch off River_online in the expt_preprocess.sh 

Example:
    bin/river_synoptic.sh 2013-01-01T00:00:00  2014-01-05T00:00:00
 
3) The source data including the model grid information and the river data (netCDF) are stored at /cluster/projects/nn9481k/GloFAS_data/data_v40 and /cluster/projects/nn9481k/GloFAS_data/TOPAZrunoff_data/. They are defined in the Python routine:river-topaz/scripts/GloFAS/hycom_river_hifre.py


4) N.B.: (hycom.2.2.98): before your model runs, the two parameter settings should be consistent: priver=0 (blkdat.input) and highfq_river  = .true. (hycom.opt). Under this setting, the high-frequency river will be active in the model. But in hycom 2.3, the "priver" setting was changed in recent months. So the relevant setting of priver requires it to be set to 2 for the high-frequency river forcing, and skipping the highfq_river=.true. setting. in hycom_opt.

Limits: At moment, all the river data are stored in one netCDF file data_all_year.nc, but it covers the years from 1979 to the end of 2024. For real-time implementation, this file should be replaced, and depending on the practical implementation case.
# updated by Jiping Xie in November 2025 
#---------------------------------------


