This directory contains routines for processing and plotting model output.


# Tools in this directory

|source file     | purpose|
|-------- | -------------|
|hycom_convert_daily2monthly_means.sh| compute monthly mean and the mean squire from daily archm files |
|interpolate_ostia_SST_to_TOPAZgrid.py| quick interpolation from OSTIA SST to TOPAZ grid; it takes as input the OSTIA filename and outputs the interpolated SST in a netCDF file |

# Plotting in this directory 

|executable     | purpose|
|-------- | -------------|
|hycom_plot_section_topaz.py  | creating section plot from ab files with an option for density plot    |
|hycom_plot_topaz_ts_from_ncfile.py | plotting  time series of spatially-averaged of 2D field from netCDF files   |
|hycom_plot_topaz_ts_from_abfile.py | plotting  time series of spatially-averaged of 2D field from ab files     |
|hycom_plot_archive_topaz_Avg |      ploting 2D  from a given ab file, and it will plot the averaged if the input are ab muliple files|


# Othe plotting tools in NERSC-HYCOM-CICE/bin

|executable     | purpose|
|-------- | -------------|
|hycom_plot_archive.py  | 2D plot  of fields from archives. Has many options to fine-tune plots |
|hycom_plot_field.py  | 2D plot  of fields from any .a file. Requires grid size as input |
|hycom_plot_section.py | section plot  of fields from archives. Has many options to fine-tune plots |
|cice_icevolume.py| Calculates ice volume from a set of CICE netcdf files. Input is actually file pattern strings interpreted by python glob module. If you dont know what this means, its probably a bit difficult to use ...|

