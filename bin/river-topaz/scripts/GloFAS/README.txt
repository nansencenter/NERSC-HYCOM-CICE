How to compute the river forcing .[ab] file from GloFAS data:

1) -----------------------Necessary python env---------------------


-Python (version 3.10 or above)
-Numpy (version 2.1.3 or above)
-Pandas (version 2.2.3)
-xarray (version 2024.11.0)
-scipy (version 1.15.2)
-scikit-image (version 0.25.2)

PS: you can also use the environment.yml file in this folder if it works

2) --------------------Necessary files and parameters-------------------

*First, go to the file compute_discharge_NoUI.py, the path to certain files must be given:

-topaz_grid_file: is the path to regional grid .a file for topaz_grid_file

-topaz_depth_file: is the path to depth grid .a file for topaz_grid_file

-ldd_file: is the path to GloFAS LDD file (auxiliary file). The GloFAS auxiliary files can be downloaded with the discharge files. Be careful that auxiliary files must be from the same GloFAS version as discharge file

-uparea_file: path to upstream area GloFAS auxiliary file

-elevation_file: path to elevation GloFAS auxiliary file

-river_correction_file: Netcdf file with correction paramters from GloFAS. Has to be computed first using program in the folder called GloFAS_correction_process.

-river_correction_file: .csv file containing potential modifications of estuary postion. If you don't plan on doing any manual modification, don't worry about it. Otherwise, head to the technical documentation for more information on how to compute this file


*Now go to the file river_forcing.py and give the followinf paths:

-input_dataset: path to NetCDF file containing GloFAS discharge. This file has to to be on the machine you're using (not on NIRD). And it needs to contain dsicharge for every day in the selected period; you may need to aggregate files from certain years or use a file that contains all years already. Be careful with the format, all the fields that are originally in the GloFAs discharge file must be here with the same name, otherwise you'll have issues. If it's still there, you can use "/nird/datalake/NS9481K/Adrien/GloFAS/data_all_year.nc", it has discharge data aggregated from 1979-2022 for GloFAS v4.0 

-output_path: the path where you want the output forcing file to be

*Still in river_forcing.py, you must set the following parameters:

-start_date: Starting date of river forcing file with format: YYYY-MM-DDT00:00:00.000000000

-end_date: End date of river forcing file with format: YYYY-MM-DDT00:00:00.000000000

-Apply_GloFAS_correction: Boolean, True to Apply GloFAS correction starting from 2001

-Edit_estuaries: Boolean, True to take into account estuary edit file and modify the positions for certain estuaries

-Propagation_cleaning_step: Boolean, True to prevent river plume to extend beyond land (can create discontinuities in the freshwater plume)

-rradius: Radius for river plume propagation based on distance to shore [m]

-alongshoreradius: Radius for river plume propagation based on distance to estuary [m]


3) --------------------Compute the file-------------------

Execute river_forcing.py and hope for the best. The output forcing .[ab] file will be in the output folder you indicated. It will be named forcing_{start}_{end}_6hourly.[ab]

