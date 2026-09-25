The scripts in this folder can compute the parameters that are used to correct GloFAS discharge data for the 3 largest rivers in Russia (Yenisei, Lena and Ob). 


1) ----------------------------------Necessary python env---------------------


-Python (version 3.10 or above)
-Numpy (version 2.1.3 or above)
-Pandas (version 2.2.3)
-xarray (version 2024.11.0)
-scipy (version 1.15.2)

2) ------------------------------Download the necessary files:------------------

-Download ArcticGRO data as csv files from https://arcticgreatrivers.org/data/ for the 3 rivers (Yenisei, Lena and Ob).

-Download yearly GloFAS discharge NetCDF files (one file per year). We need the all the files from the years that will be used in the correction parameters computation. The names of the files must be 'data_yyyy.nc'

3)----------------------------Compute the ArcticGRO Dataframes------------------------

*Go to the file compute_ArcticGRO_Dataframes.py

-In the file compute_ArcticGRO_Dataframes.py, give the path to ArcticGRO csv files in "path"

-Run the program, the dataframes will be put in the folder 'Dataframes' (amke sure that this folder exists in same folder as the python files)

-Run file compute_ArcticGRO_Dataframes.py to obtain Pandas dataframes with ArticGRO data


4) --------------------------------Compute correction parameters----------------------------------

*Go to the file compute_largest_river_correction.py

-Check the paths 'path_GRO_Dataframes' that must lead to the ArtcicGRO dataframes previously computed and 'path_GloFAS_NC_files' that must lead to GloFAS data used for correction (the names of the files must be 'data_yyyy.nc' for each year used)

-Chose the 'start_year' and 'end_year' that you want to use to compute correction (the corresponding GloFAS discharge files must be downloaded).

-Run file compute_largest_river_correction.py to compute correction paramters. The output will be a NetCDF file that will be used directly in the river discharge main code.