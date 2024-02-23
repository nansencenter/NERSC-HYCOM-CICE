Receipt of how to generate nesting conditions from an ESM:

	1. Set up the directory representing the outer model.
		- Create a directory called ESMa1.00 and inside:
			* Make a link to bin in you HYCOM-code
			* Make a directory called expt_01.0
			* Make a directory called topo
			* Copy REGION.src to this directory
			* Get the grid information from the model in questions, this can be found on the CMIP6 server: you need the file for the varible "areacello" and "deptho" and put in the same folder.
			○ From the experiment folder run script:
			../bin/Nesting_noresm/make_grid_noresm2hycom.py path_to_grid_files/areacello_file.nc 
			This will generate the regional.depth and regional.grid-files for the global mode grid.
			Place these in the topo-folder
			○ From the experiment folder run script:
			○  ../bin/isuba_gmapi.sh. $target_region to genreate the mapping from the global model to the regional hycom.
		
		
	2. Download the nesting-files needed: for physics: uo, vo, zos, thetao, so, these are annual files.
		a. In order to have files that have values at all ocean points from the coarse ESM-grid, the annual files are pre-processed before generating the boundary conditions.
../bin/Nesting_noresm/separate_and_extrapolate_files_year.sh  $varible $year, can also be called on an interactive node using  /bin/Nesting_noresm/Generate_nesting_files_year.sh
