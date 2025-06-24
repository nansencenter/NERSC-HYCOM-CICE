# Receipt of how to generate nesting conditions from an ESM:

## Set up the directory representing the outer model.
- Create a directory called ESMa1.00 in your top work directory (e.g. in Fram/Betzy --> /cluster/work/users/$USER/ESMa1.00) and inside:
- Make a link to bin in you HYCOM-code
- Make a directory called expt_01.0
- Make a directory called topo
```
mkdir ESMa1.00
cd ESMa1.00/
ln -sf $HOME/NERSC-HYCOM-CICE/bin .
mkdir expt_01.0
mkdir topo
```
- Copy REGION.src to this directory and make sure to modify the following:
  - R=ESMa1.00 
  - NHCROOT points to the HYCOM-code in your home directory.
  - ESM_ID (e.g. ESM_ID=NorESM2-MM_historical_r1i1p1f1)
  - Nesting_Files_PATH (e.g. Nesting_Files_PATH=/cluster/work/users/$USER/$R/Nesting_files/)
  - ESM_NATIVE_MESH (e.g. ESM_NATIVE_MESH=/cluster/work/users/$USER/$R/topo/areacello_Ofx_NorESM2-MM_historical_r3i1p1f1_gn.nc)
```
cp $HOME/NERSC-HYCOM-CICE/input/REGION.src .
```
- Get the grid information from the model in questions, this can be found on the CMIP6 server: you need the file for the varible "areacello" and "deptho" and put in topo-folder. Find the spatial dimentions of these files as it is input to the routine below.
- From the experiment (e.g. expt_01.0) folder run script:
```
cd ./expt_01.0
../bin/Nesting_noresm/make_grid_noresm2hycom.py ../topo/areacello_file.nc ../topo/deptho_file.nc --idm <x-dimention> --jdm <y-dimention>
```
This will generate the regional.depth and regional.grid-files for the global mode grid, copy these and move everyting to the topo-folder:
```
cp dummped_depth_ESMa1.00_01.a depth_ESMa1.00_01.a
cp dummped_depth_ESMa1.00_01.b depth_ESMa1.00_01.b
cp dummped_regional.grid.a regional.grid.a
cp dummped_regional.grid.b regional.grid.b
mv dummped_* ../topo/
mv depth_* ../topo/
mv regional.* ../topo/
```

- From the ESMa1.00-folder run the following script to generate the mapping from the global model to the regional hycom:
```
./bin/isuba_gmapi.sh <path_to_regional_model_folder>
```
			
## Physical nesting conditions 
- Download the nesting-files needed: for physics: uo, vo, zos, thetao, so, make sure you get the files on a with varibles on z-levels (if the model is not a z-level model).
- If you have not yet created bias-correction files for salinity and temperature, go to the section *Generating files with model bias used for bias correction* at the bottom. If you are using NorESM, this step is already done and the files for bias correction can be found on Fram in /cluster/projects/nn9481k/NORESM_bias/.
- ESMs generally have coarse resolution that results in empty values near the land boundary in higher resolution models when interpolated. In order to have files that have values at all ocean points from the coarse ESM-grid, the annual files are extrapolated towards the land boundary before generating the boundary conditions. 

The following command is a generic one that does the extrapolation and creates the nesting archive files in the output nest folder in the target directory. It assumes that you have performed the bias correction as described later in the document.

```
cd <experiment folder in ESM directory>
../bin/Nesting_noresm/Generate_nesting_files_year.sh <year> <region, e.g. TP2a0.10> <experiment, e.g. 010>
```
- The command above already includes the following command. It is provided here as an example. No need to do it again. 
```
../bin/Nesting_noresm/esm_to_hycom.sh ../../TZ4a0.10/expt_01.1/ ../../NORESM_Nesting/thetao_Omon_NorESM2-MM_historical_r1i1p1f1_gr_195512_extrap.nc
```

		
## Biogeochamical nesting: Download the nesting-files needed: for physics: no3, po4, si, o2, these are decadal files.
- Download the nesting-files needed: for BGC: no3, po4, si, and o2, make sure you get the files on a with varibles on z-levels (if the model is not a z-level model).
- If you have not yet created bias-correction files for bgc variables, go to the section *Generating files with model bias used for bias correction* at the bottom. If you are using NorESM, this step is already done and the files for bias correction can be found on Fram in /cluster/projects/nn9481k/NORESM_bias/.
- ESMs generally have coarse resolution that results in empty values near the land boundary in higher resolution models when interpolated. In order to have files that have values at all ocean points from the coarse ESM-grid, the annual files are extrapolated towards the land boundary before generating the boundary conditions. 

The following command is a generic one that does the extrapolation and creates the nesting archive files in the output nest folder in the target directory. It assumes that you have performed the bias correction as described later in the document.

Unlike the example above for physics files, the command can take an extra option to activate biogeochemistry post-processing. This can be performed by providing the path to the biogeochemical file directory. The following example assumes the are also located in $Nesting_Files_PATH provided in REGION.src.

```
../bin/Nesting_noresm/Generate_nesting_files_year.sh <year> <region, e.g. TP2a0.10> <experiment, e.g. 010> $Nesting_Files_PATH
```

- The command above already includes the following commands. It is provided here as an example. No need to do it again. 

```
../bin/Nesting_noresm/separate_and_extrapolate_files_year.sh  $varible $year
../bin/Nesting_noresm/noresm_to_hycom.sh ../../TZ4a0.10/expt_01.1/ ../../NORESM_Nesting/thetao_Omon_NorESM2-MM_historical_r1i1p1f1_gr_195512_extrap.nc -b <path_to_bgc_files>
```

## Generating files with model bias used for bias correction
- Find the representaive period of the climatology you are using, e.g. 2000 - 2009.
- Use the script Create_ESM_climatology.sh, set "syear" start year, "eyear" end year, and "cstr" appropriate to the file names and run the script.
- Generate a climatology from the earth system model for the represetative period. Change the ESM and scenario for the following example command. Change the start and end years for your desired output climatology. Note that the ESM decadal averages will likely have a name of this format *_gr_200001-200912.nc, you may still provide different output years for the following script.
- Assign the folder that contains the nesting files in the REGION.src file, e.g. export Nesting_Files_PATH=/cluster/work/users/$USER/$R/Nesting_files/ 
```
../bin/Nesting_noresm/Create_ESM_climatology.sh  2000 2009 _Omon_NorESM2-MM_historical_r1i1p1f1_gr_
```
- Regrid the observational climatolgy to the ESM grid, below a certain depth, seasonal rather than monthly values must be used.
- This can be done uaing the script `Create_climatology_for_bias_correction.sh', using the keywords: *temperature, salinity, oxygen, nitrate, phosphate, silicate* separately.
```
../bin/Nesting_noresm/Create_climatology_for_bias_correction.sh silicate woa=/path/to/WOA
```
	

