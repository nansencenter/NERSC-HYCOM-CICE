# Receipt of how to generate nesting conditions from an ESM:

## Set up the directory representing the outer model.
- Create a directory called ESMa1.00 and inside:
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
- Copy REGION.src to this directory and make sure R and NHCROOT are defined correctly.
```
cp $HOME/NERSC-HYCOM-CICE/input/REGION.src .
```
- Get the grid information from the model in questions, this can be found on the CMIP6 server: you need the file for the varible "areacello" and "deptho" and put in topo-folder.
- From the experiment folder run script:
```
../bin/Nesting_noresm/make_grid_noresm2hycom.py ../topo/areacello_file.nc
```
This will generate the regional.depth and regional.grid-files for the global mode grid.
- Place these in the topo-folder
- From the experiment folder run the following script to generate the mapping from the global model to the regional hycom.:
```
../bin/isuba_gmapi.sh. $target_region
```
do we need to copy blkdat to experimetn folder?
			
## Physical nesting conditions 
- Download the nesting-files needed: for physics: uo, vo, zos, thetao, so, make sure you get the files on a with varibles on z-levels (if the model is not a z-level model).
- In order to have files that have values at all ocean points from the coarse ESM-grid, the annual files are pre-processed before generating the boundary conditions:
```
../bin/Nesting_noresm/separate_and_extrapolate_files_year.sh  $varible $year
```
- This can also be called in the script:
```
../bin/Nesting_noresm/Generate_nesting_files_year.sh $year
```
- The bias corection happend in this step, see *Generating files with model bias used for bias correction* below for how to make files for bias correction.
- For NorESM, the files for bias correction can be found on Fram in /cluster/projects/nn9481k/NORESM_bias/
- Once the ESM files have been prepared the nesting condioting can be generated using the following command, which must be called from the experiment-directory of ESMa1.00:
```
../bin/Nesting_noresm/noresm_to_hycom.sh ../../TZ4a0.10/expt_01.1/ ../../NORESM_Nesting/thetao_Omon_NorESM2-MM_historical_r1i1p1f1_gr_195512_extrap.nc
```

		
## Biogeochamical nesting: Download the nesting-files needed: for physics: no3, po4, si, o2, these are decadal files.
- Download the nesting-files needed: for BGC: no3, po4, si, and o2, make sure you get the files on a with varibles on z-levels (if the model is not a z-level model).
- In order to have files that have values at all ocean points from the coarse ESM-grid, the annual files are pre-processed before generating the boundary conditions:
```
../bin/Nesting_noresm/separate_and_extrapolate_files_year.sh  $varible $year
```
- This can also be called in the script:
```
../bin/Nesting_noresm/Generate_nesting_files_year.sh $year
```
- The bias corection happend in this step, see *Generating files with model bias used for bias correction* below for how to make files for bias correction.
- For NorESM, the files for bias correction can be found on Fram in /cluster/projects/nn9481k/NORESM_bias/
- Once the ESM files have been prepared the nesting condioting can be generated using the following command, which must be called from the experiment-directory of ESMa1.00:
```
../bin/Nesting_noresm/noresm_to_hycom.sh ../../TZ4a0.10/expt_01.1/ ../../NORESM_Nesting/thetao_Omon_NorESM2-MM_historical_r1i1p1f1_gr_195512_extrap.nc -b <path_to_bgc_files>
```

## Generating files with model bias used for bias correction
- Find the representaive period of the climatology you are using.
- Generate a climatology from the earth system model for the represetative period.
- Regrid the observational climatolgy to the ESM grid, below a certain depth, seasonal rather than monthly values must be used.
- This can be done uaing the script `Create_climatology_for_bias_correction.sh'
- 

	

