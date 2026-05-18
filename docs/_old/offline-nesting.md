
#  Description of HYCOM-CICE offline nesting usage

Details about the required procedures to carry out offline nesting can be found in the HYCOM USER GUIDE (Wallcraft et. al., 2003, available at https://www.hycom.org/hycom/documentation). Here, we closely follow the standard HYCOM nesting approach based on the MERCATOR GLOBAL_ANALYSIS_FORECAST_PHY_001_024 product for physics and GLOBAL_ANALYSIS_FORECAST_BIO_001_029 product for biology (as outer model data).

Nesting in current version need to be started from source regional directory/experiment (i.e. NEMO folder here NMOa0.08/expt_01.1 for native grids and NMOb0.08/expt_01.0 for regular grids); you need to specify your target experiment directory (i.e. TOPAZ experiment directory, e.g. TP5a0.06/expt_03.0, we will use this target experiment as example in this document hereafter). REGION.src need to possess the paths of NEMO mesh/coordinate and data netcdf files. After running the following line, for example for native grid, 

../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_03.0/ /nird/projects/NS9481K/MERCATOR_DATA/PHY/2012/ext-GLORYS12V1_1dAV_20120312_20120313_grid2D_R20120314.nc 

You should expect to have your horizontally/vertically interpolated files in the TOPAZ nesting experiment folder, in this case TP5a0.06/nest/030/archv.XXX_XXX.[ab] because you specified target experiment expt_03.0. Please link the "bin" directory one back the experiment folder in all region directories, for example, if your current directory is TP5a0.06, it would be `ln -sf ~/NERSC-HYCOM-CICE/bin .`.


It is worth noting that you have to run the first year from climatology and then introduce the nesting from a restart file.


#  Quick start

Based on mesh and coordinate MERCATOR netcdf files, the archive data files for grid and bathymetry are reconstructed on the same horizontal grid cells of NEMO model output. Having these files, we need to do following basic steps to create nesting archive files:

(1) Produce mapping index [ab] files (using /bin/isuba_gmapi.sh).

(2) Interpolate horizontally the outer model fields (NEMO model outputs) to the inner subdomain (TOPAZ) and do vertical interpolation according to the chosen vertical structure located in blkdat.input (using /bin/nemo_to_hycom.sh).

(3) Modify montgomery potential values (using /bin/calc_montg1.py).

(4) Build port files (using /bin/nest_setup_ports.sh).


#  How to run on native MERCATOR grid:
Please change path to NMOa0.08/expt_01.1 (for more information look at REGION.src for path of an example, i.e., NEMO_EXAMPLE_NATIVE_FOLDER=/nird/projects/NS9481K/MERCATOR_DATA/NMOa0.08). 
After applying `../bin/isuba_gmapi.sh ../../TP5a0.06/` once, run following line:

(1)  Without biology:

../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_03.0/ /nird/projects/NS9481K/MERCATOR_DATA/PHY/2013/ext-GLORYS12V1_1dAV_20131122_20131123_grid2D_R20131127.nc 

(2)  With biology (note that biology files are all on regular grid):

../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_03.0/ /nird/projects/NS9481K/MERCATOR_DATA/PHY/2013/ext-GLORYS12V1_1dAV_20131122_20131123_grid2D_R20131127.nc -b /nird/projects/NS9481K/MERCATOR_DATA/BIO/DAILY/2013/global_analysis_forecast_bio_20131122.nc

#  How to run on regular MERCATOR grid:
Please change path to NMOb0.08/expt_01.0 (for more information look at REGION.src for path of an example, i.e., NEMO_EXAMPLE_REGULAR_FOLDER=/nird/projects/NS9481K/MERCATOR_DATA/NMOb0.08). 
After applying `../bin/isuba_gmapi.sh ../../TP5a0.06/` once, run following line:

(1) Without biology:

../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_03.0 /nird/projects/NS9481K/MERCATOR_DATA/PHY/2018/MERCATOR-PHY-24-2018-01-01-12.nc -m regular

(2) With biology:

../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_03.0 /nird/projects/NS9481K/MERCATOR_DATA/PHY/2018/MERCATOR-PHY-24-2018-01-01-12.nc -m regular -b /nird/projects/NS9481K/MERCATOR_DATA/BIO/DAILY/2018/global_analysis_forecast_bio_20180101.nc

> [!NOTE]
> The following text explains how "nemo_to_hycom.sh" works in detail. Everything between the horizontal lines below has been taken care of by "nemo_to_hycom.sh". You do **not** have to do the following steps by hand. To skip the details, you can jump directly to "Calculating montg values".

---

#  Directory structure

Following is the general structure of HYCOM-CICE directory:

	├── bin                 # Location of binaries and python routines
	├── cice                # Location of CICE code
	├── doc                 # Documentation in markdown format
	├── hycom               # Location of hycom code and utilities
	│   ├── hycom_ALL       # Location of setup/diag routines 
	│   ├── MSCPROGS        # Location of setup/diag routines developed at NERSC
	│   └── RELO            # Location of hycom source code
	├── input               # Location of some input files 
	├──  TP0a1.00           # Location of "Reference experiment"
	├──  TP2a0.10           # Location of "Reference experiment for TOPAZ2 system"
	├──  TP5a0.06           # Location of "Reference experiment for TOPAZ5 system"
	├──  NMOa0.08           # Location of "Reference experiment for nesting from NEMO of native grid"
	└──  NMOb0.08           # Location of "Reference experiment for nesting from NEMO of regular grid"


We aim to produce nesting archive [ab] files for TOPAZ region (here TP5a0.06 as an example inner model region) from the NEMO global data as outer model (i.e. NMOa0.08 directory). 
Therefore your working directory should contain two folders: Target (inner model) directory, e.g. TP5a0.06; and Source (outer model) directory, e.g. NMOa0.08. All processing programs to produce the target archive data files are issued from the (source) region directory, i.e. NMOa/b0.08. 

The following illustrates how these two directories are organised in the presence of nesting. 

## (1) before starting nesting

    └── NMOa0.08             # Region directory for NEMO files
        └── bin              # Link to bin utility folder (this directory should be in the region directory)
        └── expt_01.1        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file
Topo subdirectory in NMOa0.08 includes NEMO grid and bathymetry data files in [ab] format.

    └── TP5a0.06             # Region directory for TOPAZ5 region
        └── bin              # Link to bin utility folder including all required python codes and scripts (this directory should be in region directory)
        └── expt_03.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file

## (2) after applying nesting

    └── NMOa0.08             # Region directory for NEMO files
        └── bin              # Link to bin utility folder (this directory should be in region directory)
        └── expt_01.1        # Experiment directory
            ├── data         # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file
        └── subregion        # subregion folder
            ├── 030          # a directory named with the same experiment number as in target region
Data folder includes dummy (intermediate) archive files converted from netcdf to [ab] archive format. Subregion folder contains horizontally interpolated archive files from dummy [ab] files (in outer domain grid) to the new grid/sub-domain (here TP5a0.06 grid).

    └── TP5a0.06             # Region directory for TOPAZ5 region
        └── bin              # Link to bin utility folder (this directory should be in region directory)
        └── expt_03.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file
        └── nest             # actual directory containing all nesting archive files
            ├── 030          # a directory named with the same experiment number, i.e. 030

Nest folder includes all vertically interpolated archive data files from the horizontally interpolated files in outer region directory (i.e. NMOa0.08/subregion/030/).

# The NEMO grid

The daily product of Operational Mercator global ocean analysis and forecast system with NEMO ocean model (at 1/12 degree) is used here as outer model including daily fields of e.g. temperature, salinity, currents, sea level, mixed layer depth and ice parameters. The global ocean output files are represented on 1/12 degree horizontal resolution (no staggered grid) with regular longitude/latitude equirectangular projection. The output files cover 50 vertical levels ranging from 0 to 5500 meters.

The NEMO horizontal mesh on T-cell are converted to [ab] archive files using the coordinate/mesh netcdf files provided by MERCATOR. All pre- and post-processing programs read these [ab] files which are produced only once. 

# Building NEMO archive files

The data from MERCATOR GLOBAL_ANALYSIS_FORECAST_PHY_001_024 are input for a python code "nemo2archvz_native.py" or "nemo2archvz_regular.py" which converts temp, salin, u-vel., v-vel., and thknss fields from their native or regular grid to the same grid but in HYCOM convention, i.e. [ab] archive data files. If the program run successfully a directory called "data" are created and the archive files, archv.YYYY_ddd_hh.[ab], are copied there.

      └── NMOa0.08            
         └── REGION.src       
         └── bin             
         └── expt_01.1        
             └── data             
                   └── archv.YYYY_ddd_hh.[ab]             
         └── topo             

* HYCOM may like "hh" appeared in the name of archive files to be at "00" and MERCATOR GLOBAL_ANALYSIS_FORECAST_PHY_001_024 daily data are archived at "12". Therefore, the python code do also a time average to cope with this issue (in the case data for the next day are available).

# Horizontal interpolation to the inner domain

Our current directory, hereafter, is NMOa0.08/expt_01.1:

      └── NMOa0.08            
         └── bin             
         └── expt_01.1        
         └── topo             
         └── REGION.src       

After establishing grid and bathymetry files in "topo" folder, the user must start interpolation from the outer domain to the inner domain (e.g. TP5a0.06). In order to do this interpolation, an index mapping matrix (from outer region to inner region) is required which is used by $HYCOMALL/subregion/src/isubaregion.f. The corresponding bash script in the "bin" folder is used as follows:
       
         ../bin/isuba_gmapi.sh. $target_region

where target_region denotes ../../TP5a0.06/ according to previously shown structure. If remapping [ab] files can be created successfully, they are located in subregion folder as follows:


      └── NMOa0.08            
         └── REGION.src       
         └── bin             
         └── expt_01.1        
         └── topo             
         └── subregion             
             └── 030             
             └── TP5a0.06.gmap.[ab]             

In script "../bin/nemo_to_hycom.sh",program "../bin/archvz2hycom_biophys.sh" is executed by

    ../bin/archvz2hycom_biophys.sh $target_experiment $subregion_archz

Where "target_experiment" is "TP5a0.06/expt_03.0" and "subregion_archz" is set to "data/archv.YYYY.ddd.[ab]" generated at the previous step using python script.

use executable "$HYCOMALL/subregion/src/isubaregion" using following inputs/outputs:

   
    --- 'flnm_reg'  = target sub-region grid       filename
    --- 'flnm_map'  = target sub-region grid map   filename
    --- 'flnm_top'  = target bathymetry filename, or 'NONE'
    --- 'flnm_tin'  = input  bathymetry filename, or 'NONE'
    --- 'flnm_in'   = input  archive    filename
    --- 'flnm_out'  = output archive    filename
    --- 'cline_out' = output title line (replaces preambl(5))
         

Now, the interpolated archive files are created and located in the "subregion/030" directory: 

         └── subregion             
             └── 030             
                 └── archv.YYYY_ddd_hh_L.[ab]             
             └── TP5a0.06.gmap.[ab]             


# Vertical interpolation using the chosen vertical (isopycnal) structure

Parameters required for the vertical structure are acquired from the "blkdat.input" (located in the target experiment folder. For example, sigma, dp00s,…). We have our dummy (intermediate) archive files including temp, salin, u-vel., v-vel., thknss
 Fields at "subregion/030" which are still in z-levels. Using an executable "$HYCOMALL/relax/src/nemo_archvz_biophys", all mentioned fields are converted from z-levels to isopycnal, see "../bin/archvz2hycom_biophys.sh". By then the target archive data files are created and copied in "TP5a0.06/nest/030/archv.YYYY_ddd_hh.[ab]". 
         

    └── TP5a0.06            
        └── bin    
        └── expt_03.0        
        └── topo            
        └── REGION.src      
        └── nest             
             └── 030             
                 └── archv.YYYY_ddd_hh.[ab]             

---

# Calculating montg values 

To replace the zero values of Montgomery potential variable in the nesting archive files, which were already interpolated both horizontally and vertically, run the following code in your target experiment directory:
   python ../bin/calc_montg1.py archv_file restart_file output_path

Note that you need "any" of the restart file from your target experiment and you need to define the output_path. 
Switch to your target experiment directory first (e.g. TP5a0.06/expt_03.0), then run
   
   python ../bin/calc_montg1.py ../nest/030/archv.2013_326_00.a ./data/restart.2011_005_00_0000.a  ../nest/030/Montg/


# Generating port data

To create ports and relaxation zones as part of nesting procedure, following script is used (we are now in TP5a0.06/expt_03.0 directory):

   ../bin/nest_setup_ports.sh ${width_of_relax_zone} ${efold_time_in_day}

Where "width_of_relax_zone" and "efold_time_in_day" can be set, for example, to 20 and 20, respectively. If program runs successfully, a file called "ports.nest" and "rmu.[ab]" are created and stored in "../nest/030":
   
    └── TP5a0.06            
        └── bin    
        └── expt_03.0        
        └── topo            
        └── REGION.src      
        └── nest             
             └── 030             
                 └── archv.YYYY_ddd_hh.[ab]             
                 └── ports.nest            
                 └── rmu.[ab]   


# Using BIO-NESTING

Bio-nesting follows the same procedure as physics variables. We have for now three vaiables included, no3, pho and sil. Other variables, e.g. oxygen can be added later. 



