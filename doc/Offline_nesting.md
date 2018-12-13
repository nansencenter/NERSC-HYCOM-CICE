
#  Description of HYCOM-CICE offline nesting usage

Details about the required procedures to carry out offline nesting can be found in the HYCOM USER GUIDE (Wallcraft, 2003, available at http://hycom.rsmas.miami.edu/hycom-model/ documentation.html). Here, we closely follow the standard HYCOM nesting approach based on the MERCATOR GLOBAL_ANALYSIS_FORECAST_PHY_001_024 product (as outer model data).

Nesting in current version need to be started from source regional directory/experiment (i.e. NEMO folder here NMOa0.08/expt_01.0 for native grids and NMOa0.09 for regular grids); you need to specify your target experiment directory (i.e. TOPAZ experiment directory TP5a0.06/expt_03.0). REGION.src need to possess the paths of NEMO mesh/coordinate and data netcdf files. After runing the following line, for example for native grid, 

../bin/nemo_to_hycom.sh ../../TP5a0.06/expt_03.0/ /nird/projects/nird/NS9481K/MERCATOR_DATA/PHY/2012/ext-GLORYS12V1_1dAV_20120312_20120313_grid2D_R20120314.nc 

You should expect to have your horizontally/vertically interpolated files in the TOPAZ nesting experiment folder, in this case TP5a0.06/030/archv.XXX_XXX.[ab] because you specified target experiment expt_03.0. Please note I link the "bin" directory one back the experiment folder in all region directories for example, if my current directory is TP5a0.06 it would be ln -sf ~/NERSC-CICE/bin .



#  Quick start

Based on mesh and coordinate MERCATOR netcdf files, the archive data files for grid and bathymetry are reconstructed on the same horizontal grid cells of NEMO model output. Having these files, we need to do following basic steps to create nesting archive files:

(1) Do tiling when running on shared/distributed memory machine (/bin/tile_grid.sh).

(2) Produce mapping index [ab] files (using /bin/isuba_gmap..sh).

(3) Interpolate horizontally the outer model fields (NEMO model outputs) to the inner subdomain (TOPAZ5).

(4) Choose vertical isopycnal structure (located in blkdat.input) and do vertical interpolation of temp, salin, u-vel., v-vel., and thknness according to the chosen vertical structure.

(5) Build port files.

#  Directory structure

Following is the general structure of HYCOM-CICE directory:

	├── bin                 # Location of binaries and python routines
	├── cice                # Location of CICE code
	├── doc                 # Documentation in markdown format
	├── hycom               # Location of hycom code and utilities
	│   ├── hycom_ALL       # Location of setup/diag routines 
	│   ├── MSCPROGS        # Location of setup/diag routines  developed at NERSC
	│   └── RELO            # Location of hycom source code
	├── input               # Location of some input files 
	├──  TP0a1.00           # Location of "Reference experiment"
	├──  TP5a0.06           # Location of "Reference experiment for TOPAZ5 system"
	└──  NMOa0.08           # Location of "Reference experiment for nesting from NEMO" 				# for TOPAZ5 system"


We aim to produce nesting archive [ab] files for TOPAZ5 region (here TP5a0.06 as inner model region) from the NEMO global data as outer model (i.e. NMOa0.08 directory). 
Therefore your working directory (e.g. in Sisu: /wrk/pr2nXXXX/NERSC-HYCOM-CICE) should contain two folders: Target (inner model) directory, TPa0.06; and Source (outer model) directory, NMOa0.08. All processing programs to produce the target archive data files are issued from the (source) region directory, i.e. NMOa.08. 

The following illustrates how these two directories are organised in the presence of nesting. 

## (1) before starting nesting

    └── NMOa0.08             # Region directory for NEMO files
        └── bin             # Link to bin utility folder (this directory should be in the region directory)
        └── expt_01.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file
Topo subdirectory in NMOa0.08 includes NEMO grid and bathymetry data files in [ab] format.

    └── TP5a0.06             # Region directory for TOPAZ5 region
        └── bin              # Link to bin utility folder including all required python codes and scripts (this directory should be in region directory)
        └── expt_01.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file

## (2) after applying nesting

    └── NMOa0.08             # Region directory for NEMO files
        └── bin.             # Link to bin utility folder (this directory should be in region directory)
        └── expt_01.0        # Experiment directory
            ├── data         # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file
        └── subregion        # subregion folder
            ├── 010          # a directory named with the same experiment number
Data folder includes dummy (intermediate) archive files converted from netcdf to [ab] archive format. Subregion folder contains horizontally interpolated archive files from dummy [ab] files (in outer domain grid) to the new grid/sub-domain (here TP5a0.06 grid).

    └── TP5a0.06             # Region directory for TOPAZ5 region
        └── bin.             # Link to bin utility folder (this directory should be in region directory)
        └── expt_03.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # region configuration file
        └── nest             # actual directory containing all nesting archive files
            ├── 030          # a directory named with the same experiment number, i.e. 030

Nest folder includes all vertically interpolated archive data files from the horizontally interpolated files in outer region directory (i.e. NMOa0.08/subregion/010/).

# The NEMO grid

The daily product of Operational Mercator global ocean analysis and forecast system with NEMO ocean model (at 1/12 degree) starts on December 27, 2006. The model is used here as outer model including daily fields of temperature, salinity, currents, sea level, mixed layer depth and ice parameters. The global ocean output files are represented on 1/12 degree horizontal resolution (no staggered grid) with regular longitude/latitude equirectangular projection. The output files cover 50 vertical levels ranging from 0 to 5500 meters.

The NEMO horizontal mesh on T-cell are converted to [ab] archive files using the coordinate/mesh netcdf files provided by MERCATOR. All pre- and post-processing programs read these [ab] files which are produced only once. 

# Building NEMO archive files

The data from MERCATOR GLOBAL_ANALYSIS_FORECAST_PHY_001_024 are input for a python code "nemo2archvz.py" which converts temp, salon, u-vel., v-vel., and thinness fields from their non-native grid forms (e.g. regular grid) to the same grid but in HYCOM convention, I.e. [ab] archive data files. If the program run successfully a directory called "data" are created and the archive files, archv.YYYY_ddd_hh.[ab], are copied there.

      └── NMOa0.08            
         └── REGION.src       
         └── bin             
         └── expt_01.0        
             └── data             
                   └── archv.YYYY_ddd_hh.[ab]             
         └── topo             

* HYCOM may like "hh" appeared in the name of archive files to be at "00" and MERCATOR GLOBAL_ANALYSIS_FORECAST_PHY_001_024 daily data are archived at "12". Therefore, the python code do also a time average to cope with this issue (in the case data for the next day are available).

# Horizontal interpolation to the inner domain

Our current directory, hereafter, is NMOa0.08/expt_01.0:

      └── NMOa0.08            
         └── bin             
         └── expt_01.0        
         └── topo             
         └── REGION.src       

After establishing grid and bathymetry files in "topo" folder, the user must to start interpolation from the outer domain to the inner domain (i.e. TP5a0.06).In order to do this interpolation, an index mapping matrix (from outer region to inner region) is required which is used by $HYCOMALL/subregion/src/isubaregion.f. The corresponding bash script in the "bin" folder is used as follows:
       
         ../bin/isuba_gmapi.sh. $target_region

where target_region denotes ../../TP5a0.06 according to previously shown structure. If remapping [ab] files can be created successfully, they are located in subregion folder as follows:


      └── NMOa0.08            
         └── REGION.src       
         └── bin             
         └── expt_01.0        
         └── topo             
         └── subregion             
             └── 010             
             └── TP5a0.06.gmap.[ab]             

In script "../bin/nemo2hycom.sh",program "../bin/remap_nemo.sh" is executed by

    ../bin/remap_nemo.sh $target_experiment $subregion_archz

Where "target_experiment" is "TPa0.05/expt_03.0" and "subregion_archz" is set to "data/archv.YYYY.ddd.[ab]" generated at the previous step using python script.

use executable "$HYCOMALL/subregion/src/isubaregion" using following inputs/outputs:

   
    --- 'flnm_reg'  = target sub-region grid       filename
    --- 'flnm_map'  = target sub-region grid map   filename
    --- 'flnm_top'  = target bathymetry filename, or 'NONE'
    --- 'flnm_tin'  = input  bathymetry filename, or 'NONE'
    --- 'flnm_in'   = input  archive    filename
    --- 'flnm_out'  = output archive    filename
    --- 'cline_out' = output title line (replaces preambl(5))
         

Now, the interpolated archive files are created and located in the "subregion/010" director: 

         └── subregion             
             └── 010             
                 └── archv.YYYY_ddd_hh_L.[ab]             
             └── TP5a0.06.gmap.[ab]             


# Vertical interpolation using the chosen vertical (isopycnal) structure

Parameters required for the vertical structure are acquired from the "blkdat.input" (located in the target experiment folder. For example, sigma, dp00s,…). We have our dummy (intermediate) archive files including temp, salon, u-vel., v-vel., thkness
 Fields at "subregion/010" which are still in z-levels. Using an executable "$HYCOMALL/relax/src/nemo_archvz", all mentioned fields are converted from z-levels to isopycnal, see "../bin/archvz2hycom.sh". By this the target archive data files are created and copied in "TP5a0.0/nest/030/archv.YYYY_ddd_hh.[ab]". 
         

    └── TP5a0.06            
        └── bin.    
        └── expt_03.0        
        └── topo            
        └── REGION.src      
        └── nest             
             └── 030             
                 └── archv.YYYY_ddd_hh.[ab]             


# Generating port data

To create ports and relaxation zones as part of nesting procedure, following script is used (we are now in TP5a0.06/expt_030 directory):

   ../bin/nest_setup_ports.sh ${width_of_relax_zone} ${efold_time_in_day}

Where "width_of_relax_zone" and "efold_time_in_day" can be set, for example, to 20 and 20, respectively. If program runs successfully, a file called "ports.nest" and "rmu.[ab]" are created and stored in "../nest/030":
   
    └── TP5a0.06            
        └── bin.    
        └── expt_03.0        
        └── topo            
        └── REGION.src      
        └── nest             
             └── 030             
                 └── archv.YYYY_ddd_hh.[ab]             
                 └── ports.nest            
                 └── rmu.[ab]            
