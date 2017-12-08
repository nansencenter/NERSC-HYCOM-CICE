
#  Description of HYCOM-CICE offline nesting usage

Details about the required procedures to carry out offline nesting can be found in the HYCOM USER GUIDE (Wallcraft, 2003, available at http://hycom.rsmas.miami.edu/hycom-model/ documentation.html). Here, we closely follow the standard HYCOM nesting approach based on the MERCATOR GLOBAL_ANALYSIS_FORECAST_PHY_001_024 product (as outer model data).

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

## (1) before starting nesting procedure

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

# Horizontal interpolation to the inner domain

Our current directory, hereafter, is NMOa0.08/expt_01.0:

      └── NMOa0.08            
         └── bin             
         └── expt_01.0        
         └── topo             
         └── REGION.src       

After establishing grid and bathymetry files in "topo" folder, the user must to start interpolation from the outer domain to the inner domain (i.e. TP5a0.06).In order to do this interpolation, an index mapping matrix (from outer region to inner region) is required which is used by $HYCOMALL/subregion/src/isubaregion.f. The corresponding bash script in the "bin" folder is used as follows:
       bin/isuba_gmapi.sh. $target_region
where target_region denotes ../../TP5a0.06 according to previously shown structure.



