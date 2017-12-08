
#  Description of HYCOM-CICE offline nesting usage

Details about the required procedures to carry out offline nesting can be found in the HYCOM USER GUIDE (Wallcraft, 2003, available at http://hycom.rsmas.miami.edu/hycom-model/ documentation.html). Here, we closely follow the the standard HYCOM nesting approach for the MERCATOR GLO-PHY-024 product.


#  Directory structure

Following is the general structure of HYCOM-CICE directory hierarchy:

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
└──  NMOa0.08           # Location of "Reference experiment for nesting from NEMO for TOPAZ5 system"


We aim to produce nesting archive [ab] files for HYCOM-CICE region as inner model (here TP5a0.06) from the NEMO global data as outer model (i.e. NMOa0.08). 
Therefore your working directory should contain two folder: Target (inner model) region folder TPa0.06; and Source (outer model) region folder NMOa0.08. All processing programs to produce the target archive files are run in the (source) region folder of NMOa.08. 

The following illustrates how these two directories are usually organised 

**(1) before starting nesting procedure

    └── TP5a0.06             # Region directory for TOPAZ5 region
        └── expt_01.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # creation configuration file

    └── NMOa0.08             # Region directory for NEMO files
        └── expt_01.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # creation configuration file
**(2) after applying nesting

    └── TP5a0.06             # Region directory for TOPAZ5 region
        └── expt_03.0        # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # creation configuration file
        └── nest             # creation configuration file
            ├── 030          # a directory named with the same experiment number, i.e. 030

    └── NMOa0.08             # Region directory for NEMO files
        └── expt_01.0        # Experiment directory
            ├── data         # Experiment directory
        └── topo             # topography and grid directory
        └── REGION.src       # creation configuration file
        └── subregion        # creation configuration file
            ├── 010          # a directory named with the same experiment number

