[toc]
This is the REGION directory. Here is the typical content of this directory:

    ├── expt_01.0     # Experiment directory
    ├── force         # forcing input data for experiments
    ├── mysource      # Direcotry containing mod_hycom.F and hycom.F for running HYCOM standalone
    ├── README.md      
    ├── REGION.src    # source file containing important environment variables and paths
    ├── relax         # Relaxation input data for experiments
    └── topo          # Topography and domain input data for experiments

The directory will contain input and output files for different experiments. All experiments inside this directory will have the same model domain (described by the regional.grid files in the topo subdirectory), but may have different bathymetries.

# REGION.src

First of all, make a copy of the REGION.src file from the input folder into your working directory (TP0a1.00 for example).
Then, you need to set your PATH variable to point to $NHCROOT/bin. That way you don't have to write the full path to the commands you want to execute.

With that out of the way, here is 

Most paths to datasets etc can be set in REGION.src. It also contains the
relative paths to hycom support routines MSCprogs, shell scripts, etc. It is
possible for you to modify these if you wish, by setting variables in
REGION.src. By default the path is set to relative paths in the NERSC-HYCOM-CICE
distro when checked out of the repository.   

If you want to source code and routines in a separate location (for instance
somewhere on your your home directory), you will need to modify things in
REGION.src. 

Lets assume you have checked out the NERSC-HYCOM-CICE repository somewhere on
your home dir, but want to run the model on a work directory. First of all, move
your region dir to cesome place in the work directory. Then, edit REGION.src, and set NHCROOT to
the location of NHCROOT in your home directory. Below is an example :

   # Point to the location of NHCROOT/
   export NHCROOT=/home/nersc/knutali/hycom/nersc-hycom-cice/     # 
   export HYCOM_ALL=$NHCROOT/hycom/hycom_ALL/hycom_2.2.72_ALL/    # follows NHCROOT, but can be set differently
   export INPUTDIR=$NHCROOT/input/                                # follows NHCROOT, but can be set differently
   export BINDIR=$NHCROOT/bin/                                    # follows NHCROOT, but can be set differently
   export HYCOM_PYTHON_ROUTINES=$BINDIR                           # follows NHCROOT, but can be set differently
   export MSCPROGS=$NHCROOT/hycom/MSCPROGS/

This way your model will use the routines and source code in your home dir at
NHCROOT, but it will run on your work dir.


# Scripts for setting up new experiments and regions

The different experiments can be found in subdirectories expt_01.0, expt_01.1,
etc etc. It is possible to create new experiments based on existing ones, by
running the script expt_new.sh (found in $NHCROOT/bin).  Running the script without arguments
will give you a brief explanation of its use. As an example, the following can
be used to create a new experiment 01.1 from experiment 01.0 

    expt_new.sh 01.0 01.1

After a new experiment has been set up you should modify the configuration files to suit your new setup.

You can also set up a new region using the script region_new.sh in $NHCROOT/bin . The following will create a new region dir in ../TP3a0.12/

   region_new.sh TP4a0.12 ../

Note that this just sets up the directories in a consistent manner - you will
have to create grid files, topography files, etc on your own.

# Region names

Region naming is basically up to you. The only restriction is that the name contains no underscore (_) and no semicolons (;). I have chosen to use a naming based on a three letter identifier, a version letter and a indication of resolution in degrees. The example above uses the region name 
    
    TP4a0.12

which consists of "TP4" (for TOPAZ 4), "a" as a version number (in case you want to do modifications to domain etc), and 0.12, which is an indication of a 1/8 th degree grid (roughly).


# Scripts for setting up new region grid and bathymetry files

# Compile HYCOM standalone

For a detailed explanation on how to compile HYCOM standalone see $NHCROOT/doc/HYCOM-only-compilation.md

Quick start: To compile HYCOM as standalone (i.e. not coupled with CICE) you need to do the following:
1. In blkdat.input set ‘iceflg’ to 0.
2. Link $NHCROOT/TP0a1.00/expt_01.0/mysource to experiment directory
3. compile_model.sh -m fram -u ifort
