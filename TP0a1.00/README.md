This is the REGION directory. Here is the typical content of this directory:

    ├── expt_01.0     # Experiment directory
    ├── force         # forcing input data for experiments
    ├── README.md      
    ├── REGION.src    # source file containing important environment variables and paths
    ├── relax         # Relaxation input data for experiments
    └── topo          # Topography and domain input data for experiments

The directory will contain input and output files for different experiments. All experiments inside this directory will have the same model domain (described by the regional.grid files in the topo subdirectory), but may have different bathymetries.



# Scripts for setting up new experiments and regions

The different experiments can be found in subdirectories expt_01.0, expt_01.1, etc etc. It is possible to create new experiments based on existing ones, by running the script ../bin/expt_new.sh .  Running the script without arguments will give you a brief explanation of its use. As an example, the following can be used to create a new experiment 01.1 from experiment 01.0

    ../bin/expt_new.sh 01.0 01.1

After a new experiment has been set up you should modify the configuration files to suit your new setup.

You can also set up a new region using the script ../bin/newRegion.sh. The following will create a new region dir in ../TP3a0.12/

   ../bin/newRegion.sh TP4a0.12 ../

# Region names
Region naming is basically up to you. The only restriction is that the name contains no underscore (_) and no semicolons (;). I have chosen to use a naming based on a three letter identifier, a version letter and a indication of resolution in degrees. The example above uses the region name 
    
    TP4a0.12

which consists of "TP4" (for TOPAZ 4), "a" as a version number (in case you want to do modifications to domain etc), and 0.12, which is an indication of a 1/8 th degree grid (roughly).


# REGION.src

Most paths to datasets etc can be set in REGION.src. It also contains the relative paths to hycom support routines and python routines. It is possible for you to modify these if you wish, by setting variables in REGION.src. By default the path is set to relative paths in the NERSC-HYCOM-CICE distro. 

If you want to use a different path, where routines etc are stored in your home dir, for instance, you must set these variables to point in the right direction in REGION.src :

    export HYCOM_ALL=$mydir/../hycom/hycom_ALL/hycom_2.2.72_ALL/
    export HYCOM_PYTHON_ROUTINES=$mydir/../bin/
    export INPUTDIR=$mydir/../input/
    export MSCPROGS=$mydir/../hycom/MSCPROGS/

NB: Setting this to a different location hasnt been tested much, but should be theoretically possible  ...




