This is the REGION directory. Here is the typical content of this directory:

    ├── expt_01.0     # Experiment directory
    ├── force         # forcing input data for experiments
    ├── README.md      
    ├── REGION.src    # source file containing important environment variables and paths
    ├── relax         # Relaxation input data for experiments
    └── topo          # Topography and domain input data for experiments

The directory will contain input and output files for different experiments. All experiments inside this directory will have the same model domain (described by the regional.grid files in the topo subdirectory), but may have different bathymetries.

The different experiments can be found in subdirectories expt_01.0, expt_01.1, etc etc. It is possible to create new experiments based on existing ones, by running the script ../bin/expt_new.sh .  Running the script without arguments will give you a brief explanation of its use. As an example, the following can be used to create a new experiment 01.1 from experiment 01.0

    ../bin/expt_new.sh 01.0 01.1

After a new experiment has been set up, you should modify the configuration files to suit your new setup.
