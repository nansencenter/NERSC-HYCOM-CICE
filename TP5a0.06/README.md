[toc]
This is the region directory. Here is the typical content of this directory:

    ├── expt_01.0     # Experiment directory (for standalone HYCOM model)
    ├── expt_02.0     # Experiment directory (for coupled HYCOM-CICE model)
    └── topo          # Topography and domain input data for experiments

The directory will contain input and output files for different experiments. All experiments inside this directory will have the same model domain (described by the regional.grid files in the topo subdirectory), but may have different bathymetries.

