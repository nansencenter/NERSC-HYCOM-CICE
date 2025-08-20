[toc]
This is the region directory. Here is the typical content of this directory:

    ├── expt_01.0     # Experiment directory (for standalone HYCOM model)
    ├── expt_01.1     # Experiment directory (for coupled HYCOM-CICE model)
    ├── expt_02.3     # Experiment directory (implementation of HYCOM2.3 in HYCOM-CICE model)
    └── topo          # Topography and domain input data for experiments

The directory will contain input and output files for different experiments. All experiments inside this directory will have the same model domain (described by the regional.grid files in the topo subdirectory), but may have different bathymetries.

Note:

Under the expt_02.3 directory, there are more sample files used for two versions of the model HYCOM2.2 and HYCOM2.3. Clearly, these files could be different beacause the model's needs are not the same. 

Nesting frequence: From the blkdat.input, the nesting setting shows a difference. 'bnstfq' and 'nestfq' use '-1' in the blkdat.input_V2.3, which are different from that in blkdat.input_V2.2. 
It means the model needs the nesting files named by archm.yyyy_ddd_12.a[b] for '-1', otherwise, it would be archv.yyyy_ddd_00.a[b]. For the first implementation, the nesting boundary files are derived from the daily average. The second implementation of 'bnstfq=1' provides a clearer snapshot. For the TP5 implementation, the nesting files are interpolated by the daily source, so the recommendation is to use the setting of '-1' within the nesting boundary.

