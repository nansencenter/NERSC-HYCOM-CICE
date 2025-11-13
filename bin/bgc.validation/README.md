# BGC VALIDATION SCRIPTS

The purpose of this folder is to provide the set of validation scripts in one place, and be a guide for a generic model evaluation.
The codes evaluate the model output against:
1) in situ data
2) satellite data
3) BGC-Argo profiles

Before the scripts below, it is important to add the bgc.validation folder to your python path, e.g.:
```
export PYTHONPATH="$HOME/NERSC-HYCOM-CICE/bin/bgc.validation:$PYTHONPATH"
```

## Regional masks

The validation scripts aim to unify the process. As such, a common region definitions is essential. The commands below creates these masks as binary pickle files and stores them in the same folder. The files are already created and copied to cluster/projects (i.e. /cluster/projects/nn2993k/BGC.Validataion/). There is no need to repeat this process. The command below is here for documentation as an example:
```
python make_OM_regional_masks.py TP5a0.06 /cluster/work/users/cagyum/TP5a0.06/topo/regional.grid
```

## In situ data preprocess

The unified in situ data from various source (e.g. Copernicus, GLODAP) are stored as text files in "/cluster/projects/nn2993k/BGCDATA/prepobs_bgc/". The command below reads the text files, and creates python dictionaries that has information such as years, depth, regions. These dictionaries are saved as a binary file in the same folder. The file is already created and copied to NIRD (i.e. /cluster/projects/nn2993k/BGC.Validation/INSITU.OMmask.pckl). There is no need to repeat this process. The command below is here for documentation as an example:
```
python convert_txt_pckl_insitu_preobs_OMmask.py ~/NERSC-HYCOM-CICE/ 1997 2023
```

Once you have the INSITU.OMmask.pckl file, you can proceed with model to in situ data co-location process. You can run the following code in any directory. It will take mandatory and optional argumentns. You need to provide the region, experiment, start and end year, e.g. :
```
python $HOME/NERSC-HYCOM-CICE/bin/bgc.validation/model_vs_INSITU.py TP2 038 2006 2010
```
The code will by default assume your work directory is '/cluster/work/users/' (like in Betzy). If this is not the case, provide the optional argument --workdir=path-to-work. The code will retrieve your username automatically. If all goes well, it will look into the folder '/cluster/work/users/$USER/TP2a0.10/expt_03.8/data/' in this particular example.

The code will produce a series of output python pickle binary files in model data folder depending on the in situ data availability with names colocated...pckl. They will be used later.
