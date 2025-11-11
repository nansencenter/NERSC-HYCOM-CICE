# BGC VALIDATION SCRIPTS

The purpose of this folder is to provide the set of validation scripts in one place, and be a guide for a generic model evaluation.
The codes evaluate the model output against:
1) in situ data
2) satellite data
3) BGC-Argo profiles

## Regional masks

The validation scripts aim to unify the process. As such a common region definitions is essential. The commands below creates these masks as binary pickle files and stores them in NIRD (i.e. /nird/projects/NS9481K/BGC.Validataion/). The files are already created. There is no need to repeat this process. The command below is here for documentation as an example:
```
python make_OM_regional_masks.py TP5a0.06 /cluster/work/users/cagyum/TP5a0.06/topo/regional.grid
```