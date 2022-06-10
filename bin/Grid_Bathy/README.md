This directory contains routines for generating/processing new grids and bathymetry files.
` `

` $  -For new grid you need grid.info as input file for "hycom_grid.py" `

` $  -Then create the bathy file using "hycom_bathy.py" `

` $  -You need to topfix to manullay handle isolated wet points, etc`

` $  -Then create land mask using cice_kmt.sh`





# Generation of grids


|executable     | purpose|
|-------- | -------------|
|hycom_grid.py       | Create bathymetry based on conformal mapping or on proj4 projection. Can do  the same as old confmap routines in MSCPROGS. Ex.: hycom_grid.py confmap --filename ./grid.info|



# Generation of bathymetries

|executable     | purpose|
|-------- | -------------|
|hycom_bathy.py             | Generate hycom bathymetry for a predefined hycom grid. Ex: hycom_bathy.py --cutoff 3 regional.grid.a |
|hycom_bathy_merge.py       | Merges two bathymetries on the same grid. |
|hycom_bathy_consistency.py | Does a consistency check of a hycom bathymetri (removes isolated basins, single cell islands, etc). Ex: hycom_bathy_consistency.py depth_TP5a0.06_01.a|
|cice_kmt.py | Simple script for creating CICE land mask from hycom bathymetry |
| cice_kmt.sh | Basically a wrapper around coce_kmt.py. Run in exp. directory. Ex: cice_kmt.sh ../topo/depth_TP5a0.06_01 |
|hycom_bathy_modify.py| Can be used to modify bathymetry in points or in blocks based on input from stdin|


