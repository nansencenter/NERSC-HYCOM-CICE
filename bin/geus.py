"""GEUS Greenland solid ice discharge per flux gate, for the GloFAS and AHYPE v4b river forcing.

Reads the monthly gate discharge written by geus_preprocess.py (GEUS ice discharge, Mankoff et al.
2020, doi:10.22008/promice/data/ice_discharge). The gates are added to the rivers as extra point
sources and spread like them. The discharge is solid ice, added as liquid freshwater (1 Gt = 1e9 m3).
"""
import sys

import numpy as np
import pandas as pd
import xarray as xr

# Before the first GEUS month (May 1986), use the mean seasonal cycle of these years
CLIMATOLOGY_YEARS = (1987, 1996)


def load_geus(path):
    """Gate positions, discharge per month and the GEUS version, from geus_preprocess.py output."""
    try:
        ds = xr.open_dataset(path)
    except FileNotFoundError:
        sys.exit(f'Error: GEUS file {path} not found. Create it with bin/geus_preprocess.py, '
                 'or set geus_file = None in hycom_river_hifre.py to leave out GEUS.')
    version = 'Dataverse {}, product version {}'.format(ds.attrs.get('geus_dataverse_version', 'unknown'),
                                                        ds.attrs.get('geus_product_version', 'unknown'))
    print(f'GEUS ice discharge {version} ({ds.attrs.get("geus_source", "")}), {path}')

    q = ds['discharge'].transpose('time', 'gate').to_pandas()     # m3/s, monthly means
    q.index = q.index.to_period('M')
    q = q.fillna(0.0)                                             # gates without data in a month
    y0, y1 = CLIMATOLOGY_YEARS
    clim_years = q[(q.index.year >= y0) & (q.index.year <= y1)]
    climatology = clim_years.groupby(clim_years.index.month).mean()
    return dict(lon=ds['lon'].values, lat=ds['lat'].values, monthly=q, climatology=climatology,
                version=version)


def geus_discharge(geus, date):
    """Discharge (m3/s) of every gate on `date`: the mean of its month, or the mean seasonal cycle
    of CLIMATOLOGY_YEARS before the GEUS record starts."""
    month = pd.Timestamp(date).to_period('M')
    q = geus['monthly']
    if month in q.index:
        return q.loc[month].values
    if month < q.index[0]:
        return geus['climatology'].loc[month.month].values
    sys.exit(f'Error: {month} is after the end of the GEUS record ({q.index[-1]})')
