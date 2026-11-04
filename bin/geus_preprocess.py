#!/usr/bin/env python3
"""Monthly GEUS Greenland solid ice discharge per flux gate, for river_ahype_v4b.py.

Reads gate.nc of the GEUS ice discharge product (Mankoff et al. 2020, doi:10.5194/essd-12-1367-2020;
data doi:10.22008/promice/data/ice_discharge), which gives the discharge (Gt/yr) of ~270 gates near
the marine-terminating glacier fronts at irregular dates (about monthly before 2009, near daily
2009-2024, 12-daily since). The rates are interpolated linearly in time to daily values and averaged
per calendar month; only complete months are kept.

The GEUS Dataverse version of the download (e.g. V131, shown on the dataset page) is not in the
files, so it must be given with --dataverse-version. It is stored in the output together with the
product version and git commit from gate.nc, and printed by river_ahype_v4b.py.

The output is independent of the model grid: river_ahype_v4b.py maps each gate to the nearest coastal
ocean cell and spreads it like the AHYPE v4b outlets. The discharge is converted to a freshwater
volume flux with 1 Gt = 1e9 m3.

Usage:
   python geus_preprocess.py --dataverse-version V131 <GEUS directory with gate.nc> <output .nc>
"""
import argparse
import sys

import numpy as np
import pandas as pd
import xarray as xr

SECONDS_PER_YEAR = 365.25 * 86400


def monthly_means(discharge: pd.DataFrame) -> pd.DataFrame:
    """(time, gate) rates at irregular dates -> monthly means of the daily linear interpolation."""
    daily = discharge.resample('D').mean().interpolate('time')
    counts = daily.iloc[:, 0].resample('MS').count()
    monthly = daily.resample('MS').mean()
    complete = counts == monthly.index.days_in_month
    return monthly[complete.values]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('geusdir', help='directory of the GEUS ice discharge product (with gate.nc)')
    p.add_argument('out', help='output NetCDF file')
    p.add_argument('--dataverse-version', required=True,
                   help='GEUS Dataverse version of the downloaded files, e.g. V131')
    a = p.parse_args()

    ds = xr.open_dataset(f'{a.geusdir}/gate.nc')
    version = ds.attrs.get('product_version')
    if version is None:
        sys.exit('Error: gate.nc has no product_version attribute')
    discharge = ds['discharge'].transpose('time', 'gate').to_pandas()
    if discharge.isna().any().any():
        sys.exit('Error: gate.nc has missing discharge values')
    monthly = monthly_means(discharge)
    total = monthly.sum(axis=1)

    out = xr.Dataset(
        {'discharge': (('time', 'gate'), (monthly.values * 1e9 / SECONDS_PER_YEAR).astype('f8'),
                       {'units': 'm3 s-1', 'long_name': 'solid ice discharge as freshwater volume flux, monthly mean'}),
         'lon': ('gate', ds['mean_lon'].values.astype('f8'), {'units': 'degrees_east'}),
         'lat': ('gate', ds['mean_lat'].values.astype('f8'), {'units': 'degrees_north'}),
         'region': ('gate', ds['region'].values.astype('U2'),
                    {'long_name': 'Mouginot et al. (2019) region: NO North, NE Northeast, CE Central East, '
                                  'SE Southeast, SW Southwest, CW Central West, NW Northwest'}),
         'name': ('gate', ds['name_Mouginot'].values.astype('U32'))},
        coords={'time': monthly.index, 'gate': ds['gate'].values},
        attrs={'title': 'GEUS Greenland solid ice discharge per gate, monthly means',
               'source': 'GEUS ice discharge, doi:10.22008/promice/data/ice_discharge',
               'geus_dataverse_version': a.dataverse_version,
               'geus_product_version': version,
               'geus_source': ds.attrs.get('source', 'unknown'),
               'geus_time_range': f'{discharge.index[0]:%Y-%m-%d} to {discharge.index[-1]:%Y-%m-%d}',
               'method': 'daily linear interpolation in time, monthly means of complete months; 1 Gt = 1e9 m3',
               'history': 'created by NERSC-HYCOM-CICE bin/geus_preprocess.py'})
    out.to_netcdf(a.out)
    print(f'GEUS ice discharge Dataverse {a.dataverse_version}, product version {version} '
          f'({out.attrs["geus_source"]}), '
          f'data {out.attrs["geus_time_range"]}')
    print(f'Wrote {a.out}: {out.sizes["gate"]} gates, {monthly.index[0]:%Y-%m} to {monthly.index[-1]:%Y-%m}, '
          f'total {total.min():.0f}-{total.max():.0f} Gt/yr (mean {total.mean():.0f})')


if __name__ == '__main__':
    main()
