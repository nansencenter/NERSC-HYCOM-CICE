"""Download GloFAS v4.0 daily river discharge to extend the GloFAS input data set.

Downloads the GloFAS v4.0 reanalysis ("cems-glofas-historical", consolidated) daily mean
discharge (dis24) from the Copernicus Early Warning Data Store (EWDS), one request per
month, and makes each month match the existing data_all_year.nc (same grid, 40-90N,
same dimension and variable names), so the months can be appended with cdo mergetime:

   python GloFAS_download.py 2024-07 2024-12 ./glofas_new
   cdo mergetime data_all_year.nc ./glofas_new/data_*.nc data_all_year_extended.nc

Credentials are read by cdsapi from ~/.cdsapirc, which must contain the EWDS url and
your personal API key (from your EWDS profile page):

   url: https://ewds.climate.copernicus.eu/api
   key: <your-personal-access-token>

You also need to accept the dataset licence on the EWDS web page once.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

DATASET = "cems-glofas-historical"
REFERENCE = "/cluster/projects/nn9481k/GloFAS_data/data_v40/data_all_year.nc"
# A little larger than the reference area (40-90N); cut to the reference grid afterwards
AREA = [90, -180, 39.9, 180]  # North, West, South, East


def request(year: int, month: int) -> dict:
    """Return the EWDS request for all days of one month."""
    return {
        "system_version": ["version_4_0"],
        "hydrological_model": ["lisflood"],
        "product_type": ["consolidated"],
        "timespan": ["time_mean"],
        "variable": ["average_river_discharge_in_the_last_24_hours"],
        "year": [f"{year}"],
        "month": [f"{month:02d}"],
        "day": [f"{day:02d}" for day in range(1, 32)],  # days a month does not have are ignored
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": AREA,
    }


def match_reference(raw_file: str, out_file: str, reference: str) -> None:
    """Cut a downloaded month to the reference grid, check it and write it as out_file.

    Renames the time dimension to 'time' if needed (newer EWDS NetCDF files use
    'valid_time'), keeps only dis24 and stops if the grid does not match the reference.
    """
    ref = xr.open_dataset(reference)
    ds = xr.open_dataset(raw_file)
    if "time" not in ds.dims and "valid_time" in ds.dims:
        ds = ds.drop_vars("time", errors="ignore").rename({"valid_time": "time"})
    if "dis24" not in ds.data_vars:
        if len(ds.data_vars) != 1:
            sys.exit("Error: %s has no dis24 and more than one variable: %s" % (raw_file, list(ds.data_vars)))
        ds = ds.rename({list(ds.data_vars)[0]: "dis24"})
    ds = ds[["dis24"]].reset_coords(drop=True)
    ds = ds.sel(latitude=ref.latitude, longitude=ref.longitude, method="nearest", tolerance=1e-3)
    for c in ("latitude", "longitude"):
        if not np.allclose(ds[c].values, ref[c].values, atol=1e-6):
            sys.exit("Error: %s of %s does not match %s" % (c, raw_file, reference))
    days = pd.DatetimeIndex(ds.time.values)
    if (np.diff(days.values) != np.timedelta64(1, "D")).any():
        sys.exit("Error: %s does not contain consecutive days" % raw_file)
    ds["time"].encoding = {k: v for k, v in ref["time"].encoding.items() if k in ("units", "calendar", "dtype")}
    ds["dis24"].encoding = {"dtype": "float32", "zlib": False}
    ds.to_netcdf(out_file)
    print("Wrote %s: %s to %s, %d days" % (out_file, days[0].date(), days[-1].date(), days.size))


def main() -> None:
    """Parse arguments, download each month and match it to the reference file."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("first", help="first month, YYYY-MM")
    p.add_argument("last", help="last month (included), YYYY-MM")
    p.add_argument("outdir", help="directory for the monthly files data_YYYY_MM.nc")
    p.add_argument("--reference", default=REFERENCE, help="existing GloFAS file to match (default: %(default)s)")
    a = p.parse_args()

    if not os.path.isfile(os.path.expanduser("~/.cdsapirc")):
        sys.exit("Error: ~/.cdsapirc not found. Create it with the EWDS url and your API key (see --help).")
    import cdsapi

    os.makedirs(a.outdir, exist_ok=True)
    client = cdsapi.Client()
    for month in pd.period_range(a.first, a.last, freq="M"):
        out_file = os.path.join(a.outdir, "data_%04d_%02d.nc" % (month.year, month.month))
        if os.path.isfile(out_file):
            print("Skipping %s (exists)" % out_file)
            continue
        raw_file = out_file + ".download"
        client.retrieve(DATASET, request(month.year, month.month)).download(raw_file)
        match_reference(raw_file, out_file, a.reference)
        os.remove(raw_file)


if __name__ == "__main__":
    main()
