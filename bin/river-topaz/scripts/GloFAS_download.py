import os
# CDS API
import cdsapi
# Disable warnings for data download via API
import urllib3
urllib3.disable_warnings()
# Disable xarray runtime warnings
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)


if os.path.isfile("~/.cdsapirc"):
    cdsapi_kwargs = {}
else:
    URL = 'https://ewds.climate.copernicus.eu/api'
    KEY = '2580a944-a2bb-47eb-9a62-7031382adc84'
    cdsapi_kwargs = {
        'url': URL,
        'key': KEY,
    }


DATADIR = './river'
os.makedirs(DATADIR, exist_ok=True)
download_file = f"{DATADIR}/glofas-2012_2022.grib"
if not os.path.isfile(download_file):
    c = cdsapi.Client()
    DATASET = "cems-glofas-historical"
    REQUEST = {
        "system_version": ["version_4_0"],
        "hydrological_model": ["lisflood"],
        "product_type": ["consolidated"],
        "variable": ["river_discharge_in_the_last_24_hours"],
        "hyear": [f"{year}" for year in range(2000, 2023)],
        "hmonth": [f"{month:02d}" for month in range(1, 13)],
        "hday": [f"{day:02d}" for day in range(1, 31)],
        "data_format": "netcdf",
        "download_format": "zip",
        "area": [90, -180, 30, 180]
    }
    c.retrieve(DATASET, REQUEST).download(download_file)