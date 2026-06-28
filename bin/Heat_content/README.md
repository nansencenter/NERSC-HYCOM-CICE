# Ocean Heat Content Diagnostics

Tools for computing and visualising Ocean Heat Content (OHC) from HYCOM/TOPAZ model output and GLORYS reanalysis, over the North Atlantic and Arctic domain.

---

## Files

| File | Description |
|------|-------------|
| `yearly_meanTP2.sh` | Preprocesses raw TOPAZ2 monthly output into annual mean NetCDF files |
| `cal_ohc.py` | Core OHC calculation library |
| `diag_errinfo.py` | Diagnostic utilities (error metrics, NetCDF writing) |
| `plot_ohc.py` | Plotting routines |
| `yearly_ohc_mean.py` | Main script: computes and plots OHC time series for two model runs and GLORYS |
| `NAtlantic_Arctic_regions_generic12_TP2.nc` | Regional masks (regions 1–8) on the TP2 grid |

---

## Workflow

### Step 1 — Prepare annual mean files

Run `yearly_meanTP2.sh` from the directory where you want the output files written.
It reads monthly TOPAZ2 files and creates one annual mean per variable per year,
provided all 12 monthly files are present.

```bash
bash yearly_meanTP2.sh
```

Input files are expected at:
```
/nird/datalake/NS9481K/shuang/seaclim/output/old/3D/PHY/3DT-{var}/
    TOPAZ2_1m-m_3DT-{var}_{YYYYMM}-{YYYYMM}.nc
```

Output files are written to the current directory:
```
TP2thetao_archv{YYYY}.nc
TP2so_archv{YYYY}.nc
```

### Step 2 — Compute and plot OHC time series

```bash
python yearly_ohc_mean.py [--region N] [--depth N] [--anomaly]
```

Reads annual files for two model runs and GLORYS (1993–2024), computes OHC
in J/m² for the selected depth range(s) and region, caches intermediate results
as NetCDF, and produces time-series plots.

---

## Runtime options for `yearly_ohc_mean.py`

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `--region` | 3–8 | 0 | Regional mask to apply. `0` = whole domain. Region indices correspond to those in `NAtlantic_Arctic_regions_generic12_TP2.nc`. Regions 1 and 2 are outside the TOPAZ domain and should not be used. |
| `--depth` | 1–6 | 0 | Depth range to process. `0` = all ranges (see table below). |
| `--anomaly` | flag | off | Also produce anomaly plots (deviation from each dataset's temporal mean). |

### Depth range indices

| `--depth` | Range |
|-----------|-------|
| 1 | Full depth (0–4000 m) |
| 2 | 0–300 m |
| 3 | 0–700 m |
| 4 | 0–2000 m |
| 5 | 700–2000 m |
| 6 | 2000–4000 m |

### Examples

```bash
# All depth ranges, whole domain, total OHC only
python yearly_ohc_mean.py

# 0–700 m, region 3, with anomaly plots
python yearly_ohc_mean.py --depth 3 --region 3 --anomaly

# Full depth, whole domain, with anomaly plots
python yearly_ohc_mean.py --depth 1 --anomaly
```

---

## Output

### Cached OHC fields (NetCDF)
Written to `/cluster/work/users/annettes/OHC_test/`:

```
GOHC_{year}-{depth}[-reg{N}].nc    # GLORYS
NOHC_{year}-{depth}[-reg{N}].nc    # New_ref_run
OOHC_{year}-{depth}[-reg{N}].nc    # Old_ref_run
```

If a cache file already exists it is reused, skipping the OHC calculation for that year.

### Plots (JPG)
Written to the current working directory:

```
OHC_time_yearly_{depth}[-reg{N}].jpg       # Total OHC time series (J/m²)
OHC_anom_yearly_{depth}[-reg{N}].jpg       # Anomaly time series (if --anomaly)
```

---

## Data sources

| Dataset | Path |
|---------|------|
| GLORYS reanalysis | `/nird/datapeak/NS9481K/SEACLIM/GLORYS_annual/` |
| New_ref_run annual files | `/nird/datapeak/NS9481K/SEACLIM/New_ref_run_annual/` |
| Old_ref_run annual files | `/nird/datapeak/NS9481K/SEACLIM/Old_ref_run_annual/` |
| TP2 grid / topography | `/cluster/work/users/annettes/TZ2a0.10/topo/` |

---

## Dependencies

- `numpy`, `xarray`, `netCDF4`
- [`gsw`](https://teos-10.github.io/GSW-Python/) — Gibbs SeaWater (TEOS-10), for in-situ density
- `matplotlib`
- `abfile` — for reading HYCOM binary grid files
