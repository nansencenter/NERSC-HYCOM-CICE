## Initial conditions

Two initialization modes are available, selected via `INITFLG` in `srjob.sh`
(see [Submit a job](running.md#submit-a-job)):

1. **Restart run** (`INITFLG=""`): the model continues from HYCOM and CICE restart files
   valid at the start date. This is the standard mode for hindcast and forecast runs.
   See [Restart files](#restart-files) below for how to obtain them.

2. **Climatological initialization** (`INITFLG="--init"`): temperature and salinity (T/S)
   are read from `relax/<IEXPT>/relax_tem.[ab]` and `relax_sal.[ab]`; layer thicknesses
   are computed internally from the T/S profiles; velocities and sea surface height (SSH)
   start at zero. The start date must be in September (Arctic sea ice is at its annual
   minimum in September, making it a natural starting point). Expect a multi-year
   spin-up before the circulation is reliable. See
   [Climatologies and river forcing](#climatologies-and-river-forcing) for how to
   prepare the required files.

   > **Note:** The `relax_*` file names reflect their use as targets for climatological
   > relaxation during the run (see
   > [Climatologies and river forcing](#climatologies-and-river-forcing)). During
   > initialization, no relaxation is applied. The files are simply read once as the
   > initial T/S state.

### Restart files

HYCOM and CICE each need a restart file valid at the start date of the run
(consistent with `START` in the [srjob.sh variable table](running.md#submit-a-job)).
The naming conventions are as follows:

| Model | File | Notes |
|-------|------|-------|
| HYCOM | `restart.YYYY_DDD_00_0000.[a,b]` | `DDD` is the three-digit day of year (zero-padded) |
| CICE | `iced.YYYY-MM-DD-00000.nc` | Standard calendar date |

For TP2 hindcast runs, restart files are archived at:

```
/nird/datalake/NS9481K/shuang/TP2_output/expt_02.6/
```

For example, to start on 27 August 2016 (240th day of year):

```bash
mkdir -p $WORK/<CONFIGNAME>/expt_<EXPT_ID>/data/cice
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>/data

cp /nird/datalake/NS9481K/shuang/TP2_output/expt_02.6/restart/restart.2016_240_00_0000.a .
cp /nird/datalake/NS9481K/shuang/TP2_output/expt_02.6/restart/restart.2016_240_00_0000.b .

cd cice
cp /nird/datalake/NS9481K/shuang/TP2_output/expt_02.6/cice/iced.2016-08-27-00000.nc .
```

## Atmospheric forcing

This step is handled by the script `atmo_synoptic.sh`.

> **Note:** The job script `srjob.sh` calls `atmo_synoptic.sh` automatically (see [Submit a job](running.md#submit-a-job)). Skip this section if you
> submit via `srjob.sh` (recommended). You may still want to run `atmo_synoptic.sh`
> manually to verify the input files are in place before submitting; if so, comment out
> the `atmo_synoptic.sh` call in `srjob.sh` to prevent the job from regenerating the
> files (output filenames have no dates, so they would be silently overwritten).

Atmospheric forcing must be prepared for each run period. `START` and `END` are the run
start and end times, see the [srjob.sh variable table](running.md#submit-a-job) for the
expected format. To prepare it manually, first activate the
[Python environment](installation.md#use-the-environment), then run:

::::{dropdown} Activating the Python environment on Betzy

```{include} _snippets/betzy_python_activate.md
```

::::

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
$HOME/NERSC-HYCOM-CICE/bin/atmo_synoptic.sh era5+lw $START $END
```

::::{dropdown} What atmo_synoptic.sh does

`atmo_synoptic.sh` reads atmospheric reanalysis data from NetCDF files and converts it to
HYCOM's `.a/b` binary format, interpolated onto the model grid. It calls `hycom_atmfor.py`
internally, which requires the `hycom-cice` Python environment.

The first argument selects the forcing dataset. `era5+lw` uses ERA5 6-hourly data with
downwelling longwave radiation (STRD) included directly. Other supported options are
`erai`, `noresm`, and `ec_op`.

The input data paths are configured in XML files under `$NHCROOT/input/` (e.g.
`era5.xml`). On Betzy, the ERA5 data is at `/cluster/projects/nn9481k/ERA5_6h/`. This
path can be overridden by setting `ERA5_PATH` in the environment before running the
script.

The script must be run from the experiment directory. It reads `EXPT.src` and
`REGION.src` to determine output paths and grid settings.

::::

Running `atmo_synoptic.sh` creates the following files under `$WORK/<CONFIGNAME>/force/synoptic/<IEXPT>/`:

| File | Contents |
|------|----------|
| `airtmp.[ab]` | 2 m air temperature |
| `dewpt.[ab]` | 2 m dew point temperature |
| `mslprs.[ab]` | Mean sea level pressure |
| `precip.[ab]` | Total precipitation rate |
| `radflx.[ab]` | Net longwave radiation at the surface |
| `shwflx.[ab]` | Downwelling shortwave radiation |
| `vapmix.[ab]` | Water vapour mixing ratio |
| `wndewd.[ab]` | Eastward (u) wind component at 10 m |
| `wndnwd.[ab]` | Northward (v) wind component at 10 m |

Since output goes to an experiment-specific directory (`force/synoptic/<IEXPT>/`), running
`atmo_synoptic.sh` from different experiments within the same configuration will not
cause conflicts.

Since forcing files do not include dates in their filenames, leave a stamp in the
directory to record what period they cover:

```bash
cd $WORK/<CONFIGNAME>/force/synoptic/<IEXPT>
rm -f stamp_*
touch stamp_${START}-${END}
```

## Open boundary forcing

The open boundaries of a regional ocean model must be forced by time-varying fields (e.g., temperature, salinity, velocity, and optionally BGC tracers) from an external ocean product. The source can be a reanalysis such as GLORYS12 or output from another ocean model simulation that covers the domain boundaries. This approach is called **offline nesting**, and is the approach supported by NERSC-HYCOM-CICE: information flows in one direction only, from the external product into the regional model. The regional model has no influence on the boundary conditions it receives. This is in contrast to *online* (two-way) nesting, where the regional model (the inner nest) and an outer nest model run simultaneously and exchange information at their shared boundary. Online nesting is not supported here.
Further detail on offline nesting in HYCOM can be found in the [HYCOM User Guide](https://www.hycom.org/hycom/documentation) and the [HYCOM examples wiki](https://github.com/HYCOM/HYCOM-examples/wiki/GOMb0.08).

Offline nesting is activated in HYCOM-CICE by setting `nestfq > 0` (interval in days between 3D nesting archive reads) or `bnstfq > 0` (interval in days between barotropic nesting archive reads) in `blkdat.input`. Three groups of files are required:

- **Boundary configuration files** — define the open boundary geometry and relaxation coefficients for the boundary nudging zone (`ports.input`, `rmu`, `rmutr`).
- **Offline nesting archive files** — pre-interpolated boundary conditions from the external product, stored in HYCOM's `archv` format. Must span the run period; HYCOM interpolates in time between snapshots.
- **Sponge layers** *(optional)* — spatially varying biharmonic diffusion fields (`thkdf4`, `veldf4`) that damp noise near the open boundaries.



### Boundary configuration files

The following files must be present under `$WORK/<CONFIGNAME>/nest/<IEXPT>/`:

| File | Contents |
|------|----------|
| `ports.input` | Open boundary section definitions |
| `rmu.a` / `rmu.b` | 2D relaxation coefficient field for physics nudging at open boundaries | 
| `rmutr.a` / `rmutr.b` | 2D relaxation coefficient field for BGC tracer nudging at open boundaries |

For TP2, you can copy these files from the following reference experiment:

```bash
DIR_NST=/nird/datalake/NS9481K/shuang/nest/TP2_expt023
mkdir -p $WORK/<CONFIGNAME>/nest/<IEXPT>
cd $WORK/<CONFIGNAME>/nest/<IEXPT>
cp ${DIR_NST}/ports.input .
cp ${DIR_NST}/rmu.a ${DIR_NST}/rmu.b .
cp ${DIR_NST}/rmutr.a ${DIR_NST}/rmutr.b .
```

::::{dropdown} Open boundary sections (`ports.input`)

Each port in `ports.input` defines one open boundary segment. The file lists `nports` followed by five
values per port:

| Field | Description |
|-------|-------------|
| `kdport` | Orientation: `1`=N, `2`=S, `3`=E, `4`=W |
| `ifport` | First i-index of the boundary segment |
| `ilport` | Last i-index (`= ifport` for N or S boundaries) |
| `jfport` | First j-index of the boundary segment |
| `jlport` | Last j-index (`= jfport` for E or W boundaries) |

For N/S boundaries the segment spans a range of i-indices at a fixed j; for E/W
boundaries it spans a range of j-indices at a fixed i.

```{figure} TP2_domain.png
:alt: TP2 domain (400×380 grid)

TP2 domain (400×380 grid). Gray areas are the land mask. The two open boundaries are the
Pacific boundary (past the Bering Strait, northern edge) and the Atlantic boundary
(southern edge). Technically a third open boundary exists at the western portion of the
northwest corner, but it has been closed (set to land) for simplicity — visible as a
small stripe of land in the figure.
```

Example `ports.input` for TP2 (two open boundary passages):

```
2  'nports'
1  'kdport'   5  'ifport'   102  'ilport'   380  'jfport'   380  'jlport'
2  'kdport'   174  'ifport'  348  'ilport'     2  'jfport'     2  'jlport'
```

- **Port 1** (`kdport=1`, North): along the northern edge of the domain (`j=380`), spanning
  the western portion of the grid (`i=5–102`). This is the Pacific open boundary, past the
  Bering Strait.
- **Port 2** (`kdport=2`, South): along the southern edge of the domain (`j=2`), spanning
  the central portion of the grid (`i=174–348`). This is the Atlantic open boundary.

At these boundaries the model reads fields like temperature, salinity, and velocity from the
parent model and relaxes toward them according to the `rmu` field.

::::

::::{dropdown} Relaxation coefficients (`rmu` / `rmutr`)

`rmu` and `rmutr` are 2D fields (μ, s⁻¹) of the size of the grid, controlling how strongly the model is nudged
toward the nesting boundary conditions at each grid point: 0 in the interior, ramping
linearly to 1/(e-folding time) at the boundary. `rmu` applies to physics (e.g., temperature,
salinity, velocity); `rmutr` applies to BGC tracers.

HYCOM uses paired `.a`/`.b` files throughout. The `.b` file is a plain-text header
containing variable names, grid dimensions, and min/max summaries: it describes the
contents but holds no spatial data. The `.a` file holds the actual binary data: the full
2D field as 32-bit floats, with a value at every grid point. The two
files must always be kept together.

Example `.b` file for `rmu` and `rmutr` (TP2, zone width=20 cells, e-folding time=20 days;
both files are identical for TP2):

```
Relaxation mask
Relaxation mask created by topo_ports.py. rel zone width=20, efold time=20 days

i/jdm =  400   380
 rmu:  min,max =               0     5.787037e-07
```

Note that the internal variable name is `rmu` in both files. The max value
(5.787037×10⁻⁷ s⁻¹) equals 1/(20 days × 86400 s/day). The `.a` file contains the full
400×380 array, one float per grid point.

::::

::::{dropdown} Generating nesting files from scratch

For TP2 the files can simply be copied from a reference experiment (see above). If you
need to generate them from scratch (e.g. to adjust the relaxation zone width or
e-folding time), run `nest_setup_ports.sh` from the experiment directory with the zone
width and e-folding time as arguments (TP2 values: width=20 cells, e-folding=20 days):

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
$HOME/NERSC-HYCOM-CICE/bin/nest_setup_ports.sh 20 20
```

The script reads `regional.grid.a/b` and `regional.depth.a/b` from `topo/`, detects
open boundary segments automatically from where ocean cells sit at the domain edge, and
writes `ports.input` and `rmu.a/b` to `nest/<IEXPT>/`.

Copy `rmu` to `rmutr` if the same relaxation parameters apply to BGC tracers:

```bash
cd $WORK/<CONFIGNAME>/nest/<IEXPT>
cp rmu.a rmutr.a
cp rmu.b rmutr.b
```

> **Note:** `regional.grid.a/b` and `regional.depth.a/b` are pre-existing configuration
> files in `$WORK/<CONFIGNAME>/topo/`. Generating them for a new grid is a separate
> offline step not covered here.

::::

### Offline nesting archive files

HYCOM uses `archv` files as its standard model state snapshot format: each file stores
the full model state (3D fields such as temperature, salinity, velocity, and layer
thickness, plus 2D fields such as sea surface height and Montgomery potential) at a
single point in time on the model grid and hybrid vertical coordinate system. In the nesting
context, the same format is used for the pre-interpolated boundary conditions from the
external product. The fields cover the entire inner model domain (not just the boundary
cells) because nudging is applied volumetrically over the boundary nudging zone defined
by `nest/rmu`.

HYCOM reads a new archive at the interval set by `nestfq` (or `bnstfq`) and
interpolates in time between consecutive snapshots. For TP2, daily files are used.
These files must be present in `$WORK/<CONFIGNAME>/nest/<IEXPT>/`, covering the full run
period. The expected file names are:

| File | Contents |
|------|----------|
| `archv.YYYY_DDD_HH.a` / `.b` | Physics boundary conditions (`DDD` = zero-padded day of year, `HH` = hour, typically `00`) |
| `archv_fabm.YYYY_DDD_HH.a` / `.b` | BGC boundary conditions (only needed when running with FABM) |

Nesting files for TP2 are archived at `/nird/datalake/NS9481K/shuang/nest/TP2_expt023/`.
Use `stage_nesting_files.sh` to copy or extract the files needed for a given date range:

```bash
$HOME/NERSC-HYCOM-CICE/bin/stage_nesting_files.sh \
    /nird/datalake/NS9481K/shuang/nest/TP2_expt023 \
    $WORK/<CONFIGNAME>/nest/<IEXPT> \
    <START> \
    <END>
```
Here, `START` and `END` are the run start and end times, see the [srjob.sh variable table](running.md#submit-a-job) for the
expected format. Two optional flags are supported (in any order):

| Flag | Effect |
|------|--------|
| `--no-fabm` | Skip BGC (FABM) archive files |
| `--skip-existing` | Skip dates where all expected files are already present in the destination directory |


:::::{dropdown} Generating nesting archive files from GLORYS12/CMEMS

For TP2, pre-generated files are available at the path above. If you need to generate
nesting archive files from scratch (for a different time period, source product, or
configuration) use the workflow below. Two products from the Copernicus Marine Environment Monitoring Service (CMEMS) are used as source data, both
on a regular lat/lon grid:

- **Physics**: GLORYS12 reanalysis (CMEMS `GLOBAL_MULTIYEAR_PHY_001_030`)
- **BGC**: CMEMS global ocean biogeochemistry reanalysis/forecast (`GLOBAL_REANALYSIS_BIO_001_029` for historical, `GLOBAL_ANALYSIS_FORECAST_BIO_001_029` for near-real-time)

**Step 0 — Copy source region files** *(once)*

The source region directory `NMOb0.08/` needs two things in place before any step can run:

- **`NMOb0.08/REGION.src`, `NMOb0.08/expt_01.0/`** (including `blkdat.input`) — these
  are part of the repository. If you cloned the repo before these were added, pull the
  latest version from the `develop` branch:
  ```bash
  cd $HOME/NERSC-HYCOM-CICE && git pull origin develop
  ```
- **`NMOb0.08/topo/`** — grid and bathymetry files that are too large for the repository.
  Both Step 1 and Step 2 read from this directory, so copy them once from NIRD:
  ```bash
  cp -r /nird/datapeak/NS9481K/MERCATOR_DATA/NMOb0.08/topo \
      $HOME/NERSC-HYCOM-CICE/NMOb0.08/
  ```

**Step 1 — Create interpolation matrix** *(once per source/destination pair)*

This step reads the grid coordinates of both the source (GLORYS12) and destination (HYCOM)
regions and uses a conjugate gradient method to find, for each point on the destination
grid, the corresponding fractional indices on the source grid. The result is a mapping
array stored in `gmap.[ab]` that Step 2 uses to horizontally interpolate every field.
This step only depends on the grid geometry, not the ocean fields, so it only needs to
be run once per source/destination grid pair.

```bash
cd $HOME/NERSC-HYCOM-CICE/NMOb0.08/expt_01.0
$HOME/NERSC-HYCOM-CICE/bin/isuba_gmapi.sh $WORK/<CONFIGNAME>/
```

`isuba_gmapi.sh` is a bash script that calls a compiled Fortran binary in
`$HOME/NERSC-HYCOM-CICE/hycom/hycom_ALL/hycom_2.2.72_ALL`. No Python environment is
needed. If you followed the [compilation guide](compilation.md#compile-hycom-all),
the binary will already be in place.

Output: `$HOME/NERSC-HYCOM-CICE/NMOb0.08/subregion/TP2a0.10.gmap.[ab]`

**Step 2 — Interpolate fields**

`nemo_to_hycom.sh` horizontally and vertically interpolates the external fields onto
the HYCOM grid and hybrid vertical coordinate, and accepts the following options:

| Option | Description |
|--------|-------------|
| `-d` | Destination experiment path (mandatory) |
| `-n` | Path pattern of physics input NetCDF files (mandatory) |
| `-g` | Source grid type; always set to `regular` for current GLORYS12 and BGC products |
| `-b` | BGC input path pattern; also activates BGC boundary creation (optional) |
| `-i` | Search radius for wet-point lookup (optional, default: 50 grid cells) |
| `-m` | Path to the GLORYS12 mesh file (contains grid, bathymetry, mask); (optional) |
| `-c` | Path to the GLORYS12 coordinates file (contains vertical layer thicknesses `e3t`); required when `-m` is set (optional) |
| `-h` | Print usage information and exit |

The GLORYS12 and BGC files on NIRD follow these naming conventions:

- **Physics** (`-n`): `MERCATOR-PHY-24-YYYY-MM-DD-12.nc` (one file per day, timestamped at 12:00 UTC)
- **BGC** (`-b`): `global_analysis_forecast_bio_YYYYMMDD.nc` (one file per day)

> **Note:** BGC tracers are interpolated onto the vertical layer structure derived from
> the physics fields and must be processed together in a single pass. **`-n` and `-b`
> must always cover the same date**; the script does not verify this, so make sure
> they match.

`nemo_to_hycom.sh` calls both a Python script and a compiled Fortran binary. Before
running, load the HPC environment and activate the Python environment:

:::{dropdown} Loading the HPC environment and Python environment on Betzy

```bash
source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh
```

```{include} _snippets/betzy_python_activate.md
```

:::

In the examples below, omit `-b` for physics-only runs.

*Single day (2018-01-01), from the login node:*

```bash
cd $HOME/NERSC-HYCOM-CICE/NMOb0.08/expt_01.0
$HOME/NERSC-HYCOM-CICE/bin/nemo_to_hycom.sh \
    -d $WORK/<CONFIGNAME>/expt_<EXPT_ID>/ \
    -n "/nird/datapeak/NS9481K/MERCATOR_DATA/PHY/2018/MERCATOR-PHY-24-2018-01-01-12.nc" \
    -g regular \
    -b "/nird/datapeak/NS9481K/MERCATOR_DATA/BIO/DAILY/2018/global_analysis_forecast_bio_20180101.nc"
```

Processing a single day takes approximately 20 minutes on the login node. For
multi-day or multi-year runs, use a compute node instead:

::::{dropdown} Interactive session on a compute node

Compute nodes on Betzy cannot access NIRD, so copy the source files to `$WORK`
from the login node first.

The following example processes January 2018. Adjust the rsync patterns and date
range for your period.

**Copy source files** *(from the login node)*

```bash
mkdir -p $WORK/input/GLORYS12/PHY/2018
rsync -av --include="MERCATOR-PHY-24-2018-01-*.nc" --exclude="*" \
    /nird/datapeak/NS9481K/MERCATOR_DATA/PHY/2018/ \
    $WORK/input/GLORYS12/PHY/2018/
mkdir -p $WORK/input/GLORYS12/BIO/2018
rsync -av --include="global_analysis_forecast_bio_201801*.nc" --exclude="*" \
    /nird/datapeak/NS9481K/MERCATOR_DATA/BIO/DAILY/2018/ \
    $WORK/input/GLORYS12/BIO/2018/
mkdir -p $WORK/input/GLORYS12
rsync -av \
    /nird/datapeak/NS9481K/MERCATOR_DATA/REGULAR_GRID_COORD/GLO_MFC_001_24_MESH.nc \
    /nird/datapeak/NS9481K/MERCATOR_DATA/REGULAR_GRID_COORD/GLO_MFC_001_24_COORD.nc \
    $WORK/input/GLORYS12/
```

**Start an interactive session and run the script**

```bash
srun --nodes=1 --time=00:30:00 --qos=devel --account=nn2993k --pty bash
```

Once the session starts:

```bash
source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh
```

```{include} _snippets/betzy_python_activate.md
```

```bash
cd $HOME/NERSC-HYCOM-CICE/NMOb0.08/expt_01.0
START=2018-01-01
END=2018-01-31
DATE=$START
while [[ "$DATE" < "$END" || "$DATE" == "$END" ]]; do
    YYYY=$(date -d "$DATE" +%Y)
    YYYYMMDD=$(date -d "$DATE" +%Y%m%d)
    $HOME/NERSC-HYCOM-CICE/bin/nemo_to_hycom.sh \
        -d $WORK/<CONFIGNAME>/expt_<EXPT_ID>/ \
        -n "$WORK/input/GLORYS12/PHY/${YYYY}/MERCATOR-PHY-24-${DATE}-12.nc" \
        -g regular \
        -b "$WORK/input/GLORYS12/BIO/${YYYY}/global_analysis_forecast_bio_${YYYYMMDD}.nc" \
        -m "$WORK/input/GLORYS12/GLO_MFC_001_24_MESH.nc" \
        -c "$WORK/input/GLORYS12/GLO_MFC_001_24_COORD.nc"
    DATE=$(date -d "$DATE + 1 day" +%Y-%m-%d)
done
```

::::

::::{dropdown} Submission script (multi-day or multi-year runs)

Compute nodes on Betzy cannot access NIRD, so copy the source files to `$WORK`
from the login node first.

**Step 1 — Copy source files** *(from the login node)*

Adjust the `for` loop to cover all years you need:

```bash
for YYYY in 2018 2019 2020; do
    mkdir -p $WORK/input/GLORYS12/PHY/${YYYY}
    rsync -av /nird/datapeak/NS9481K/MERCATOR_DATA/PHY/${YYYY}/ \
        $WORK/input/GLORYS12/PHY/${YYYY}/
    mkdir -p $WORK/input/GLORYS12/BIO/${YYYY}
    rsync -av /nird/datapeak/NS9481K/MERCATOR_DATA/BIO/DAILY/${YYYY}/ \
        $WORK/input/GLORYS12/BIO/${YYYY}/
done
mkdir -p $WORK/input/GLORYS12
rsync -av \
    /nird/datapeak/NS9481K/MERCATOR_DATA/REGULAR_GRID_COORD/GLO_MFC_001_24_MESH.nc \
    /nird/datapeak/NS9481K/MERCATOR_DATA/REGULAR_GRID_COORD/GLO_MFC_001_24_COORD.nc \
    $WORK/input/GLORYS12/
```

**Step 2 — Submit the job**

Adjust `START`, `END`, and `#SBATCH` settings as needed:

| Period | `START` | `END` |
|--------|---------|-------|
| One month (January 2018) | `2018-01-01` | `2018-01-31` |
| Several months (Jan–Mar 2018) | `2018-01-01` | `2018-03-31` |
| One year (2018) | `2018-01-01` | `2018-12-31` |
| Several years (2015–2017) | `2015-01-01` | `2017-12-31` |

```bash
#!/bin/bash
#SBATCH --job-name=nemo2hycom
#SBATCH --account=nn2993k
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

source ${HOME}/NERSC-HYCOM-CICE/environment/betzy_env.sh
```

```{include} _snippets/betzy_python_activate.md
```

```bash
cd $HOME/NERSC-HYCOM-CICE/NMOb0.08/expt_01.0
START=2018-01-01
END=2018-12-31
DATE=$START
while [[ "$DATE" < "$END" || "$DATE" == "$END" ]]; do
    YYYY=$(date -d "$DATE" +%Y)
    YYYYMMDD=$(date -d "$DATE" +%Y%m%d)
    $HOME/NERSC-HYCOM-CICE/bin/nemo_to_hycom.sh \
        -d $WORK/<CONFIGNAME>/expt_<EXPT_ID>/ \
        -n "$WORK/input/GLORYS12/PHY/${YYYY}/MERCATOR-PHY-24-${DATE}-12.nc" \
        -g regular \
        -b "$WORK/input/GLORYS12/BIO/${YYYY}/global_analysis_forecast_bio_${YYYYMMDD}.nc" \
        -m "$WORK/input/GLORYS12/GLO_MFC_001_24_MESH.nc" \
        -c "$WORK/input/GLORYS12/GLO_MFC_001_24_COORD.nc"
    DATE=$(date -d "$DATE + 1 day" +%Y-%m-%d)
done
```

::::

:::::

Output: `$WORK/<CONFIGNAME>/nest/<IEXPT>/archv.YYYY_DDD_HH.[ab]`

**Step 3 — Fix Montgomery potential** *(always required)*

Step 2 writes `montg1 = 0` in the nesting archives — the interpolation code does not
compute it. The correct value must be computed from the layer thicknesses and requires
`psikk` (the Montgomery potential at the deepest isopycnal interface) and `thkk` (the
virtual bottom layer thickness) from a restart file of the destination model.

Run `calc_montg1.py` once per archive file, from the destination experiment directory:

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
python $HOME/NERSC-HYCOM-CICE/bin/calc_montg1.py \
    ../nest/<IEXPT>/archv.YYYY_DDD_00.a \
    ./data/restart.YYYY_DDD_00_0000.a \
    ../nest/<IEXPT>/Montg/
```

Any restart from the same model configuration works — it does not need to be the
restart used to initialise your particular run. The script writes corrected archive
files to the `Montg/` subdirectory. Repeat for every archive file in the nest directory.

:::{dropdown} What calc_montg1.py does

The script recomputes only the `montg1` (Montgomery potential) field in each archive,
leaving all other fields unchanged.

**What it reads from the archive**

For each layer it reads temperature (`temp`), salinity (`salin`), and layer thickness
(`thknss`). It integrates `thknss` from the surface downward to build the cumulative
pressure profile `p[k]`, and computes the potential density anomaly `thstar[k]` from T
and S using the equation of state selected by `thflag` in `blkdat.input`. It also reads
`srfhgt` (sea surface height).

**What it reads from the restart**

It reads `psikk` and `thkk` — the Montgomery potential and thickness of the virtual
bottom layer, which sits below the deepest active isopycnal layer. For configurations
with a terrain-following bottom layer, both are non-zero. For TP2 (z + isopycnal, no
terrain-following), `thkk = 0` everywhere, but `psikk` is still non-zero — it
represents the bottom Montgomery potential determined by the bathymetry and reference
density structure. Because `psikk` does not depend on the instantaneous ocean state,
any restart from the same configuration produces the correct result.

**How montg1 is computed**

`montg1` is split into two parts:

- `montg1c`: the part independent of mean bottom pressure `pbavg`, assembled from
  `psikk`, `thkk`, and the density/pressure profile.
- `montg1pb`: the part proportional to `pbavg`, assembled from the same profile.

`pbavg` is then solved by requiring the surface value to match `srfhgt`:

```
pbavg = (srfhgt − montg1c) / (montg1pb + thref)
montg1 = montg1pb × pbavg + montg1c
```

:::

### Sponge layers

A common use of nesting is to add a sponge layer near the open boundaries: enhanced
biharmonic diffusion that smooths the solution near the boundary. This is done by
setting `thkdf4` and `veldf4` to negative values in `blkdat.input`, which tells the model
to read spatially varying 2D fields from file rather than using a uniform scalar. The
following files must then be present under `$WORK/<CONFIGNAME>/relax/<IEXPT>/`:

> **Note:** The sponge layer (enhanced diffusion), the boundary nudging zone (`nest/rmu`),
> and the climatological relaxation zone (`relax_rmu`) are three independent mechanisms
> with independently defined spatial extents. They are generated by separate scripts and
> need not cover the same area.

| File | Description |
|------|-------------|
| `thkdf4.[a,b]` | 2D biharmonic diffusion coefficient for layer thickness |
| `veldf4.[a,b]` | 2D biharmonic diffusion coefficient for velocity |

> **Note:** Spatially varying `thkdf4` and `veldf4` fields are not limited to sponge
> layers. Any spatially varying biharmonic diffusion pattern can be prescribed this way.
> The sponge layer is just the most common use case in a nesting context.

These files are generated as follows. First, create a scratch directory and link in the grid and topography files:

```bash
mkdir -p $WORK/<CONFIGNAME>/relax/<IEXPT>/SCRATCH
cd $WORK/<CONFIGNAME>/relax/<IEXPT>/SCRATCH
ln -sf $WORK/<CONFIGNAME>/topo/regional.grid.a .
ln -sf $WORK/<CONFIGNAME>/topo/regional.grid.b .
ln -sf $WORK/<CONFIGNAME>/topo/regional.depth.a .
ln -sf $WORK/<CONFIGNAME>/topo/regional.depth.b .
```

Copy the script locally:

```bash
cp $HOME/NERSC-HYCOM-CICE/bin/hycom_creat_thkdf4_veldf4.py .
```

Open the script and set the parameters for your configuration:

| Parameter | TP2 | TP5 | Description |
|-----------|-----|-----|-------------|
| `mnn` | `0.05` | `0.02` | Minimum diffusion value |
| `mxx` | `0.125` | `0.125` | Maximum diffusion value |

Then activate the [Python environment](installation.md#use-the-environment) and run the script with the sponge layer width as its argument (20 grid cells for TP2):

::::{dropdown} Activating the Python environment on Betzy

```{include} _snippets/betzy_python_activate.md
```

::::


```bash
python hycom_creat_thkdf4_veldf4.py 20
```

Copy the output to the parent directory:

```bash
cp thkdf4.a thkdf4.b veldf4.a veldf4.b ../
```

::::{dropdown} What hycom_creat_thkdf4_veldf4.py does

The script reads the grid (`regional.grid.a/b`) and bathymetry (`regional.depth.a/b`) to
identify open boundary cells, then constructs 2D diffusion coefficient fields that ramp
linearly from a minimum value (`mnn`) in the interior to a maximum value (`mxx`) at the
open boundaries, over the specified number of grid cells. The result is written as
`thkdf4.a/b` (diffusion for layer thickness) and `veldf4.a/b` (diffusion for velocity).

::::

## Climatologies and river forcing

HYCOM has two separate boundary nudging systems that operate simultaneously:

- **Nesting nudging** (`nest/rmu` + nesting archives): nudges T, S, and interface
  heights toward time-varying fields from the parent model (NEMO) near the open
  boundary. See [Offline nesting archive files](#offline-nesting-archive-files).
- **Climatological relaxation** (`relax_rmu` + climatologies): nudges T, S, and interface
  heights toward monthly climatology in a broader zone behind the boundary. Controlled
  by `relax` in `blkdat.input`.

This step is handled by `create_ref_case.sh`, which generates the files needed for
climatological relaxation and climatological initialization:

- T/S and interface climatologies (`relax_sal`, `relax_tem`, `relax_int`) — nudging
  targets for climatological relaxation, and initial T/S state for
  [climatological initialization](forcing.md#initial-conditions)
- Climatological relaxation mask (`relax_rmu`) — 2D field of 1/e-folding times (1/s),
  non-zero where climatological relaxation is active; gates nudging of T, S, and
  interface heights toward the climatologies above
- river forcing
- light attenuation (kpar) forcing
- MPI tile definition files

Copy the script to your experiment directory:

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
cp $HOME/NERSC-HYCOM-CICE/bin/create_ref_case.sh .
```

Open `create_ref_case.sh` and update:

- `Icore` — number of MPI tiles in the i-direction (e.g. `29` for TP2 on Betzy)
- `Jcore` — number of MPI tiles in the j-direction (e.g. `26` for TP2 on Betzy)
- `iceclim` — set to `0` (no ice climatology) or `1` (initialize from ice climatology)

`Icore` and `Jcore` control how the domain is decomposed for MPI parallelism. Choose
them so that `Icore × Jcore × (ocean fraction)` is close to your available core count,
while keeping individual tiles large enough to be efficient (at least ~10 grid cells in
each direction). On Betzy, nodes have 128 cores and jobs should use a multiple of 128
cores. For TP2 (~67% ocean), `Icore=29, Jcore=26` gives `NMPI=504`, which is close to
512 (4 nodes).

Also find the `tile_grid.sh` line in `create_ref_case.sh` and add the `-s 1` flag (tile
size variation factor; `1` encourages near-uniform tiles, default `9.5` allows more
variation), so it reads:

```bash
tile_grid.sh -s 1 -${Icore} -${Jcore} ${T} > $EDIR/log/ref_tiling.out 2>&1
```

Before running, source the HPC environment to make the compiled binaries available.
For BGC runs (`ntracr > 0`), also activate the Python environment (needed for the river
forcing scripts):

::::{dropdown} Source HPC environment — Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_hpc_env.md
```

::::

::::{dropdown} Activate Python environment — Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_python_activate.md
```

::::

Then run the script:

```bash
bash create_ref_case.sh
```

After the script completes, check the tile definition file created in `topo/partit/`. The
four-digit number in the suffix is `NMPI` (the actual number of ocean MPI tasks, less than
`Icore × Jcore` because land-only tiles are excluded). Set this value in `EXPT.src`. For
TP2 on Betzy (`Icore=29`, `Jcore=26`, topography version `04`), this is `NMPI=504`.

::::{dropdown} Files generated by create_ref_case.sh

| Location | Contents | Condition |
|----------|----------|-----------|
| `relax/<IEXPT>/` | T/S and interface climatologies (`relax_sal`, `relax_tem`, `relax_int`; z-level and hybrid-level) — climatological relaxation targets and climatological initialization | always |
| `relax/<IEXPT>/` | Climatological relaxation mask (`relax_rmu`) — 2D field of 1/e-folding times gating nudging toward the climatologies above | always |
| `relax/<IEXPT>/` | BGC relaxation climatology (z-level and hybrid-level) | `ntracr ≠ 0` |
| `relax/<IEXPT>/` | CO2 relaxation climatology | `ntracr ≠ 0` |
| `relax/<IEXPT>/` | Sea ice cover climatology | `iceclim=1` |
| `force/rivers/<IEXPT>/` | River forcing | always |
| `force/seawifs/` | kpar forcing (diffuse light attenuation from SeaWiFS; only read by the model when `jerlv0=0` in `blkdat.input`) | always |
| `topo/partit/` | MPI tile definition files (e.g. `depth_TP2a0.10_04.0504`) | always |

All paths are relative to `$WORK/<CONFIGNAME>/`.

> **Note:** Most output goes to IEXPT-specific paths (`relax/<IEXPT>/`,
> `force/rivers/<IEXPT>/`), so running `create_ref_case.sh` from different experiments
> does not cause conflicts. Two files are shared across all experiments in a
> configuration: the z-level relaxation climatology (`relax/<CLIM_CHOICE>/`) and kpar
> forcing (`force/seawifs/`). These would be overwritten, but with identical content
> since they do not depend on experiment settings. Tile files (`topo/partit/`) are also
> shared; they are named by tile count, so different `Icore`/`Jcore` settings produce
> separate files rather than overwriting, but `NMPI` in `EXPT.src` must match the tile
> file used.
>
> If all files already exist from a previous run of `create_ref_case.sh` (e.g. for a
> different experiment in the same configuration), you can skip rerunning the script
> entirely. In that case, symlink or copy the IEXPT-specific files into the appropriate
> `relax/<IEXPT>/` and `force/rivers/<IEXPT>/` folders for the new experiment.

::::

::::{dropdown} Optional: custom river forcing

The default river forcing generated by `create_ref_case.sh` is sufficient for most
experiments. A more complete river forcing combining AHYPE/EHYPE, TRIP, and Greenland
ice melt can be built manually using the scripts in `$HOME/NERSC-HYCOM-CICE/bin/`:

- **AHYPE/EHYPE** (Arctic/European HYPE — Hydrological Predictions for the Environment):
  a hydrological model providing spatially distributed river discharge. Produced by
  `river_nersc.sh`, which places discharge along the coast using a `rivers.dat` file.
- **TRIP** (Total Runoff Integrating Pathways): routes ERA40/ERAI reanalysis runoff to
  river mouths using a river network database. Produced by `river_trip.sh`.

All scripts must be run from the experiment directory. They write to
`$WORK/<CONFIGNAME>/force/rivers/`, overwriting the default files from
`create_ref_case.sh`.

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
```

Steps:

1. Run `river_nersc.sh` to produce AHYPE/EHYPE river runoff.
2. Run `river_trip.sh` to produce TRIP river runoff.
3. Run `merge_TRIP_Hype_abfile_4Greenlnad.py` to merge all sources into a single
   `rivers.a/b`, including Greenland ice melt.

> **TODO:** Clarify whether `triver` in `blkdat.input` needs to be set to `1` when using
> the merged river forcing that includes Greenland ice melt.

> **TODO:** The Greenland files (`griver.a/b`) are on the TP5 grid (780×800) and must be regridded to TP2 before merging. How?

::::

