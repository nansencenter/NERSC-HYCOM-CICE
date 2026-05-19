## Initial conditions

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

> **Note:** This step is handled automatically by `srjob.sh` (see [Running the Model](running.md)). Skip it if you use
> `srjob.sh` to submit jobs (recommended).

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

The open boundaries of a regional ocean model must be forced by time-varying fields (e.g., temperature, salinity, velocity, and optionally BGC tracers) from an external ocean product. The source can be a reanalysis such as GLORYS12, output from another ocean model run (e.g., NEMO), or any other product that covers the domain boundaries. This approach is called **offline nesting**, and is the approach supported by NERSC-HYCOM-CICE: information flows in one direction only, from the external product into the regional model. The regional model has no influence on the boundary conditions it receives. This is in contrast to *online* (two-way) nesting, where the regional model (the inner nest) and an outer nest model run simultaneously and exchange information at their shared boundary. Online nesting is not supported.

Offline nesting is activated in HYCOM-CICE by setting `nestfq > 0` (interval in days between 3D nesting archive reads) or `bnstfq > 0` (interval in days between barotropic nesting archive reads) in `blkdat.input`. Three groups of files are required:

- **Boundary configuration files** — define the open boundary geometry and relaxation coefficients (`ports.input`, `rmu`, `rmutr`).
- **Offline nesting archive files** — the external time-varying boundary conditions, interpolated onto the model grid and stored as HYCOM archive files (`archv.YYYY_DDD_HH.[ab]`). Must span the run period; HYCOM interpolates in time between snapshots.
- **Sponge layers** *(optional)* — spatially varying biharmonic diffusion fields (`thkdf4`, `veldf4`) that damp noise near the open boundaries.



### Boundary configuration files

The following files must be present under `$WORK/<CONFIGNAME>/nest/<IEXPT>/`:

| File | Contents |
|------|----------|
| `ports.input` | Open boundary section definitions |
| `rmu.a` / `rmu.b` | 2D relaxation coefficient field for physics |
| `rmutr.a` / `rmutr.b` | 2D relaxation coefficient field for BGC tracers |

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
e-folding time) use `topo_ports.py`. The script reads `regional.grid.a/b` (grid
coordinates) and `regional.depth.a/b` (bathymetry and land mask), detects open boundary
segments automatically from where ocean cells sit at the domain edge, and writes both
`ports.input` and `rmu.a/b`.

> **Note:** `regional.grid.a/b` and `regional.depth.a/b` are pre-existing configuration
> files in `$WORK/<CONFIGNAME>/topo/`. Generating them for a new grid is a separate
> offline step not covered here.

Create a scratch directory and link in the required files:

```bash
mkdir -p $WORK/<CONFIGNAME>/nest/<IEXPT>/SCRATCH
cd $WORK/<CONFIGNAME>/nest/<IEXPT>/SCRATCH
ln -sf $WORK/<CONFIGNAME>/topo/regional.grid.a .
ln -sf $WORK/<CONFIGNAME>/topo/regional.grid.b .
ln -sf $WORK/<CONFIGNAME>/topo/regional.depth.a .
ln -sf $WORK/<CONFIGNAME>/topo/regional.depth.b .
```

Run the script with the bathymetry file, zone width, and e-folding time as arguments
(TP2 values: width=20 cells, e-folding=20 days):

```bash
python $HOME/NERSC-HYCOM-CICE/pythonlibs/modeltools/modeltools/_old/scripts/topo_ports.py \
    regional.depth 20 20
```

This writes `ports.input.tmp` (review and rename to `ports.input`) and `rmu.a/b`. Copy
`rmu` to `rmutr` if the same parameters apply to BGC tracers:

```bash
mv ports.input.tmp ../ports.input
cp rmu.a rmu.b ../
cp rmu.a ../rmutr.a
cp rmu.b ../rmutr.b
```

::::

### Offline nesting archive files

Archive files from the external product must be present in
`$WORK/<CONFIGNAME>/nest/<IEXPT>/`, covering the full run period. They provide the
time-varying boundary conditions (temperature, salinity, velocity, and optionally BGC
tracers) that the model relaxes toward at the open boundaries. HYCOM reads a new archive
at the interval set by `nestfq` (or `bnstfq`) and interpolates in time between
consecutive snapshots. For TP2, daily files are used.

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
expected format. The script automatically stages one extra day before and after the run period so HYCOM can interpolate boundary
conditions across the full run.

Add `--no-fabm` as a fifth argument to skip BGC archive files. 


::::{dropdown} Generating nesting archive files from an external product

For TP2, pre-generated files are available at the path above. If you need to generate
nesting archive files from scratch — for a different time period, source product, or
configuration — use the four-step workflow below. It interpolates fields from an external
ocean product (e.g., GLORYS12 or a NEMO run) onto the inner model grid. Further detail
can be found in the [HYCOM User Guide](https://www.hycom.org/hycom/documentation) and the
[HYCOM examples wiki](https://github.com/HYCOM/HYCOM-examples/wiki/GOMb0.08).

All scripts are run from the source region directory (e.g.,
`$HOME/NERSC-HYCOM-CICE/NMOa0.08/expt_01.1/` for a native NEMO grid,
`$HOME/NERSC-HYCOM-CICE/NMOb0.08/expt_01.0/` for a regular grid).

**Step 1 — Create interpolation matrix** *(once per source/destination pair)*

Computes the horizontal grid mapping from the source (NEMO/GLORYS) grid to the
destination (HYCOM) grid:

```bash
$HOME/NERSC-HYCOM-CICE/bin/isuba_gmapi.sh $WORK/<CONFIGNAME>/
```

Output: `subregion/<IEXPT>/<CONFIGNAME>.gmap.[ab]`

**Step 2 — Interpolate fields** *(repeat for every input file)*

Horizontally and vertically interpolates the external fields onto the HYCOM grid and
isopycnal vertical coordinate:

```bash
# Native GLORYS/NEMO grid, physics only:
$HOME/NERSC-HYCOM-CICE/bin/nemo_to_hycom.sh \
    -d $WORK/<CONFIGNAME>/expt_<EXPT_ID>/ \
    -n /path/to/input/ext-GLORYS12V1_1dAV_20070302_...nc \
    -g native

# Regular grid, with BGC:
$HOME/NERSC-HYCOM-CICE/bin/nemo_to_hycom.sh \
    -d $WORK/<CONFIGNAME>/expt_<EXPT_ID>/ \
    -n /path/to/input/MERCATOR-PHY-24-2018-01-01-12.nc \
    -g regular \
    -b /path/to/input/global_analysis_forecast_bio_20180101.nc
```

| Option | Description |
|--------|-------------|
| `-d` | Destination experiment path (mandatory) |
| `-n` | Path/pattern of input NetCDF files (mandatory) |
| `-g` | Grid type: `native` or `regular` |
| `-b` | BGC input file; also activates BGC boundary creation |
| `-i` | Search radius for wet-point lookup (default: 50 grid cells) |

Output: `$WORK/<CONFIGNAME>/nest/<IEXPT>/archv.YYYY_DDD_HH.[ab]`

**Step 3 — Fix Montgomery potential** *(repeat for every archive file)*

The interpolated archive files have zero Montgomery potential, which must be corrected
for each file. The script uses a restart file from the destination experiment only as a
reference for the bottom boundary condition — any available restart file works. Run from
the destination experiment directory:

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
python $HOME/NERSC-HYCOM-CICE/bin/calc_montg1.py \
    ../nest/<IEXPT>/archv.YYYY_DDD_00.a \
    ./data/restart.YYYY_DDD_00_0000.a \
    ../nest/<IEXPT>/Montg/
```

**Step 4 — Create port and relaxation files**

See [Boundary configuration files](#boundary-configuration-files) for how to generate
`ports.input` and `rmu.[ab]`. Alternatively, `nest_setup_ports.sh` can generate both
in one step with explicit zone width and e-folding time arguments:

```bash
$HOME/NERSC-HYCOM-CICE/bin/nest_setup_ports.sh 20 20
```

::::

### Sponge layers

A common use of nesting is to add a sponge layer near the open boundaries: enhanced
biharmonic diffusion that smooths the solution in the relaxation zone. This is done by
setting `thkdf4` and `veldf4` to negative values in `blkdat.input`, which tells the model
to read spatially varying 2D fields from file rather than using a uniform scalar. The
following files must then be present under `$WORK/<CONFIGNAME>/relax/<IEXPT>/`:

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

This step is handled by `create_ref_case.sh`: a grab-bag script that generates

- relaxation climatologies
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
| `relax/<IEXPT>/` | Physics relaxation climatology (z-level and hybrid-level) | always |
| `relax/<IEXPT>/` | BGC relaxation climatology (z-level and hybrid-level) | `ntracr ≠ 0` |
| `relax/<IEXPT>/` | CO2 relaxation climatology | `ntracr ≠ 0` |
| `relax/<IEXPT>/` | Relaxation mask | always |
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

