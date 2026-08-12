# MSCPROGS post-processing tools

MSCPROGS is a collection of Fortran (and C) programs for post-processing and
analysing HYCOM model output. It is a general-purpose toolkit covering grid projection, field extraction, temporal
averaging, section transports, Lagrangian ice drift, grid utilities, and more.
The tools live in `hycom/MSCPROGS/`.

All of these post-processing and analysis steps can alternatively be done with [xhycom](https://xhycom.readthedocs.io/en/latest/index.html), a Python-based toolkit for working with HYCOM output.

For compilation instructions, see [Compile MSCPROGS](compilation.md#compile-mscprogs-libhycnersca).

:::{important}
There are two ways to post-process and analyse HYCOM output — pick one:

- **MSCPROGS** (this page) — the original Fortran/C toolkit, described below.
- **[xhycom](xhycom.md)** — a Python-based toolkit covering the same
  post-processing and analysis tasks.

You don't need both. See [xhycom](xhycom.md) if you'd rather work in Python.
:::

## Running the tools

### Environment setup

Before running any MSCPROGS tool, load the HPC modules and set `$MSCPROGS` and `$PATH`.

::::{dropdown} Source HPC environment — Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_hpc_env.md
```

::::

::::{dropdown} Source HPC environment — Olivia (NRIS/Sigma2)

```{include} _snippets/olivia_hpc_env.md
```

::::

```bash
export MSCPROGS=${HOME}/NERSC-HYCOM-CICE/hycom/MSCPROGS
export PATH=${MSCPROGS}/bin:${MSCPROGS}/bin_setup:${PATH}
```

### Working directory

The natural place to run MSCPROGS tools is the experiment data directory:

```bash
CONFIGNAME=<CONFIGNAME>   # e.g. TP2a0.10
EXPT_ID=<EXPT_ID>         # e.g. 01.0
cd $WORK/${CONFIGNAME}/expt_${EXPT_ID}/data/
```

[expt_postprocess.sh](running.md#submit-a-job) copies `regional.grid.*` and `regional.depth.*` here
automatically after each model run, so the grid files are already in place.

You still need to copy the relevant [input files](#input-files) into the working
directory. Example input files are provided in `$MSCPROGS/Input/`.

:::{tip}
You can also run MSCPROGS tools directly from the scratch directory
(`expt_<EXPT_ID>/SCRATCH/`) to analyse output before the run is finalised.
`expt_preprocess.sh` already copies `regional.grid.*` and `regional.depth.*`
there, so only the input files need to be copied in — the same step as above.
:::

### Typical workflow

The general pattern for any MSCPROGS tool is:

1. Navigate to the data directory and copy the needed input files (above).
2. Run the tool — most print usage information when called with no arguments:
   ```bash
   hyc2proj           # prints usage
   m2nc               # prints usage
   m2section          # prints usage
   ```
3. Process the output NetCDF files with your analysis scripts.

## `hyc2proj` and `hyc2stations` — projection and interpolation

`hyc2proj` projects HYCOM isopycnal-layer output onto a regular grid and
interpolates vertically to fixed depth levels, writing CF-compliant NetCDF.
Four output projections are supported:

| Projection | `proj.in` keyword | Notes |
|---|---|---|
| Regular lon/lat | `regular` | Uniform degree spacing; grid defined by start/end/spacing in degrees |
| Mercator | `mercator` | Uniform spacing in projected coordinates (non-uniform in degrees); defined by start/end bounds + number of points; conformal |
| Polar stereographic | `polar_stereographic` | Defined by projection coordinates + number of points; optionally rotate vectors onto the PS grid |
| Model native | `native` | No horizontal interpolation; extract sub-domain by index range |

`hyc2proj` requires the following files in the working directory:

- `regional.grid.a/.b` and `regional.depth.a/.b` — already present in `data/` and `SCRATCH/`
- `grid.info` — conformal mapping parameters for the model grid; copy from `topo/`:
  ```bash
  cp $WORK/<CONFIGNAME>/topo/grid.info .
  ```
- `proj.in` — target projection and grid (see [Input files](#input-files))
- `depthlevels.in` — vertical depth levels to interpolate to
- an `extract.*` file — fields to extract; the tool auto-selects the file by name based on the input file type (e.g. `extract.archm` for archm files), so the name must not be changed

Horizontal interpolation is bilinear; vertical interpolation defaults to cubic
spline. `staircase` and `linear` options are also available (faster but lower
quality).

```bash
cp $WORK/<CONFIGNAME>/topo/grid.info .
cp $MSCPROGS/Input/proj.in.regular_grid proj.in     # choose and edit a sample proj.in
cp $MSCPROGS/Input/depthlevels.in .
cp $MSCPROGS/Input/extract.archm .
hyc2proj archm.*.a                                   # process one or more .a files
```

One NetCDF output file is created per input file. The name is derived from the
input filename with the extension changed to `.nc` and the date embedded — for
example `archm_20100101_00.nc` for an archm file. The original `.a`/`.b` files
are not modified.

:::{note}
`hyc2proj` can be slow on a login node for non-native projections. As a
reference, a single TP2 archm file with 4 three-dimensional fields and 40
target depth levels takes roughly 2 minutes on a regular 0.125° grid covering
the Arctic (−38°–20°E, 65°–88°N). For many files, consider running in an
interactive job or a short batch script.
:::

`hyc2stations` is a companion routine that interpolates HYCOM fields onto
specified spatial positions (station locations or cross-section points) rather
than a grid. It reads `stations.in` instead of `proj.in`, together with
`depthlevels.in` and an `extract.*` file. Output is one NetCDF file per station
group, named `<inputbase>_group<groupname>.nc`. A list of all generated files is
written to `hyc2stations.filelist`. A helper script `setupstations.py` generates
`stations.in` entries from start/end coordinates along a rhumb line.

### Batch post-processing

For a single run or a few files, running `hyc2proj` interactively (above) is fine. For many
files — e.g. a multi-year reanalysis — submit it as a SLURM job. A convenient pattern is to
parallelise by year: stage each year's archives into its own `Y<year>/` directory together
with a copy of the required input files, then run one `hyc2proj` task per year.

The required files (see
[`hyc2proj`](#hyc2proj-and-hyc2stations-projection-and-interpolation)) must be present in
each per-year directory: `regional.grid.*`, `regional.depth.*`, `grid.info`, `proj.in`,
`depthlevels.in`, the matching `extract.*` file, and the compiled `hyc2proj` executable.

:::{note}
The parallel-by-year wrapper does not work with `detvflux` — leave it off in
`extract.archm`/`extract.archv`.
:::

::::{dropdown} Example: parallel-by-year SLURM script (Olivia)

The header (`#SBATCH` options) and module block are machine-specific. On Olivia, use a CPU
partition and load the modules established when
[compiling MSCPROGS](compilation.md#compile-mscprogs-libhycnersca).

```bash
#!/bin/bash
#SBATCH --account=nn2993k
#SBATCH --job-name=hyc2proj
#SBATCH --nodes=1
#SBATCH --ntasks=32            # one task per year
#SBATCH --mem-per-cpu=10G
#SBATCH --time=05:00:00        # ~3 h per model year, as a guide

# --- Olivia HPC modules ---
source $HOME/NERSC-HYCOM-CICE/environment/olivia_env.sh
ulimit -s 2000000

yr_s=1993
yr_e=2024
rid=archm.

# Stage each year's archives + input files into Y<year>/
for year in $(seq $yr_s $yr_e); do
    dir="Y$year"; mkdir -p "$dir"
    cp ./regional.* ./grid.info ./extract.archm ./depthlevels.in ./proj.in ./hyc2proj "$dir"/
    mv ${rid}${year}_*.a ${rid}${year}_*.b "$dir"/ 2>/dev/null
done

# Run one hyc2proj per year in parallel
for year in $(seq $yr_s $yr_e); do
    ( cd "Y$year" && srun --exclusive -n1 -c1 ./hyc2proj --vertint linear ${rid}${year}_*.a ) &
done
wait
```

::::

## `m2nc` / `m2t` — 2D field extraction to NetCDF

Reads HYCOM `.ab` files and writes selected 2D fields to `tmp1.nc`. Fields to
extract are specified in the relevant `extract.*` file. Multiple input files
produce multiple time records in the output. `m2t` additionally writes a
`tmp1.dat` file in [Tecplot](https://www.tecplot.com) ASCII format (a
commercial CFD visualization tool; the format is also readable by ParaView and
similar tools).

```bash
cp $MSCPROGS/Input/extract.archm .
m2nc archm.*.a
```

:::{note}
No vertical interpolation — output stays on isopycnal layers.
:::

## `h2nc` (ExtractNC3D) — 3D field extraction to NetCDF

```bash
h2nc <file(s)>
```

Converts HYCOM `.ab` files to NetCDF with 3D fields. Like `m2nc` it uses the
`extract.*` files to select fields, but retains the full vertical structure.
Output is written to `tmp1.nc`.

## `hycave` — temporal averaging

Computes a thickness-weighted time mean across the supplied files. `filetype` is
one of `restart`, `nersc_daily`, `nersc_weekly`, `archv`, or `archm`. Requires
`regional.grid` and `regional.depth` in the working directory. Output is a new
`.ab` file pair with a fixed placeholder date — for example
`AVEDAILY_9999_999_9999_999.ab` for `nersc_daily` and `AVE.archv.9999_999_99.ab`
for `archv`. Input files are not modified.

```bash
hycave nersc_daily TP2daily*.a
```

## `m2section`, `m2transport`, `m2transport2` — section transports

Extracts data along user-defined sections and computes transports. Sections are
defined in `sections.in` (see [Input files](#input-files)). All three scripts
call `section_intersect` first to identify model grid points along each section
using a great-circle algorithm.

```bash
cp $MSCPROGS/Input/sections.in .
cp $MSCPROGS/Input/extract.archm .              # auto-selected by file type
m2section archm.*.a                            # output: section001.nc, section002.nc, ...
m2transport archm.*.a                          # output: transport_net001.dat, ...
```

- **`m2section`** — extracts fields along sections to `section001.nc`,
  `section002.nc`, … Uses `extract.*` files to select fields.
- **`m2transport`** — fast barotropic volume and ice transport across sections.
  Output to ASCII files `transport_net001.dat`, `transport_net002.dat`, …
- **`m2transport2`** — full-depth transport with finer control. Transport
  criteria (depth range, water mass temperature range, etc.) are set in
  `transport.in`. Also writes `transports2.nc`.

`m2transport2` additionally supports scalar (e.g. heat) transports when
`scalartransport.in` is present. Each line specifies:

```
'<variable>'  <offset>  <scale_factor>
```

For example, heat transport relative to 6 °C:

```
'temp'  6.  3987000
```

where `3987000 ≈ cₚ × ρ`.

:::{note}
Positive transport direction: stand at the first section point and look toward
the second — positive is to the right. Sections follow great circles, not
rhumb lines.
:::

## `icedrift2` — Lagrangian ice drift

```bash
icedrift2            # reads Input/icedrift.in
icedrift2 <args>     # CERSAT mode: takes CERSAT ice drift positions as input
```

Integrates ice particle trajectories forward in time using Runge–Kutta (2nd
order), reading drift velocities from HYCOM daily output. Useful for validation
against IABP buoy data or as input for data assimilation. Output is written to
`drift.uf` and `drift_diag.uf` (in CERSAT mode the names follow the input drift
file, e.g. `<driftfile>.uf` and `<driftfile>_diag.uf`). A text summary is
written to `tsdrift.asc`.

## `nestbat` — nesting bathymetry

Adjusts a local-domain depth file so bathymetry at the nesting boundary is
consistent with the parent (global) grid, smoothly transitioning across a
user-specified zone. The `nestbat-2.2` version is used by the `nestbat.sh` setup
script. Run from an empty working directory containing the required global and
local grid files; the program is interactive.

## `Tides_FES2014` — tidal forcing from FES atlas

Generates tidal boundary conditions for HYCOM from the FES2014 tidal atlas.
Requires the external FES2014 C library and the GNU C compiler; see
`src/Tides_FES2014/README.md` for build instructions.

## Other tools

| Tool | Executable | Purpose |
|---|---|---|
| `Barstrf` | `hyc2bastrf` | Barotropic streamfunction |
| `Conf_grid` | — | Generate conformal-mapping model grids |
| `ConfmapRoutines` | — | Convert between grid indices and lon/lat |
| `Curviint` | `curviint` | Interpolate restart files from one grid to another |
| `DateTools` | — | Convert between ordinal days and calendar dates |
| `DProfile` | `dprofile` | Extract vertical profiles from model files |
| `Ensstat` | — | Ensemble and time-domain statistics |
| `FindLayer` | `findlayer` | Layer thickness between variable thresholds (e.g. salinity 35–36 psu) |
| `GP` / `GPdens` | — | Extract and analyse grid-point time series (activated via `infile.in`) |
| `GridToLL` | `gridtoll` | Interpolate fields from `tmp1.nc` to a regular lon/lat grid |
| `MkEnsemble` | `mkensemble` | Create ensemble by perturbing restart files |
| `Relax` | `rmunew` | Generate relaxation masks |
| `RelaxToNetCDF` | — | Interpolate climatology files to projection/depth levels |
| `River_Forcing` | `rivers` | Spread point-source river discharge along the model coast |
| `SSHFromState` | `ssh_from_restart` | Compute SSH from restart files |
| `Tides_CSR` | `csr2mod_GE` | Tidal boundary forcing from CSR tidal atlas |
| `TRIP` | `trip_*` | River forcing from the TRIP database + ERA40/ERA-i runoff |
| `ZONAL` | `zonal`, `mosf` | Zonal averages and meridional overturning streamfunction |

## Post-processing runs from NIRD

NIRD access differs between machines:

| | Login node | SVC node | Compute node |
|---|---|---|---|
| Betzy | read-write | — | not mounted |
| Olivia | not mounted | read-write | read-only |

The recommended workflow is to **copy** (not symlink) archives onto the scratch filesystem
before processing — the [parallel-by-year `hyc2proj` script](#batch-post-processing) *moves*
files into per-year directories, which fails on a read-only source, and streaming large
3-D archives live over NIRD is slow.

Stage the archive files (`archm.*` / `archv.*`, both `.a` and `.b`) together with the run's
`regional.*` and `grid.info` into a scratch directory:

- **Betzy**: copy from the login node to `/cluster/work/users/$USER/postproc/<run>/`
- **Olivia**: copy from a service (SVC) node to `/cluster/work/projects/nn2993k/$USER/postproc/<run>/`

Then run the post-processing tool from the staging directory. On Olivia, tools that only
read archives — `m2nc`, `hycave`, `m2section` — can also read directly from NIRD inside a
batch job, skipping the staging step.

Move the resulting `.nc` files to `/cluster/projects/nn2993k/$USER/...` (or back to NIRD)
before the 21-day scratch purge, then delete the staged archives.

:::{tip}
For a large multi-year run, stage and process one year at a time to keep the scratch
footprint bounded.
:::

## Input files

Most tools share a common set of input files, stored as examples in
`$MSCPROGS/Input/`. Copy the relevant file to the working directory and edit
as needed.

### `extract.*` files

Controls which fields are extracted. The tool auto-detects the file type from
the input filename and opens the corresponding extract file by name — so the
file must be present in the working directory under exactly that name (do not
rename it). There is a separate template for each file type:

| Template | File type | Notes |
|---|---|---|
| `extract.archm` | `archm` — time-averaged archive | Most commonly used; averaging period set by `meanfq` in `blkdat.input` |
| `extract.archv` | `archv` — instantaneous snapshot archive | Written at `diagfq` frequency |
| `extract.archv_wav` | `archv_wav` — wave archive | |

`Input/` also contains `extract.daily`, `extract.weekly`, `extract.restart`,
and `extract.sec` — these are for older TOPAZ-format output files and are not
used in standard NERSC HYCOM-CICE runs.

The `extract.sec` file in `Input/` is an additional example for section tools;
like all other tools, `m2section` auto-selects `extract.archm` or `extract.archv`
based on file type — there is no special `extract.sec` lookup.

```
4                      # version number — do not change
pres     intf          # field used to locate layer interfaces — do not change
2 F F T                # flags: (unused) sphere rotate index-velocities — do not change
1 270 1 210            # section-plot index range (ia ib ja jb) — do not change
1 1   XX XX   X        # 3-D layer range header — do not change
saln    01 01   T      # 3-D field: extract all layers
temp    01 01   T
pres    01 01   T
utot    01 01   T
vtot    01 01   T
ssh     00 00   T      # 2-D field (layer indices are 00 00 for 2-D fields)
fice    00 00   T
hice    00 00   F      # F = skip this field
```

Each field line has four columns:

| Column | Description |
|--------|-------------|
| `fieldname` | Must match the field name in the `.b` file exactly. Run `cat *.b` to see what is available in your output. |
| `layer_start` | First isopycnal layer to extract. Use `00` for 2-D fields. |
| `layer_end` | Last isopycnal layer to extract. Use `00` for 2-D fields. `layer_start layer_end` is a literal range: `01 01` extracts only layer 1, `01 26` extracts layers 1–26. |
| `T` / `F` | `T` to extract this field, `F` to skip it. |

The layer range is respected by `m2nc` and `h2nc`. For `hyc2proj` the range is
**ignored** — it always reads all isopycnal layers and then interpolates them
onto the z-levels in `depthlevels.in`, so `01 01` is conventional placeholder.

:::{note}
Field names differ between file types — `saln` / `temp` in daily output, `salin`
/ `thknss` in archv. Always copy the extract file that matches your input file type.
:::

### `proj.in`

Defines the target projection for `hyc2proj`. Example files for all four
projections are in `Input/`.

**Regular lon/lat** (`proj.in.regular_grid`):
```
regular    # projection type
-38.0      # longitude of SW corner (degrees E)
 20.0      # longitude of NE corner (degrees E)
0.125      # longitude grid spacing (degrees)
 65.0      # latitude of SW corner (degrees N)
 88.0      # latitude of NE corner (degrees N)
0.125      # latitude grid spacing (degrees)
```

**Mercator** (`proj.in.mercator`):
```
mercator
-100.0     # longitude W edge (degrees E)
  20.0     # longitude E edge (degrees E)
   800     # number of grid points in E-W direction
 -20.0     # latitude S edge (degrees N)
  70.0     # latitude N edge (degrees N)
   400     # number of grid points in S-N direction
    0.     # central longitude
```

**Polar stereographic** (`proj.in.stereographic`):
```
polar_stereographic
-3800.     # left  edge of projection area (km from pole)
 3800.     # right edge of projection area (km)
 12.5      # grid spacing in x-direction (km)
-5500.     # bottom edge of projection area (km)
 5500.     # top    edge of projection area (km)
 12.5      # grid spacing in y-direction (km)
-45.0      # central longitude of projection
 90.0      # central latitude  of projection
 0.01      # conversion factor (km → m if 0.001, km → km if 1.0)
 F         # rotate vectors to output grid (F = keep East/North components)
```

**Model native** (`proj.in.native`):
```
native
    1      # first i-index
  100      # last  i-index
    1      # i-index stride
    1      # first j-index
  100      # last  j-index
    1      # j-index stride
```
Extracts a sub-domain without horizontal interpolation.

### `depthlevels.in`

Specifies the target vertical depth levels (in metres) for `hyc2proj` and
`hyc2stations`. The first line is the number of levels; one depth per line
follows.

```
40        # number of depth levels
0.0
2.0
4.0
...
4000.0
```

### `sections.in`

Defines named sections as lon/lat endpoint pairs for `m2section`,
`m2transport`, and `m2transport2`.

```
1.2 F          # file version; F = quiet output (T = verbose)
#####
FRAM_STRAIT    # section name (used in output filenames)
-20. 79.       # start point: lon, lat (degrees)
15.  79.       # end  point:  lon, lat (degrees)
#
DENMARK_STRAIT
-37.  66.1
-22.68 66.
```

Sections follow great circles between the two endpoints. The transport sign
convention is: positive to the right when standing at the start point and
facing the end point.

### `transport.in`

Used by `m2transport2` to filter transport by depth range, salinity range,
or temperature range. Each line defines one transport component:

```
SVINOY-ATL-N   SVINOY   SALINITY   n   35.00   39.00
SVINOY         SVINOY   DEPTH      n    0.00 1000.00
```

| Column | Description |
|--------|-------------|
| Label | Output name used in result files |
| Section | Section name — must match an entry in `sections.in` |
| Criterion | `DEPTH`, `SALINITY`, or `TEMPERATURE` |
| Direction | `n` (no direction filter) |
| Min | Lower bound of the criterion range |
| Max | Upper bound of the criterion range |

### `stations.in`

Defines groups of station positions for `hyc2stations`. Each group starts
with `#GROUPNAME` and is followed by lon/lat pairs, one per line.

```
#FRAM_STRAIT
-20.00  79.00
-10.00  79.25
  0.00  79.50
 15.00  79.00
#DENMARK_STRAIT
-37.00  66.10
...
```

The helper script `setupstations.py` generates these point lists by
interpolating at regular spacing along a rhumb line between two endpoints.
