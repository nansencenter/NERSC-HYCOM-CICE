# MSCPROGS post-processing tools

MSCPROGS is a collection of Fortran (and C) programs for post-processing and
analysing HYCOM model output. It is a general-purpose toolkit covering grid projection, field extraction, temporal
averaging, section transports, Lagrangian ice drift, grid utilities, and more.
The tools live in `hycom/MSCPROGS/`.

For compilation instructions, see [Compile MSCPROGS](compilation.md#compile-mscprogs-libhycnersca).

## Running the tools

### Environment setup

Before running any MSCPROGS tool, load the HPC modules and set `$MSCPROGS`.
See the [HPC environment](installation.md#hpc-environment) section for details.

::::{dropdown} Source HPC environment — Betzy (NRIS/Sigma2)

```{include} _snippets/betzy_hpc_env.md
```

::::

```bash
source $WORK/<CONFIGNAME>/REGION.src    # sets $MSCPROGS and other paths
```

`REGION.src` also adds `$MSCPROGS/bin` and `$MSCPROGS/bin_setup` to your `PATH`, so all MSCPROGS executables are immediately available.

### Working directory

The natural place to run MSCPROGS tools is the experiment data directory:

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>/data/
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

`hyc2proj` requires three input files in the working directory, further described in  [Input files](#input-files) below:

- `proj.in` — target projection and grid
- `depthlevels.in` — vertical depth levels to interpolate to
- an `extract.*` file — fields to extract (e.g. `extract.daily`, `extract.archv`)

Horizontal interpolation is bilinear; vertical interpolation defaults to cubic
spline. `staircase` and `linear` options are also available (faster but lower
quality).

```bash
cp $MSCPROGS/Input/proj.in.regular_grid proj.in     # choose and edit a sample proj.in
cp $MSCPROGS/Input/depthlevels.in .
cp $MSCPROGS/Input/extract.daily .
hyc2proj TP2daily*.a                                 # process one or more .a files
```

Output is written to NetCDF files in the working directory.

`hyc2stations` is a companion routine that interpolates HYCOM fields onto
specified spatial positions (station locations or cross-section points) rather
than a grid. It reads `stations.in` instead of `proj.in`, together with
`depthlevels.in` and an `extract.*` file. Output is one NetCDF file per station
group. A helper script `setupstations.py` generates `stations.in` entries from
start/end coordinates along a rhumb line.

## `m2nc` / `m2t` — 2D field extraction to NetCDF

Reads HYCOM `.ab` files and writes selected 2D fields to `tmp1.nc`. Fields to
extract are specified in the relevant `extract.*` file. Multiple input files
produce multiple time records in the output. `m2t` additionally writes a Tecplot
file.

```bash
cp $MSCPROGS/Input/extract.daily .
m2nc TP2daily*.a
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

## `hycave` — temporal averaging

Computes a thickness-weighted time mean across the supplied files. `filetype` is
one of `restart`, `nersc_daily`, `nersc_weekly`, or `archv`. Requires
`regional.grid` and `regional.depth` in the working directory. Output files use
the placeholder date `9999_999`.

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
cp $MSCPROGS/Input/extract.daily extract.sec   # m2section uses extract.sec
m2section TP2daily*.a                          # output: section001.nc, section002.nc, ...
m2transport TP2daily*.a                        # output: transport_net001.dat, ...
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
against IABP buoy data or as input for data assimilation.

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

## Input files

Most tools share a common set of input files, stored as examples in `Input/`.

### `extract.*` files

Control which fields are extracted. There are separate files for each HYCOM file
type: `extract.restart`, `extract.archv`, `extract.daily`, `extract.weekly`.
From line 6 onward, each line specifies a field name (matching the `.b` file
header), a layer range, and a flag (`T`/`F`) to include or skip the field.
Field names can be looked up directly in the `.b` file of any HYCOM output.

### `proj.in`

Defines the target projection for `hyc2proj` and `RelaxToNetCDF`. Example files
for all four projections are in `Input/`.

### `depthlevels.in`

Specifies the vertical depth levels for `hyc2proj` and `hyc2stations`.

### `sections.in`

Defines named sections as pairs of lon/lat endpoints for the `m2section` /
`m2transport` / `m2transport2` scripts.

### `transport.in` and `scalartransport.in`

Used by `m2transport2` to define transport criteria (depth range, temperature
range, etc.) and optional scalar transport variables.

### `stations.in`

Defines station groups for `hyc2stations`. The helper script
`setupstations.py` can generate entries along a rhumb line between two points.
