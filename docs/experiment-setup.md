In the following, `<CONFIGNAME>`, `<EXPT_ID>`, and `<IEXPT>` are placeholders for user-defined
values, see the table in the [directory structure section](overview.md#directory-structure) for a
description and examples of each.

## Set up the work directory

Copy the model configuration template and create a symlink to the utility scripts:

```bash
CONFIGNAME=<CONFIGNAME>   # e.g. TP2a0.10

cd $WORK
cp -r $HOME/NERSC-HYCOM-CICE/${CONFIGNAME} .
cd $WORK/${CONFIGNAME}
ln -sf $HOME/NERSC-HYCOM-CICE/bin .
```

:::{note}
The `bin` symlink is required because scripts in `bin/` call each other using relative
`../bin/` paths, assuming they are run from an experiment subdirectory. Without the
symlink, those internal calls would fail.
:::

`expt_preprocess.sh` looks up `nest/`, `relax/`, and `force/` relative to the `<CONFIGNAME>/`
directory — redirect them with symlinks onto the scratch filesystem. Set `WDIR` to the scratch
path for your machine:

::::{dropdown} WDIR — scratch path by machine
::::{tab-set}
:::{tab-item} Betzy
```bash
WDIR=$USERWORK/${CONFIGNAME}
```
:::
:::{tab-item} Olivia
```bash
WDIR=/cluster/work/projects/nn2993k/$USER/${CONFIGNAME}
```
:::
::::
::::

Then create the directories and symlinks:

```bash
mkdir -p $WDIR $WDIR/nest

cd $WORK/${CONFIGNAME}
mv relax $WDIR/relax
mv force $WDIR/force
mkdir -p $WDIR/force/synoptic

ln -sf $WDIR/nest  nest
ln -sf $WDIR/relax relax
ln -sf $WDIR/force force
```

:::{important}
`nest/`, `relax/`, and `force/` hold input data that is expensive to regenerate. Because the
scratch filesystem is purged after 21 days, stage this data back from NIRD (or re-create it)
before restarting a run if these directories have aged out.
:::


## Configure REGION.src

Copy the template into the configuration directory:

```bash
CONFIGNAME=<CONFIGNAME>   # e.g. TP2a0.10

cd $WORK/${CONFIGNAME}
cp $HOME/NERSC-HYCOM-CICE/input/REGION.src .
```

Open `REGION.src` and update the following lines to match your setup:

```bash
export R=<CONFIGNAME>
...
export NHCROOT=${HOME}/NERSC-HYCOM-CICE
```

## Create the experiment directory

Use `expt_new.sh` to create a new experiment from the template experiment `01.0`:

```bash
CONFIGNAME=<CONFIGNAME>   # e.g. TP2a0.10
EXPT_ID=<EXPT_ID>         # e.g. 01.0

cd $WORK/${CONFIGNAME}
bin/expt_new.sh 01.0 ${EXPT_ID}
```

This creates `expt_<EXPT_ID>/` with default configuration files that you will
adjust in the following steps.

:::{note}
**When you must create a new experiment** — small parameter tweaks (e.g. viscosity)
can be made in an existing experiment. But the following changes require a new
experiment directory to avoid inconsistencies between the restart files, topography,
and relaxation fields:

- **Topography version change** — the bottom pressure field (`pbot`) stored in restart
  files is tied to the topography; mismatched versions will cause inconsistencies.
- **Vertical grid change** — any of the vertical discretisation parameters in
  `blkdat.input` (`kdm`, `nhybrd`, `nsigma`, `dp00*`, `ds00*`, sigma layer densities)
  must be consistent with the relaxation and climatology files.
- **HYCOM version change** — new model versions often introduce changes to
  `blkdat.input`; a fresh experiment directory lets you update the configuration
  without overwriting a working one.
:::

::::{dropdown} What does expt_new.sh do?

1. **Validates the environment** — the script must be run from either an experiment directory
   (containing `EXPT.src`) or the configuration root (containing `REGION.src`). It also
   checks that `topo/regional.grid.b` exists, which means the grid must already be in place.

2. **Copies the template experiment** — all files from `expt_<old>/` are copied to the new
   directory with `rsync`, skipping the `data/`, `log/`, and `SCRATCH/` subdirectories
   (those are created fresh and empty). This includes configuration files such as
   `blkdat.input`, `EXPT.src`, and the job submission script `srjob.sh`.

3. **Updates `EXPT.src`** — the `X=` and `E=` lines are rewritten to reflect the new
   experiment ID (e.g. `X="02.6"` and `E="026"`).

4. **Updates `blkdat.input`** — the `iexpt`, `idm`, and `jdm` fields are updated
   automatically: `iexpt` is set from the new experiment number, and `idm`/`jdm` are read
   directly from `topo/regional.grid.b`, so they always match the actual grid.

5. **Renames PBS job scripts** — any `pbsjob*.sh` files have their `-N` job-name line
   updated to `<CONFIGNAME>_X<EXPT_ID>` for easier identification in the job queue.

::::

Copy `hycom_opt` into the experiment directory:

```bash
CONFIGNAME=<CONFIGNAME>   # e.g. TP2a0.10
EXPT_ID=<EXPT_ID>         # e.g. 01.0

cd $WORK/${CONFIGNAME}/expt_${EXPT_ID}
cp $HOME/NERSC-HYCOM-CICE/TP5a0.06/expt_01.0/hycom_opt .
```

`hycom_opt` is a Fortran namelist read by the model at runtime. Open it and update the
values for your experiment. For TP2, the file should look like this:

```
&hycom_nml
  write_arche = .false.
  sss_underice= .false.
  highfq_river= .false.
  sssrmx_scalar=.5
/
```

| Option | Description |
|--------|-------------|
| `write_arche` | Write ESMF archive files at each coupling step. Useful for debugging the ESMF coupling interface; leave `.false.` for normal runs. |
| `sss_underice` | Apply SSS relaxation under sea ice. When `.false.`, SSS relaxation is suppressed where ice cover ≥ 15%. |
| `highfq_river` | Use time-varying (high-frequency) river forcing instead of the climatological river forcing. Requires `priver=0` in `blkdat.input`; the two options are mutually exclusive. |
| `sssrmx_scalar` | Maximum SSS anomaly (psu) at which relaxation is still applied. Relaxation is suppressed where the model–climatology difference exceeds this value. `99.` (template default) means no cap; `.5` limits relaxation to within 0.5 psu of climatology. |

## Configure blkdat.input

`blkdat.input` controls core model parameters. The easiest starting point is to copy it
from an existing experiment for your configuration. For TP2, reference files are available
at:

- **With nesting boundary**:
  `/nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/blkdat.input`
- **With climatology boundary**:
  `/nird/datalake/NS9481K/shuang/TP2_setup/exp02.6_seaclim_ref/blkdat.input_clim`

Copy it into place:

```bash
CONFIGNAME=<CONFIGNAME>   # e.g. TP2a0.10
EXPT_ID=<EXPT_ID>         # e.g. 01.0

cd $WORK/${CONFIGNAME}/expt_${EXPT_ID}
cp <path/to/reference/blkdat.input> .
```

Each field occupies its own line; change the leading numeric value only. At minimum,
update these fields for each new experiment:

| Field | Description | Notes |
|-------|-------------|-------|
| `iexpt` | Experiment number ×10 | e.g. `026` for expt `02.6` |
| `ntracr` | Number of BGC tracers | `0` = none, `1` = ECOSMO |

::::{dropdown} Full blkdat.input parameter reference (with example values)

**Grid and configuration**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `iversn` | `23` | HYCOM version number ×10 |
| `iexpt` | `026` | Experiment number ×10 |
| `idm` | `400` | Longitudinal array size |
| `jdm` | `380` | Latitudinal array size |
| `itest` | `-1` | Grid point for detailed diagnostics (-1=none) |
| `jtest` | `-1` | Grid point for detailed diagnostics (-1=none) |
| `kdm` | `50` | Number of layers |
| `nhybrd` | `50` | Number of layers in the hybrid generator (0=all isopycnal) |
| `nsigma` | `0` | Number of terrain-following sigma levels at the surface |

**Vertical grid spacing**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `dp00` | `1.2` | Deep z-level minimum thickness (m) |
| `dp00x` | `250.0` | Deep z-level maximum thickness (m) |
| `dp00f` | `1.1` | Deep z-level stretching factor (1.0=constant) |
| `ds00` | `1.2` | Shallow z-level minimum thickness (m) |
| `ds00x` | `75.0` | Shallow z-level maximum thickness (m) |
| `ds00f` | `1.1` | Shallow z-level stretching factor (1.0=constant) |
| `dp00i` | `1.2` | Deep isopycnal spacing minimum thickness (m) |
| `isotop` | `10.0` | Shallowest depth for isopycnal layers (m; <0=from file) |
| `oneta0` | `0.1` | Minimum 1+eta, must be >0 |

**Layer densities**

One isopycnal target density (`sigma`) per layer — `nhybrd` values in total. Each layer
tracks water of that target density when the water column is sufficiently stratified. The
choice of target densities determines the vertical coordinate structure and is
configuration-specific; take these values from a reference file for your configuration.

In the TP2 configuration, the first 10 entries are assigned unrealistically light densities
(sigma 0.1–1.0). Real ocean water never reaches those densities, so these layers remain
near the surface and act as quasi-fixed surface layers. Layers 11–50 carry realistic
isopycnal target densities ranging from 24.05 to 28.12 (sigma units).

**Hybrid generator**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `hybraf` | `0` | Robert-Asselin flag (0=F, 1=T) |
| `hybrlx` | `8.0` | Inverse relaxation coefficient (time steps) |
| `hybiso` | `0.01` | Use PCM if layer is within this value of target density |
| `hybthn` | `1.0` | Ratio of layer thicknesses to select the thinner |
| `hybthk` | `0.0` | Thick-thin-thick ratio to expand the thin layer |
| `hybmap` | `3` | Remapper flag (0=PCM, 1=PLM, 2=PPM) |
| `hybflg` | `0` | Generator flag (0=T&S, 1=th&S, 2=th&T) |

**Initialization**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `saln0` | `35.0` | Initial salinity (psu); only used when `iniflg` < 2 |
| `locsig` | `0` | Locally-referenced potential density for stability (0=F, 1=T) |
| `kapref` | `0` | Thermobaric reference state (-1=input, 0=none, 1–3=constant) |
| `thflag` | `0` | Reference pressure flag (0=Sigma-0, 2=Sigma-2) |
| `thbase` | `25.0` | Reference density (sigma units) |
| `vsigma` | `0` | Spatially varying isopycnal target densities (0=F, 1=T) |
| `iniflg` | `2` | Initial state flag (0=level, 1=zonal, 2=climatology) |
| `jerlv0` | `0` | Initial Jerlov water type (1–5; 0=use KPAR) |

**Time stepping**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `yrflag` | `3` | Days-in-year flag (0=360, 1=366, 2=366J1, 3=actual) |
| `baclin` | `600.0` | Baroclinic time step (s); must divide 86400 |
| `batrop` | `30.0` | Barotropic time step (s); must divide baclin/2 |
| `ra2fac` | `0.125` | Robert-Asselin time filter weight |
| `wbaro` | `0.125` | Barotropic time smoothing weight |
| `btrlfr` | `1` | Leapfrog barotropic time step (0=F, 1=T) |
| `btrmas` | `0` | Barotropic mass conservation (0=F, 1=T) |

**Advection and dynamics**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `advflg` | `0` | Thermal advection flag (0=T&S, 1=th&S, 2=th&T) |
| `advtyp` | `2` | Scalar advection type (0=PCM, 1=MPDATA, 2=FCT2, 4=FCT4) |
| `momtyp` | `2` | Momentum advection type (2=2nd order, 4=4th order) |
| `slip` | `-1.0` | Boundary condition (+1=free-slip, -1=non-slip) |
| `visco2` | `0.05` | Deformation-dependent Laplacian viscosity factor |
| `visco4` | `0.2` | Deformation-dependent biharmonic viscosity factor |
| `facdf4` | `0.0` | Speed-dependent biharmonic viscosity factor |

**Diffusion**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `veldf2` | `0.005` | Laplacian momentum diffusion velocity (m/s) |
| `veldf4` | `0.05` | Biharmonic momentum diffusion velocity (m/s; <0=from file) |
| `thkdf2` | `0.0` | Laplacian thickness diffusion velocity (m/s) |
| `thkdf4` | `0.05` | Biharmonic thickness diffusion velocity (m/s; <0=from file) |
| `temdf2` | `0.005` | Laplacian temperature/salinity diffusion velocity (m/s) |
| `temdfc` | `1.0` | Temperature diffusion conservation (0=density, 1=temperature) |
| `vertmx` | `2.e-5` | Momentum diffusion velocity at MICOM mixed layer base (m/s) |

**Bottom friction**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `cbar` | `0.1` | RMS flow speed for linear bottom friction (m/s) |
| `cb` | `2.e-3` | Quadratic bottom friction coefficient |
| `drglim` | `0.0` | Limiter for explicit friction (1.0=none, 0.0=implicit) |
| `thkbot` | `10.0` | Bottom boundary layer thickness (m) |

**Mixed layer and interface**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `sigjmp` | `0.01` | Minimum density jump across interfaces (kg/m³) |
| `tmljmp` | `0.2` | Equivalent temperature jump across mixed layer (°C) |
| `sefold` | `30.0` | SSS relaxation e-folding time |
| `thkmls` | `15.0` | Reference mixed-layer thickness for SSS relaxation (m) |
| `thkmlt` | `0.0` | Reference mixed-layer thickness for SST relaxation (m) |
| `thkriv` | `6.0` | Nominal thickness of river inflow (m) |
| `thkcdw` | `3.0` | Thickness for near-surface currents in ice-ocean stress (m) |
| `thkfrz` | `10.0` | Maximum thickness of near-surface freezing zone (m) |

**KPP mixed layer**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `mlflag` | `1` | Mixed layer scheme (0=none, 1=KPP, 2–3=KT, 4=PWP, 5=MY, 6=GISS) |
| `pensol` | `1` | Penetrating solar radiation (0=F, 1=T) |
| `thkmin` | `19.2` | Minimum mixed-layer thickness (m) |
| `dypflg` | `1` | Diapycnal mixing flag (0=none, 1=KPP, 2=explicit) |
| `mixfrq` | `64` | Time steps between diapycnal mixing calculations |
| `diapyc` | `1.e-7` | Diapycnal diffusivity × buoyancy frequency (m²/s²) |
| `rinfty` | `0.7` | KPP: maximum gradient Richardson number (shear instability) |
| `ricr` | `0.3` | KPP: critical bulk Richardson number |
| `bldmin` | `0.0` | KPP: minimum surface boundary layer thickness (m) |
| `bldmax` | `1200.0` | KPP: maximum surface boundary layer thickness (m) |
| `cekman` | `0.7` | KPP/KT: Ekman depth scale factor |
| `cmonob` | `1.0` | KPP: Monin-Obukhov depth scale factor |
| `bblkpp` | `1` | KPP: activate bottom boundary layer (0=F, 1=T) |
| `shinst` | `1` | KPP: activate shear instability mixing (0=F, 1=T) |
| `dbdiff` | `1` | KPP: activate double diffusion mixing (0=F, 1=T) |
| `nonloc` | `1` | KPP: activate nonlocal boundary layer mixing (0=F, 1=T) |
| `botdiw` | `0` | Activate bottom-enhanced internal wave mixing (0=F, 1=T) |
| `difout` | `1` | Output viscosity/diffusivity coefficients in archive (0=F, 1=T) |
| `difm0` | `50.e-4` | KPP: max viscosity due to shear instability (m²/s) |
| `difs0` | `50.e-4` | KPP: max diffusivity due to shear instability (m²/s) |
| `difmiw` | `3.e-5` | KPP: background/internal wave viscosity (m²/s) |
| `difsiw` | `1.e-5` | KPP: background/internal wave diffusivity (m²/s) |
| `dsfmax` | `10.e-4` | KPP: salt fingering diffusivity factor (m²/s) |
| `hblflg` | `2` | KPP: boundary layer interpolation flag (0=const., 1=linear, 2=quad.) |
| `niter` | `2` | KPP: iterations for semi-implicit solution (2 recommended) |
| `langmr` | `0` | KPP: activate Langmuir turbulence factor (0=F, 1=T) |

**Sea ice**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `iceflg` | `2` | Sea ice model (0=none, 1=energy loan, 2=coupled/ESMF) |
| `tfrz_0` | `0.0` | Ice melting point at S=0 psu (°C) |
| `tfrz_s` | `-0.054` | Gradient of ice melting point (°C/psu) |
| `frzifq` | `0.25` | E-folding time back to freezing point (days; <0=time steps) |
| `ticegr` | `2.0` | Temperature gradient inside ice (°C/m); 0=use surface temp |
| `hicemn` | `0.1` | Minimum ice thickness (m) |
| `hicemx` | `10.0` | Maximum ice thickness (m) |
| `ishelf` | `0` | Ice shelf flag (0=none, 1=ice shelf over ocean) |

**Tracers**

| Field | Example | Description |
|-------|---------|-------------|
| `ntracr` | `1` | Number of tracers (0=none; negative=initialize) |
| `trcflg` | `9` | Tracer flags (one digit per tracer, most significant replicated) |
| `mtracr` | `0` | Number of diagnostic tracers |
| `tsofrq` | `64` | Time steps between anti-drift offset calculations |
| `tofset` | `0.0` | Temperature anti-drift offset (°C/century) |
| `sofset` | `0.0` | Salinity anti-drift offset (psu/century) |

**Forcing**

| Field | Example | Description |
|-------|---------|-------------|
| `clmflg` | `12` | Climatology frequency (6=bimonthly, 12=monthly) |
| `wndflg` | `4` | Wind stress input (0=none, 1=u/v-grid, 2–3=p-grid, 4–5=10m wind) |
| `ocnscl` | `1.0` | Scale factor for ocean current in relative wind (0=absolute wind) |
| `ustflg` | `2` | ustar forcing flag (1–2=wind speed, 3=input, 4=stress) |
| `flxflg` | `6` | Thermal forcing flag (0=none, 3=net flux, 1–2/4–6=SST-based) |
| `empflg` | `6` | E-P forcing flag (0=none, 3=net E-P, 1–2/4–6=SST-based) |
| `emptgt` | `0.0` | E-P balance target (mm/week into ocean) |
| `empbal` | `1` | E-P balance flag (0=none, 1=offset, 2=scale) |
| `dswflg` | `0` | Diurnal shortwave flag (0=none, 1=daily-to-diurnal correction) |
| `albflg` | `2` | Ocean albedo flag (0=none, 1=constant, 2=L&Y) |
| `sssflg` | `1` | SSS relaxation flag (0=none, 1=climatology) |
| `sssbal` | `0` | SSS relaxation balance flag (0=none, 1=offset, 2=scale) |
| `lwflag` | `-1` | Longwave flag (0=none, 1=climatology, 2=atmosphere) |
| `sstflg` | `0` | SST relaxation flag (0=none, 1=clim, 2=atmos, 3=observed) |
| `icmflg` | `0` | Ice mask flag (0=none, 1=clim, 2=atmos, 3=obs/coupled) |
| `prsbas` | `0` | MSL pressure offset: input field + prsbas (Pa) |
| `mslprf` | `1` | MSL pressure forcing flag (0=F, 1=T) |
| `stroff` | `0` | Net stress offset flag (0=F, 1=T) |
| `flxoff` | `0` | Net flux offset flag (0=F, 1=T) |
| `flxsmo` | `0` | Smooth surface fluxes (0=F, 1=T) |

**Lateral boundary and nesting**

| Field | Example | Description |
|-------|---------|-------------|
| `lbflag` | `0` | Lateral barotropic boundary flag (0=none, 1=port, 2=input) |
| `lbmont` | `0` | Barotropic nesting archives have sshflg=2 |
| `relax` | `1` | Activate lateral boundary nudging (0=F, 1=T) |
| `trcrlx` | `1` | Activate lateral boundary tracer nudging (0=F, 1=T) |
| `priver` | `1` | Rivers as precipitation (0=F, 1=T) |
| `epmass` | `0` | Treat E-P as mass exchange (0=F, 1=T) |

**Incremental update**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `incflg` | `0` | Incremental update flag (0=no, 1=yes, 2=full-velocity) |
| `incstp` | `100` | Number of time steps for full update (1=direct insertion) |
| `incupf` | `1` | Number of days of incremental updating input |

**Stochastic forcing**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `stfflg` | `0` | Stochastic anomaly forcing (0=no, 1=TS, 2=TSV) |
| `stfrdt` | `50.0` | Stochastic T anomaly e-folding depth (m) |
| `stfrds` | `50.0` | Stochastic S anomaly e-folding depth (m) |
| `stfrdv` | `50.0` | Stochastic V anomaly e-folding depth (m) |

**Tides**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `tidflg` | `0` | Tidal forcing flag (0=none, 1=open boundary, 2=boundary+body) |
| `tidein` | `0` | Tide field input flag (0=no, 1=yes, 2=SAL) |
| `tidcon` | `11111111` | Tidal constituents: one digit per (Q1 K2 P1 N2 O1 K1 S2 M2) |
| `tidsal` | `0.06` | Scalar self-attraction and loading factor (<0=from file) |
| `tiddrg` | `0` | Tidal drag flag (0=no, 1=scalar, 2=tensor) |
| `thkdrg` | `500.0` | Bottom boundary layer thickness for tidal drag (m) |
| `drgscl` | `0.0` | Scale factor for tidal drag (0=no tidal drag) |
| `tidgen` | `1` | Generic time (0=F, 1=T) |
| `tidrmp` | `3.0` | Tidal ramp time (days) |
| `tid_t0` | `0.0` | Origin for tidal ramp time (model day) |

**Synthetic floats**

| Field | TP2 value | Description |
|-------|-----------|-------------|
| `fltflg` | `0` | Synthetic float flag (0=no, 1=yes) |
| `nfladv` | `4` | Advect every nfladv baroclinic time steps (even, ≥4) |
| `nflsam` | `1` | Output (0=every nfladv steps; >0=number of days) |
| `intpfl` | `0` | Horizontal interpolation (0=2nd order+n.n., 1=n.n. only) |
| `iturbv` | `0` | Add horizontal turbulent advection velocity (0=no, 1=yes) |
| `ismpfl` | `1` | Sample water properties at float (0=no, 1=yes) |
| `tbvar` | `4.63e-6` | Horizontal turbulent velocity variance scale (m²/s²) |
| `tdecri` | `0.4` | Inverse decorrelation time scale (1/day) |

**Output frequencies**

| Field | Example | Description |
|-------|---------|-------------|
| `sshflg` | `0` | SSH output flag (0=SSH, 1=SSH+steric SSH) |
| `dsurfq` | `0.0` | Days between surface diagnostics |
| `diagfq` | `0.0` | Days between full diagnostics |
| `proffq` | `0.0` | Days between profile diagnostics |
| `tilefq` | `0.0` | Days between tile diagnostics |
| `meanfq` | `1.0` | Days between time-averaged diagnostics |
| `rstrfq` | `5.0` | Days between restart output |
| `bnstfq` | `0.0` | Days between barotropic nesting archive input |
| `nestfq` | `0.0` | Days between 3-D nesting archive input |
| `cplifq` | `-6` | Days (or time steps) between sea ice coupling |

**Output file names**

| Field | Example | Description |
|-------|---------|-------------|
| `nmrsti` | `restart` | Restart input filename (dot and date appended; empty=old HYCOM name) |
| `nmrsto` | `restart` | Restart output filename (dot and date appended; empty=old HYCOM name) |
| `nmarcs` | `archv_s` | Surface archive filename |
| `nmarcv` | `archv` | Full archive filename |
| `nmarcm` | `archm` | Mean archive filename |

::::



## Configure EXPT.src

`EXPT.src` sets experiment-level parameters sourced at runtime and during compilation.
Open `$WORK/<CONFIGNAME>/expt_<EXPT_ID>/EXPT.src` and update:

| Line to find | Value to set | Notes |
|--------------|--------------|-------|
| `X=` | `"<EXPT_ID>"` e.g. `"02.6"` | Experiment identifier (dot notation) |
| `E=` | `"<IEXPT>"` e.g. `"026"` | Experiment identifier (no dot) |
| `T=` | `"04"` | Topography version |
| `export NMPI=` | e.g. `504` for TP2 on Betzy | Number of ocean MPI tiles |
| `export MXBLCKS=` | e.g. `9` | Maximum ice blocks per MPI process |
| `export COMPILE_BIOMODEL=` | `"yes"` or `"no"` | BGC coupling on/off |

:::{note}
For TP2, the available topography versions are: `01` (initial interpolation), `02` (adds
Ob river channel), `03` (identical to `02`), `04` (blends `02`/`03` with NEMO topography
at the nesting boundary). Higher numbers are newer refinements, not changes in resolution.
Use `04` for current experiments.
:::

:::{note}
`NMPI` is determined by the tiling step in
[Preparing External Files](external-files.md) (which requires
[hycom_ALL](compilation.md#compile-hycom-all) to be compiled first). For TP2 on Betzy,
with topography version `04` and `Icore=29, Jcore=26`, the result is `504`, so you can
set that now. For other configurations, leave it as a placeholder and fill it in after
running `create_ref_case.sh`. Compilation of HYCOM-CICE does not use `NMPI`.
:::

:::{note}
If the model reports that ice blocks exceed the maximum, increase `MXBLCKS` to the value
recommended in the error message.
:::

::::{dropdown} Other variables in EXPT.src

| Variable | Description |
|----------|-------------|
| `export V=` | HYCOM version; determines which source directory is used when compiling |
| `export SIGVER=` | Equation of state version; must be consistent with `thflag` in `blkdat.input` |
| `export K=` | Number of layers — auto-derived from `blkdat.input`, no need to edit |
| `export P=` | Experiment directory path — set automatically from the script location |
| `export D=` | Permanent data directory (`P/data`) — set automatically |
| `export S=` | Scratch directory (`P/SCRATCH`) — set automatically |

::::

::::{dropdown} Both machines — redirect `SCRATCH` and `data` onto the scratch filesystem

Keep the experiment tree (configuration and `build/`) on the non-purged `$WORK`
(`/cluster/projects/nn2993k/$USER`), and put the two large directories — scratch and output —
on the fast, purged scratch filesystem. Override the auto-set `S=` and `D=` lines in `EXPT.src`.
The scratch path differs by machine:

::::{tab-set}
:::{tab-item} Betzy
```bash
export S=$USERWORK/<CONFIGNAME>/expt_<EXPT_ID>/SCRATCH
export D=$USERWORK/<CONFIGNAME>/expt_<EXPT_ID>/data
```
:::
:::{tab-item} Olivia
```bash
export S=/cluster/work/projects/nn2993k/$USER/<CONFIGNAME>/expt_<EXPT_ID>/SCRATCH
export D=/cluster/work/projects/nn2993k/$USER/<CONFIGNAME>/expt_<EXPT_ID>/data
```
:::
::::

This is safe: `expt_preprocess.sh` only creates (`mkdir -p`) and enters (`cd`) `$S` and `$D` —
it never deletes them — and the overrides propagate automatically when `expt_new.sh` copies
`EXPT.src` to a new experiment. `P=` and `build/` stay on `$WORK`, so the configuration and
compiled executable survive the purge.

:::{note}
`expt_preprocess.sh` also accesses `relax/` via `${D}/../../relax/` (two levels up from `$D`).
When `D=` is overridden to a path on the scratch filesystem, `${D}/../../` resolves to the
scratch `<CONFIGNAME>/` subtree — so `relax/` must be on scratch as well, which is exactly what
the symlinks described above ([Set up the work directory](#set-up-the-work-directory)) provide.
:::

:::{important}
`data/` is now on the purged filesystem, so **archive completed output to NIRD on a rolling
basis** (e.g. per model year as it finishes) from a service node — the scratch filesystem is
purged by file age, so early output can age out while a long run is still going.
:::

::::

