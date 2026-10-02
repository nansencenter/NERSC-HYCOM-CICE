(submit-a-job)=
## Submit a job

The job script `srjob.sh` handles three things automatically:

- **[Atmospheric forcing generation](forcing.md#atmospheric-forcing)** — downloads and
  prepares ERA5 forcing for the run period
- **Stage run files** — runs `expt_preprocess.sh` to stage all files the model needs into
  the scratch directory (`expt_<EXPT_ID>/SCRATCH/`). Specifically, it:
  - Reads and validates settings from `blkdat.input` (time steps, flags, experiment number)
  - Copies the MPI partition file (`topo/partit/depth_R_T.NNNN`, where `NNNN` is the
    MPI task count) to `patch.input` in the scratch directory
  - Symlinks or copies all forcing files (atmospheric, river, kpar, relaxation
    climatologies, diffusion fields, nesting files)
  - Verifies that atmospheric forcing covers the full run period
  - Copies restart files (HYCOM and CICE) from the data directory
  - Copies the `hycom_cice` executable from `build/`
  - Writes the `limits` file with start and stop times via `hycom_limits.py`
  - Moves old output files to `KEEP/`

  It exits with a non-zero code if any required file is missing, aborting the job before
  the model starts.

:::{warning}
As all the files needed are staged into the scratch directory (`expt_<EXPT_ID>/SCRATCH/`), where the model runs, `expt_preprocess.sh` should be run for every new build or change in the files mentioned above. Otherwise, the previous versions of the files remain staged. By default, this step is handled automatically by the job script `srjob.sh`.
:::

  ::::{dropdown} Run manually

  ```bash
  cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
  ../bin/expt_preprocess.sh ${START} ${END} $INITFLG
  ```

  Use the same `START`, `END`, and `INITFLG` values as in `srjob.sh` (see the table
  below for formats). A successful run prints `No fatal errors. Ok to start model set up in ..`.

  Running `expt_preprocess.sh` manually beforehand is useful to verify that all required
  files are in place before submitting the job. If you do, comment out the
  `expt_preprocess.sh` call in `srjob.sh` to avoid staging the files twice.
  ::::

- **Job submission** — submits the model to the SLURM queue

- **Post-processing** — runs `expt_postprocess.sh` after the model finishes. It moves
  restart files, daily mean archives (`archv.*`), and CICE output from the scratch
  directory to `data/` and `data/cice/`, and writes `log/hycom.stop` with `GOODRUN` or
  `BADRUN` depending on whether the model reached a normal stop. The `data/` directory
  is then ready for analysis with the [MSCPROGS post-processing tools](mscprogs.md).

Open `srjob.sh` in a text editor and update:

| Variable | Value | Description |
|----------|-------|-------------|
| `NMPI` | `"<NMPI>"` e.g. `504` | Update NMPI if needed |
| `START` | `"YYYY-MM-DDT00:00:00"` | Run start time |
| `END` | `"YYYY-MM-DDT00:00:00"` | Run end time |
| `INITFLG` | `""`, `"--init"` or `"--init-ice"` | `"--init"` for a cold start; `""` to restart from files; `"--init-ice"` for a HYCOM restart with a CICE cold start (see notes below) |
| `#SBATCH --time` | `"HH:MM:SS"` | Wall-clock time limit |

> **`INITFLG="--init"` (cold start):** No restart files are needed.
> Temperature and salinity (T/S) are set directly from the climatological fields in
> `relax/` (see [Climatologies and river forcing](forcing.md#climatologies-and-river-forcing));
> velocities and sea surface height (SSH) start at zero. The model then spins up under realistic atmospheric
> forcing. The start date must be in September (the month of Arctic sea ice minimum).
> See [Initial conditions](forcing.md#initial-conditions) for details on all options.

> **`INITFLG=""` (restart from files):** HYCOM and CICE read from restart files in
> `data/` at the start date. The open boundary forcing is determined by `blkdat.input`,
> not by this flag, so three distinct scenarios share this setting:
>
> - **Continuing a spin-up** (e.g. the previous job hit the wall-time limit): leave
>   `blkdat.input` unchanged. Climatological boundaries remain active. Update `START`
>   to the date of the last restart file in `data/` and resubmit.
> - **Starting a hindcast or forecast from a spun-up state**: update `blkdat.input` to
>   GLORYS settings (`relax=0`, `nestfq=1`, `bnstfq=1`, `lbflag=2`) and stage the
>   GLORYS nesting files before submitting. See
>   [Open boundary forcing](forcing.md#open-boundary-forcing).
> - **Continuing a hindcast or forecast** (e.g. the previous job hit the wall-time
>   limit): leave `blkdat.input` unchanged. GLORYS boundaries remain active. Update
>   `START` to the date of the last restart file in `data/` and resubmit.

> **`INITFLG="--init-ice"` (HYCOM restart, CICE cold start):** HYCOM reads its restart
> file from `data/` at the start date, as for `INITFLG=""`, while CICE is initialised
> from `ice_initial.nc` as for `INITFLG="--init"`; no CICE restart file is needed. Use this
> when the ocean initial state comes from a restart file that was not written together
> with a CICE restart, e.g. one built from a GLORYS state. The start date must be in September.
> See [Initial conditions](forcing.md#initial-conditions).

:::{note}
For reference when setting `#SBATCH --time`: a 1-year TP2 run with BGC on 4 Betzy nodes (504 cores) takes approximately 5–6 hours of wall time.
:::

Then submit:

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
sbatch srjob.sh
```

Check job status with:

```bash
squeue -u $USER
```

A useful tool for generating SLURM job scripts is available at:
<https://open.pages.sigma2.no/job-script-generator/>


## Check results

When the job finishes, restart files and daily mean files are moved to the `data/` directory.
Confirm successful completion by checking the stop file:

```bash
cat $WORK/<CONFIGNAME>/expt_<EXPT_ID>/log/hycom.stop
```

The file should contain `GOODRUN`. If not, inspect the log files under `log/` for
error messages.

## Restarting after a crash

What happens to output files depends on how the run ended:

- **Model exits with an error but the SLURM job completes normally** (e.g. a model
  error or failed assertion): `expt_postprocess.sh` still runs, so restart and
  archive files are moved to `data/` and `log/hycom.stop` contains `BADRUN`.
  The restart files written before the crash are already in `data/` and can be
  used directly.

- **Job is killed by SLURM** (wall-time limit, out-of-memory, node failure, or
  `scancel`): `expt_postprocess.sh` never runs. Output files remain in `SCRATCH/`.
  Run postprocessing manually before resubmitting:

  ```bash
  cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
  ../expt_postprocess.sh
  ```

  This moves the restart and archive files written before the crash to `data/`,
  where `srjob.sh` will pick them up on the next run. If you skip this step and
  resubmit directly, `expt_preprocess.sh` will move the files to `SCRATCH/KEEP/`
  instead — they are not lost, but you will need to copy them to `data/` manually.

In both cases, update `START` in `srjob.sh` to the date of the last restart file
written, and set `INITFLG=""` (restart run).

## Cycled spin-up

:::{warning}
The cycled spin-up is **experimental** and not yet tested in production runs. The
standard spin-up with climatological boundaries (see
[Initial conditions](forcing.md#initial-conditions)) remains the default. Nothing
described here changes the behaviour of `srjob.sh` or existing experiments.
:::

As an alternative to spinning up with climatological boundaries, the model can be run
repeatedly over a fixed period with GLORYS boundaries (e.g. 1993–1997 five times, or
1993–2002 twice). `srjob_cycle.sh` does this; it is in the TP2 (`TP2a0.10/expt_01.0`) and
TP5 (`TP5a0.06/expt_02.3`) template experiments, which differ only in the SLURM settings
and `NMPI`. `blkdat.input` keeps the GLORYS nesting settings throughout; no new
`blkdat.input` option is needed. The first cycle can start on 1 September 1993 from the
GLORYS state (see [Initial state for 1 September 1993](#initial-state-for-1-september-1993)).

The model always runs on real dates inside the period, so atmospheric forcing, nesting
files and river forcing are used unchanged. At the end of each cycle, `bin/cycle_wrap.py`
copies the HYCOM and CICE restart files valid at `CYCLE_END` to the next cycle as restart
files valid at `CYCLE_START`, shifting only their time information (`nstep`, `dtime` in the
HYCOM `.b` header; `istep1`, `time` and the calendar attributes of the CICE restart).

| Variable | Example | Description |
|----------|---------|-------------|
| `SPINUP_START` | `"1993-09-01T00:00:00"` | Start of the first cycle; may lie inside the period |
| `CYCLE_START` | `"1993-01-01T00:00:00"` | Start of cycles 2, 3, … |
| `CYCLE_END` | `"1998-01-01T00:00:00"` | End of every cycle (wrap date, exclusive) |
| `NCYCLES` | `5` | Number of cycles |
| `SEGMENT_MONTHS` | `12` | Length of one preprocess/run/postprocess segment |
| `SEGMENTS_PER_JOB` | `8` | Segments per job before it resubmits itself (`RESUBMIT="yes"`) |

- Each cycle writes to its own data directory `data/cycle_NN`, so archives and
  restart files of different cycles do not overwrite each other. `D` in `EXPT.src`
  must contain `${SPINUP_CYCLE:+/cycle_${SPINUP_CYCLE}}`, e.g.
  `export D=$P/data${SPINUP_CYCLE:+/cycle_${SPINUP_CYCLE}}`; `SPINUP_CYCLE` is only set by `srjob_cycle.sh`, so other job scripts still write to `data/`.
- The first cycle needs the HYCOM restart for `SPINUP_START` in `data/cycle_01`. If there is
  no CICE restart for `SPINUP_START`, CICE is cold-started from `ice_initial.nc` in the
  experiment directory (`INITFLG="--init-ice"`, see
  [Initial conditions](forcing.md#initial-conditions)). For 1 September 1993, see below.
- The nesting files must cover `CYCLE_END` itself, since HYCOM interpolates in time
  between daily files.
- Both models must write a restart at the end of every segment: `rstrfq` in `blkdat.input`
  must be positive (a negative value suppresses the end-of-run restart) and `dump_last = .true.`
  in `ice_in`. Restart dates need not line up with `CYCLE_END`. The job stops if a segment
  ends without both restart files.
- The job continues from the latest HYCOM/CICE restart pair in the cycle data
  directories, so it can simply be resubmitted after a crash. If the job was killed by
  SLURM, run postprocessing with the cycle of the interrupted segment first, e.g.
  `SPINUP_CYCLE=02 ../bin/expt_postprocess.sh`.
- The boundary and atmospheric forcing jump from `CYCLE_END` back to `CYCLE_START` at each
  wrap. Atmospheric CO2 (`co2_annmean_gl.txt`) follows the model date and is also reset
  each cycle, which matters for BGC.

### Initial state for 1 September 1993

The ocean starts from the GLORYS12 state of the nesting file for 1993-09-01, the sea ice
from GLORYS12 concentration and TOPAZ4b thickness. The files are shared on NIRD:

```
/nird/datapeak/NS9481K/SPINUP_INIT/
├── README.md
├── sources/
│   ├── glorys12_ice_19930901.nc    # GLORYS12 siconc, sithick (Copernicus Marine)
│   └── topaz4b_ice_19930901.nc     # TOPAZ4b siconc, sithick, sisnthick (Copernicus Marine)
└── TP2a0.10/                       # topography 04, 50 layers, HYCOM 2.3
    ├── ice_initial_19930901.nc
    └── restart.1993_244_00_0000.[ab]
```

For TP2a0.10 with topography 04, 50 layers and HYCOM 2.3, copy the ready-made files:

```bash
EXPT_ID=<EXPT_ID>         # e.g. 03.0
N=/nird/datapeak/NS9481K/SPINUP_INIT/TP2a0.10
mkdir -p $WDIR/expt_${EXPT_ID}/data/cycle_01
cp $N/restart.1993_244_00_0000.[ab] $WDIR/expt_${EXPT_ID}/data/cycle_01/
cp $N/ice_initial_19930901.nc $WORK/TP2a0.10/expt_${EXPT_ID}/ice_initial.nc
```

where `WDIR` is the scratch path (see [Restart files](forcing.md#restart-files)); the
restart must end up in the experiment's `D` for cycle 01.

For another configuration (e.g. TP5a0.06, or another topography or number of layers),
build the two files from the experiment directory; this needs the nesting file for
1993-09-01 and a restart file of the same configuration as template:

```bash
CONFIGNAME=<CONFIGNAME>   # e.g. TP5a0.06
IEXPT=<IEXPT>             # nesting dir, e.g. 023
T=<T>                     # topography version, e.g. 05
N=/nird/datapeak/NS9481K/SPINUP_INIT

# HYCOM restart restart.1993_244_00_0000.[ab]
python $HOME/NERSC-HYCOM-CICE/bin/archv2restart.py ../nest/${IEXPT}/archv.1993_244_00 \
    <template restart without .a/.b> --outdir <data dir of cycle 01>

# CICE initial state
python $HOME/NERSC-HYCOM-CICE/bin/cice_ice_initial.py 1993-09-01T00:00:00 ../nest/${IEXPT}/archv.1993_244_00 \
    --griddir ../topo --kmt ../topo/kmt_${CONFIGNAME}_${T}.nc \
    --aice-source glorys --hi-source topaz \
    --glorys-file $N/sources/glorys12_ice_19930901.nc --topaz-file $N/sources/topaz4b_ice_19930901.nc \
    -o ice_initial.nc
```

- `archv2restart.py` takes `pbot`, `psikk` and `thkk` from the template restart. `montg1` in
  the nesting files must be consistent with its `psikk`/`thkk` (check with `calc_montg1.py`,
  see [Nesting files](forcing.md#nesting-files)), otherwise the barotropic pressure is wrong.
  NaN values in the nesting file are filled with the nearest ocean value. Options
  `--zero-velocity` and `--zero-ssh` start the ocean with zero velocities or SSH instead.
- `cice_ice_initial.py` can also use PIOMAS (`--aice-source piomas`, `--hi-source piomas`)
  or other combinations of the three sources. Interpolation is done on points, so it does
  not depend on the structure of the source or target grid. The GLORYS land mask is read
  from the GLORYS physics file of the same day in `MERCATOR_DATA/PHY`.
- Please add files built for a new configuration to `SPINUP_INIT/<CONFIGNAME>/` and
  describe them in its README.

## Visualization

If you prefer to work with xarray for analysis and visualization, check out [xhycom](https://xhycom.readthedocs.io/en/latest/index.html).

Otherwise, sample Jupyter notebooks for plotting HYCOM output are also available in the [TP2_setup repository](https://github.com/nansencenter/TP2_setup).

These notebooks read files from `$WORK/<CONFIGNAME>/expt_<EXPT_ID>/data/`.

