# Running the model

The model is run using one of the job scripts in the experiment directory. Choose the one for your scenario:

| Script | Use case |
|--------|----------|
| [`srjob.sh`](#submit-a-job) | Single run segment — hindcast, forecast, or continuation spin-up |
| [`srjob_cycle.sh`](#cycled-spin-up) | Cycled spin-up — repeats a period with GLORYS boundaries |

Each script runs `expt_preprocess.sh` at the start of every segment. This script:

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
the model starts. After the model finishes, `expt_postprocess.sh` moves restart files,
daily mean archives (`archv.*`), and CICE output from `SCRATCH/` to `data/` and
`data/cice/`, and writes `log/hycom.stop` with `GOODRUN` or `BADRUN`. The `data/`
directory is then ready for analysis with the [post-processing tools](postprocessing.md).

:::{warning}
Because all required files are staged into `SCRATCH/` before each run,
`expt_preprocess.sh` must be re-run whenever the build or any input file changes.
The job scripts do this automatically; if you run it manually beforehand, comment out
the call in the job script to avoid staging twice.
:::

:::{dropdown} Run expt_preprocess.sh manually

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
../bin/expt_preprocess.sh ${START} ${END} $INITFLG
```

Use the same `START`, `END`, and `INITFLG` values as in the job script. A successful run
prints `No fatal errors. Ok to start model set up in ..`. Useful to verify that all
required files are in place before submitting.

:::

(submit-a-job)=
## Submit a job (`srjob.sh`)

Open `srjob.sh` in a text editor and update:

| Variable | Value | Description |
|----------|-------|-------------|
| `NMPI` | `"<NMPI>"` e.g. `504` | Update NMPI if needed |
| `START` | `"YYYY-MM-DDT00:00:00"` | Run start time |
| `END` | `"YYYY-MM-DDT00:00:00"` | Run end time |
| `INITFLG` | `""` or `"--init"` | `"--init"` for a cold start; `""` to restart from files |
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


:::{note}
For reference when setting `#SBATCH --time`: a 1-year TP2 run with BGC on 4 Betzy nodes (504 cores) takes approximately 5–6 hours; a 1-year TP5 run without BGC on 5 nodes (636 cores) takes approximately 6–7 hours.
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

### Restarting after a crash

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
  ../bin/expt_postprocess.sh
  ```

  This moves the restart and archive files to `data/`, where `srjob.sh` will pick
  them up. If you skip this step and resubmit directly, `expt_preprocess.sh` will
  move the files to `SCRATCH/KEEP/` instead — they are not lost, but you will need
  to copy them to `data/` manually.

In both cases, update `START` in `srjob.sh` to the date of the last restart file
written, and set `INITFLG=""`.

## Cycled spin-up (`srjob_cycle.sh`)

`srjob_cycle.sh` runs a fixed period with GLORYS boundaries repeatedly (e.g. 1993–1997
five times). At the end of each cycle, `bin/cycle_wrap.py` shifts the restart files
valid at `CYCLE_END` back to `CYCLE_START` so the next cycle starts seamlessly. For the
initial files and nesting file requirements, see
[Cycled spin-up](forcing.md#cycled-spin-up-initial-files).

Open `srjob_cycle.sh` and update:

| Variable | Example | Description |
|----------|---------|-------------|
| `SPINUP_START` | `"1993-09-01T00:00:00"` | Start of the first cycle; may lie inside the period |
| `CYCLE_START` | `"1993-01-01T00:00:00"` | Start of cycles 2, 3, … |
| `CYCLE_END` | `"1998-01-01T00:00:00"` | End of every cycle (wrap date, exclusive) |
| `NCYCLES` | `5` | Number of cycles |
| `CYCLES_PER_JOB` | `1` | Number of cycles run in one job before it resubmits itself (see [Restarting after a crash](#restarting-after-a-crash)) |

:::{note}
A full cycle must fit within the wall-time limit (`#SBATCH --time`); raise `CYCLES_PER_JOB` if several fit. On Betzy, TP2 with BGC takes 5–6 h/year and TP5 without BGC 6–7 h/year, so a 5-year cycle takes about 25–35 h (see the timing note in [Submit a job](#submit-a-job)).
:::

Before submitting, check:

- `D` in `EXPT.src` already contains `${SPINUP_CYCLE:+/cycle_${SPINUP_CYCLE}}` in the
  template experiments, so each cycle writes to its own `data/cycle_NN/` directory.
  If you set up your experiment from scratch, verify this is in place.
- Both models must write a restart at the end of every cycle: `rstrfq` in
  `blkdat.input` must be positive and `dump_last = .true.` in `ice_in`.

Then submit:

```bash
cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
sbatch srjob_cycle.sh
```

### Restarting after a crash

`srjob_cycle.sh` always continues from the latest HYCOM/CICE restart pair in the cycle
data directories, so in most cases you can simply resubmit.

The resubmitted job finishes the interrupted cycle from that restart. The partial cycle
counts towards `CYCLES_PER_JOB`: with `CYCLES_PER_JOB=1` the job resubmits itself once
the interrupted cycle is done, and the next cycle starts in a new job. With
`CYCLES_PER_JOB=2` it finishes the interrupted cycle and then runs one more full cycle.
That is always shorter than two full cycles, so it fits the wall-time limit if two full
cycles do.

- **Model exits with an error but the SLURM job completes normally**: `expt_postprocess.sh`
  still runs and moves files to `data/cycle_NN/`. Resubmit directly.

- **Job is killed by SLURM** (wall-time limit, out-of-memory, node failure, or
  `scancel`): `expt_postprocess.sh` never runs. Output files remain in `SCRATCH/`.
  Run postprocessing manually with the cycle number of the interrupted cycle before
  resubmitting:

  ```bash
  CONFIGNAME=<CONFIGNAME>   # e.g. TP2a0.10
  EXPT_ID=<EXPT_ID>         # e.g. 01.0
  cd $WORK/${CONFIGNAME}/expt_${EXPT_ID}
  SPINUP_CYCLE=<NN> ../bin/expt_postprocess.sh
  ```

  This moves the restart and archive files to `data/cycle_NN/`, where
  `srjob_cycle.sh` will pick them up on the next run.

## Check results

When the job finishes, restart files and daily mean files are moved to `data/`.
Confirm successful completion:

```bash
cat $WORK/<CONFIGNAME>/expt_<EXPT_ID>/log/hycom.stop
```

The file should contain `GOODRUN`. If not, inspect the log files under `log/` for
error messages.

## Visualization

If you prefer to work with xarray for analysis and visualization, check out [xhycom](https://xhycom.readthedocs.io/en/latest/index.html).

Otherwise, sample Jupyter notebooks for plotting HYCOM output are also available in the [TP2_setup repository](https://github.com/nansencenter/TP2_setup).

These notebooks read files from `$WORK/<CONFIGNAME>/expt_<EXPT_ID>/data/`.

