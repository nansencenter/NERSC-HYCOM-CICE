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
| `INITFLG` | `""` or `"--init"` | `"--init"` for a cold start; `""` to restart from files (see note below) |
| `#SBATCH --time` | `"HH:MM:SS"` | Wall-clock time limit |

> **`INITFLG="--init"` (cold start):** No restart files are needed.
> Temperature and salinity (T/S) are set directly from the climatological fields in
> `relax/` (see [Climatologies and river forcing](forcing.md#climatologies-and-river-forcing));
> velocities and sea surface height (SSH) start at zero. The model then spins up under realistic atmospheric
> forcing. The start date must be in September (the month of Arctic sea ice minimum).
> See [Initial conditions](forcing.md#initial-conditions) for details on both options.

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

## Visualization

If you prefer to work with xarray for analysis and visualization, check out [xhycom](https://xhycom.readthedocs.io/en/latest/index.html).

Otherwise, sample Jupyter notebooks for plotting HYCOM output are also available in the [TP2_setup repository](https://github.com/nansencenter/TP2_setup).

These notebooks read files from `$WORK/<CONFIGNAME>/expt_<EXPT_ID>/data/`.

