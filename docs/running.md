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
  `BADRUN` depending on whether the model reached a normal stop.

Open `srjob.sh` in a text editor and update:

| Variable | Value | Description |
|----------|-------|-------------|
| `START` | `"YYYY-MM-DDT00:00:00"` | Run start time |
| `END` | `"YYYY-MM-DDT00:00:00"` | Run end time |
| `INITFLG` | `""` or `"--init"` | `""` for a restart run; `"--init"` to initialize from climatology |
| `#SBATCH --time` | `"HH:MM:SS"` | Wall-clock time limit |

> **`INITFLG="--init"` (climatological initialization):** No restart files are needed.
> Temperature and salinity (T/S) are set directly from the climatological fields in
> `relax/` (see [Climatologies and river forcing](forcing.md#climatologies-and-river-forcing));
> velocities and sea surface height (SSH) start at zero. The model then spins up under realistic atmospheric
> forcing. The start date must be in September (the month of Arctic sea ice minimum).
> See [Initial conditions](forcing.md#initial-conditions) for details on both options.

> **Note:** As a reference, a 7-day TP2 run on 4 Betzy nodes (504 cores) takes approximately 10 minutes.

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

## Visualization

Sample Jupyter notebooks for plotting HYCOM output are available in the
[TP2_setup repository](https://github.com/nansencenter/TP2_setup):

- [plot_2D_TP2_temp.ipynb](https://github.com/nansencenter/TP2_setup/blob/main/plot_2D_TP2_temp.ipynb) — 2D surface maps
- [plot_section_TP2.ipynb](https://github.com/nansencenter/TP2_setup/blob/main/plot_section_TP2.ipynb) — vertical cross-sections

These notebooks read files from `$WORK/<CONFIGNAME>/expt_<EXPT_ID>/data/`.
