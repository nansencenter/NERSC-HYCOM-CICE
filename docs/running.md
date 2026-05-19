## Submit a job

The job script `srjob.sh` handles three things automatically:

- **[Atmospheric forcing generation](forcing.md#atmospheric-forcing)** — downloads and
  prepares ERA5 forcing for the run period
- **SCRATCH setup** — runs `expt_preprocess.sh` to stage files needed at runtime

  ::::{dropdown} Run manually
  ```bash
  cd $WORK/<CONFIGNAME>/expt_<EXPT_ID>
  ../bin/expt_preprocess.sh ${START} ${END} $INITFLG
  ```
  Use the same `START`, `END`, and `INITFLG` values as in `srjob.sh`, see the table below for formats.
  A successful run prints a message ending with `No fatal errors. Ok to start model set up in ..`.
  ::::

- **Job submission** — submits the model to the SLURM queue

Open `srjob.sh` in a text editor and update:

| Variable | Value | Description |
|----------|-------|-------------|
| `START` | `"YYYY-MM-DDT00:00:00"` | Run start time |
| `END` | `"YYYY-MM-DDT00:00:00"` | Run end time |
| `INITFLG` | `""` or `"--init"` | `""` for a restart run; `"--init"` to initialize from climatology |
| `#SBATCH --time` | `"HH:MM:SS"` | Wall-clock time limit |

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
