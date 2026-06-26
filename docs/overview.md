## Model components

This setup couples several models and libraries to provide ocean, sea ice, and biogeochemical simulations.

1. **[NERSC-HYCOM-CICE](https://github.com/nansencenter/NERSC-HYCOM-CICE/tree/develop)** (develop branch)  
   Coupled ocean–sea ice–biogeochemistry framework developed at NERSC.

2. **[ECOSMO_operational / nersc](https://github.com/nansencenter/nersc)** (main branch)  
   Biogeochemical module optimized for operational use at NERSC. Cloned as the `nersc` package
   and used within the FABM framework.

3. **[GOTM](https://github.com/gotm-model/code)**  
   1D column ocean model. Only the light turbulence module is used, installed as part of FABM.

4. **[FABM](https://github.com/fabm-model/fabm)**  
   Fortran 2003 coupler for biogeochemical models. Provides the framework that links ECOSMO,
   ERSEM, and GOTM to HYCOM-CICE.

5. **[ERSEM](https://github.com/pmlmodelling/ersem)**  
   Marine biogeochemical and ecosystem model. Only the carbon chemistry module is used in
   this installation.

6. **NERSC Python library**  
   Included in the NERSC-HYCOM-CICE repository under `pythonlibs/`. Used for pre- and
   post-processing.

## Directory structure

The key convention is a separation between source and run directories: `HOME` holds
cloned repositories and compiled libraries, while `WORK` holds configuration files
and model output. This separation also matters for data safety: on many HPC systems the
work filesystem is subject to automatic purging, so source code, compiled binaries, and
configuration files should always live on home.

The same `$HOME` and `$WORK` variables are used throughout this documentation. `$HOME` is
always set for you. `$WORK` is set for you **on Betzy**, but **not on Olivia** — there you set
it yourself. Add this to your `~/.bashrc` on Olivia so that every `$WORK` reference in this
documentation resolves the same way on both machines:

```bash
# Olivia only — Betzy already sets $WORK for you
export WORK=/cluster/projects/<PROJECT>
```

With `$WORK` defined, a single directory tree applies on both machines. The filesystem each
role maps to is summarised in the table below, with per-machine specifics in the storage
dropdowns.

```
${HOME}/${USER}/
├── NERSC-HYCOM-CICE/
│   ├── pythonlibs/                  # NERSC Python libraries
│   ├── bin/                         # NERSC utility scripts
│   └── hycom/
│       ├── RELO/src_2.2.98.../      # HYCOM source files
│       └── MSCPROGS/
│           ├── bin/
│           └── lib/libhycnersc.a
├── ersem/                           # ERSEM source (carbon chemistry module)
├── fabm/                            # FABM source and compiled libraries
└── nersc/                           # NERSC-specific ECOSMO model source

${WORK}/${USER}/
└── <CONFIGNAME>/                    # e.g. TP2a0.10
    ├── bin -> (symlink to repo bin)
    ├── REGION.src
    ├── expt_<EXPT_ID>/              # e.g. expt_04.2
    │   ├── EXPT.src
    │   ├── SCRATCH/
    │   ├── data/
    │   └── build/src_.../hycom_cice
    ├── nest/<IEXPT>/                # e.g. nest/042
    ├── relax/<IEXPT>/
    └── force/synoptic/<IEXPT>/
```
:::{note}
This shows the finished directory structure for reference. The following sections build it up step by step.
:::

`$HOME` and `$USER` are set automatically on every machine; `$WORK` is set automatically on
Betzy and via your `~/.bashrc` on Olivia (above). `<CONFIGNAME>`, `<EXPT_ID>`, and `<IEXPT>`
are placeholders for user-defined values, with examples in the table below.

| Variable | Description | Example (Betzy) | Example (Olivia) |
|---|---|---|---|
| `HOME` | Root of the source/home filesystem | `/cluster/home` | `/cluster/home` |
| `WORK` | Root of the run/config tree | `/cluster/work/users` (system-set) | `/cluster/projects/<PROJECT>` (set in `~/.bashrc`) |
| `USER` | Your username on the system | `nlo043` | `nlo043` |
| `<PROJECT>` | Compute/storage allocation code | `nn9481k` | `nn9481k` |
| `<CONFIGNAME>` | Grid/configuration name | `TP2a0.10` | `TP2a0.10` |
| `<EXPT_ID>` | Experiment identifier (dot notation) | `04.2` | `04.2` |
| `<IEXPT>` | Integer experiment identifier (used in directory names) | `042` | `042` |

:::{note}
On Olivia, only the large `SCRATCH/` and `data/` directories are redirected off `$WORK` onto
the purged work filesystem (and archived to NIRD); everything else in the tree is identical.
See [Olivia storage areas](#olivia-storage-layout) below.
:::

::::{dropdown} Betzy storage areas (NRIS/Sigma2)

Betzy's storage areas (see the
[Sigma2 storage documentation](https://documentation.sigma2.no/files_storage/clusters.html)
for the authoritative description):

| Area | Path | Capacity | Backup | Purged | Role |
|---|---|---|---|---|---|
| Home | `/cluster/home/$USER` (`$HOME`) | 20 GiB / 100k files | yes (snapshots) | no | Dotfiles and small personal configs only |
| Projects | `/cluster/projects/<PROJECT>` | 1 TiB default (up to 20 TiB) | no | no | Persistent shared project data |
| User work | `/cluster/work/users/$USER` (`$WORK`) | no fixed quota | no | **yes** — a weekly scan removes files older than 21 days (up to 42 if space allows) | The run tree: configuration, `build/`, `SCRATCH/`, and `data/` output |
| Job scratch | `/cluster/work/jobs/$SLURM_JOB_ID` (`$SCRATCH`) | — | no | deleted when the job finishes | Per-job temporary space |

The cloned source and compiled libraries live under `$HOME`; the run tree lives under `$WORK`
(`/cluster/work/users/$USER`). Because `$WORK` is purged after 21 days, archive output you want
to keep to your project area or NIRD.

::::

<a name="olivia-storage-layout"></a>
::::{dropdown} Olivia storage areas (NRIS/Sigma2)

Olivia has **no per-user work area**, so `$WORK` is set to the project area
(`/cluster/projects/<PROJECT>`, configured in your `~/.bashrc` above) and the large run
directories are redirected onto the work filesystem. See the
[Sigma2 storage documentation](https://documentation.sigma2.no/files_storage/clusters.html)
for the authoritative description.

`projects` and `work` are **different Lustre filesystems** with different sizes, striping,
and purge policies — the distinction drives the layout:

| Area | Path | Capacity (shared by the allocation) | Striping | Backup | Purged | What it is for |
|---|---|---|---|---|---|---|
| Home | `/cluster/home/$USER` | ~20 GiB | — | yes | no | Dotfiles and small personal configs only — too small for the environment or build artifacts |
| Projects | `/cluster/projects/<PROJECT>` | 10 TiB hard cap | `stripe_count=1` (single OST) | no | **no** | **Small + permanent**: source, compiled binaries, container environment, experiment **configuration** and `build/` — *not* large output |
| Work | `/cluster/work/projects/<PROJECT>` | ~1.1 PB (hundreds of TiB free) | `stripe_count=4` | no | **yes, 21 days** | **Large + transient**: `SCRATCH/` and freshly-written `data/` output |
| NIRD | `/nird/datalake/<NS-PROJECT>`, `/nird/datapeak/<NS-PROJECT>` | large archive tier | — | — | no | **Large + keep long-term**: archived model output. Read-write on service (SVC) nodes, **read-only** on compute nodes, **absent** on login nodes |

**No per-user work area.** Unlike Betzy's `/cluster/work/users/$USER`, Olivia provides only the
per-project `/cluster/work/projects/<PROJECT>`, shared by everyone on the allocation. That is
why `$WORK` points at the non-purged project area and only the large directories are redirected
onto the work filesystem.

**Why the bulk run I/O (`SCRATCH/` and `data/`) goes on `work`, not `projects`.** The **10 TiB
cap is shared across the whole allocation**, so a churning `SCRATCH/` plus large multi-year
output can fill it and block *every collaborator* from writing — whereas `work` is a petascale
scratch with hundreds of TiB free. And `projects` defaults to `stripe_count=1` (writes
serialise), while `work` stripes 4-wide — what HYCOM's large `.a` archive and restart writes
want. This is why NRIS labels `work` "job data" and `projects` "user software".

**Redirecting the large directories.** With `$WORK=/cluster/projects/<PROJECT>`, the
`<CONFIGNAME>/` tree — configuration and `build/` — already lives on the non-purged filesystem.
Redirect the two large directories onto `work` by overriding `S=` (scratch) and `D=` (data) in
`EXPT.src` (see [Configure EXPT.src](experiment-setup.md#configure-exptsrc)). Because `work` is
purged by file age — so output from early in a long run is at risk *while the run is still
going* — archive completed `data/` to NIRD on a **rolling** basis (e.g. per model year as it
finishes) from a service node, not just at the end.

This redirect is safe because `expt_preprocess.sh` only creates (`mkdir -p`) and enters (`cd`)
the `$S` and `$D` directories — it never deletes them — so the overridden `S=`/`D=` paths are
honoured transparently. The purge then only ever touches `work`; the configuration and `build/`
on `projects` are untouched, and the output is yours to move to NIRD before it ages out.

:::{tip}
NIRD is not visible from the login nodes. To stage collaborator data from NIRD, copy it
on a service (SVC) node — where NIRD is read-write — or read it directly inside a batch
job, where compute nodes see it read-only. See
[Post-processing collaborator runs](mscprogs.md#post-processing-runs-from-nird) for the
recommended staging workflow.
:::

::::

