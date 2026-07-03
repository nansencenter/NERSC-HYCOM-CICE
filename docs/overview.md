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
cloned repositories and compiled libraries, while `WORK` holds configuration files,
compiled binaries, and model output. `$WORK` is always set to the **non-purged** projects
filesystem so that configuration files and compiled binaries survive automatic purge cycles.
The large, I/O-heavy directories (`SCRATCH/`, `data/`, `nest/`, `relax/`, `force/`) are
redirected onto the fast scratch filesystem instead: `projects` is a shared allocation with
a limited quota, and filling it with run data blocks all collaborators from writing (see the
storage dropdowns below for machine-specific paths and capacities).

Add this to your `~/.bashrc` on **all machines**:

```bash
export WORK=/cluster/projects/nn2993k/$USER
```

On Betzy, the system also sets `$WORK` to the purged user-work area — overriding it in
`~/.bashrc` points to the non-purged projects filesystem instead, keeping configuration and
compiled binaries safe. On Olivia, there is no system-set `$WORK`, so this line is required.

With `$WORK` defined, a single directory tree and a single set of setup steps apply on both
machines. The filesystem each role maps to is summarised in the table below, with per-machine
specifics in the storage dropdowns.

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

${WORK}/
└── <CONFIGNAME>/                    # e.g. TP2a0.10
    ├── bin -> (symlink to repo bin)
    ├── REGION.src
    ├── expt_<EXPT_ID>/              # e.g. expt_04.2
    │   ├── EXPT.src
    │   ├── SCRATCH/                 # on scratch via S= in EXPT.src
    │   ├── data/                    # on scratch via D= in EXPT.src
    │   └── build/src_.../hycom_cice
    ├── nest/<IEXPT>/                # e.g. nest/042   symlink → scratch
    ├── relax/<IEXPT>/               #                 symlink → scratch
    └── force/synoptic/<IEXPT>/      #                 symlink → scratch
```
:::{note}
This shows the finished directory structure for reference. The following sections build it up step by step.
:::

`$HOME` and `$USER` are set automatically on every machine; `$WORK` is set via your
`~/.bashrc` on all machines (above). `<CONFIGNAME>`, `<EXPT_ID>`, and `<IEXPT>` are
placeholders for user-defined values, with examples in the table below.

| Variable | Description | Value |
|---|---|---|
| `HOME` | Home filesystem (source code) | `/cluster/home` |
| `WORK` | Non-purged config tree (set in `~/.bashrc`) | `/cluster/projects/nn2993k/$USER` |
| `USER` | Your username | `nlo043` |
| `<CONFIGNAME>` | Grid/configuration name | `TP2a0.10` |
| `<EXPT_ID>` | Experiment identifier (dot notation) | `04.2` |
| `<IEXPT>` | Integer experiment identifier (used in directory names) | `042` |

:::{note}
On all machines, the large run directories — `SCRATCH/`, `data/`, `nest/`, `relax/`, and
`force/` — are redirected off `$WORK` onto the fast, purged scratch filesystem. `SCRATCH/`
and `data/` are redirected by overriding `S=` and `D=` in `EXPT.src`; `nest/`, `relax/`,
and `force/` by placing symlinks in the `<CONFIGNAME>/` directory. Everything else —
configuration files and `build/` — stays on `$WORK`. Archive completed `data/` to NIRD on
a rolling basis; `nest/`, `relax/`, and `force/` must be staged back from NIRD if lost to a
purge. See the storage dropdowns below for machine-specific paths.
:::

::::{dropdown} Betzy storage areas (NRIS/Sigma2)

Betzy's storage areas (see the
[Sigma2 storage documentation](https://documentation.sigma2.no/files_storage/clusters.html)
for the authoritative description):

| Area | Path | Capacity | Backup | Purged | Role |
|---|---|---|---|---|---|
| Home | `/cluster/home/$USER` (`$HOME`) | 20 GiB / 100k files | yes (snapshots) | no | Dotfiles and small personal configs only |
| Projects | `/cluster/projects/nn2993k` | 1 TiB default (up to 20 TiB) | no | no | **Small + permanent**: experiment configuration and `build/` (`$WORK`) |
| User work | `/cluster/work/users/$USER` (`$USERWORK`) | no fixed quota | no | **yes** — files older than 21 days (up to 42 if space allows) | **Large + transient**: `SCRATCH/`, `data/`, `nest/`, `relax/`, `force/` |
| Job scratch | `/cluster/work/jobs/$SLURM_JOB_ID` (`$SCRATCH`) | — | no | deleted when the job finishes | Per-job temporary space |
| NIRD | `/nird/datalake/<NS-PROJECT>`, `/nird/datapeak/<NS-PROJECT>` | large archive tier | — | no | **Large + keep long-term**: archived model output. Mounted on login nodes (read-write); not available on compute nodes |

::::


<a name="olivia-storage-layout"></a>
::::{dropdown} Olivia storage areas (NRIS/Sigma2)

Olivia's storage areas (see the
[Sigma2 storage documentation](https://documentation.sigma2.no/files_storage/clusters.html)
for the authoritative description). Unlike Betzy, Olivia has **no per-user work area** — the
work area is per-project and shared.

`projects` and `work` are **different Lustre filesystems** with different sizes, striping,
and purge policies — the distinction drives the layout:

| Area | Path | Capacity (shared by the allocation) | Striping | Backup | Purged | What it is for |
|---|---|---|---|---|---|---|
| Home | `/cluster/home/$USER` | ~20 GiB | — | yes | no | Dotfiles and small personal configs only — too small for the environment or build artifacts |
| Projects | `/cluster/projects/nn2993k` | 10 TiB hard cap | `stripe_count=1` (single OST) | no | **no** | **Small + permanent**: experiment configuration and `build/` (`$WORK`) |
| Work | `/cluster/work/projects/nn2993k` | ~1.1 PB (hundreds of TiB free) | `stripe_count=4` | no | **yes, 21 days** | **Large + transient**: `SCRATCH/`, `data/`, `nest/`, `relax/`, `force/` |
| NIRD | `/nird/datalake/<NS-PROJECT>`, `/nird/datapeak/<NS-PROJECT>` | large archive tier | — | — | no | **Large + keep long-term**: archived model output. Read-write on service (SVC) nodes, **read-only** on compute nodes, **absent** on login nodes |

On Olivia, `projects` defaults to `stripe_count=1` (serialised writes on a single OST),
while `work` stripes 4-wide — better suited to HYCOM's large `.a` archive and restart writes.

:::{tip}
NIRD is not visible from the login nodes. To stage collaborator data from NIRD, copy it
on a service (SVC) node — where NIRD is read-write — or read it directly inside a batch
job, where compute nodes see it read-only. See
[Post-processing collaborator runs](mscprogs.md#post-processing-runs-from-nird) for the
recommended staging workflow.
:::

::::

